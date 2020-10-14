
subroutine compute_rotating_spotted_photosphere(time_arr, n_time_points, star_prop, spot_list, &
                                                spot_params, cth_coeffs, interpolated_filter_coeffs, wv, dlnp_ph, &
                                                flpk_ph, flnp_ph, dlnp_sp, flpk_sp, flnp_sp, &
                                                flnp_fc, flt_immaculate_ph, s_index_im, flt_bolometric, s_index, &
                                                filling_factor_sp, filling_factor_fc, spectra, n_wv, n_spots)

    implicit double precision (a-h,o-z)
    
    integer, intent(in) :: n_wv, n_spots, n_time_points
    double precision, intent(in) :: time_arr, star_prop(8), cth_coeffs, interpolated_filter_coeffs(n_wv)
    double precision, intent(in) :: wv(n_wv), flpk_ph(n_wv), dlnp_ph(n_wv,18), flnp_ph(n_wv)
    double precision, intent(in) :: flpk_sp(n_wv), dlnp_sp(n_wv,18), flnp_sp(n_wv), flt_immaculate_ph(n_wv)
    double precision, intent(in) :: flnp_fc(n_wv), s_index_im
    double precision, intent(in) :: spot_list(n_spots,5), spot_params(4)
    double precision, intent(out) :: flt_bolometric(n_time_points,2), s_index(n_time_points,2), filling_factor_sp(n_time_points,2)
    double precision, intent(out) :: filling_factor_fc(n_time_points,2), spectra(n_time_points, n_wv+1)
    dimension cth_coeffs(*), time_arr(*)
    double precision :: flt_t(n_wv), x_i, r_s, w_star, evo_rate, diff_rotation, amu, are_sp, are_fc
    !declarations for external subroutines
    double precision :: flt_convoluted(n_wv), Q, dlp_ph, dlp_sp, sflp_ph, sflp_sp
    double precision :: flt_excess_ph, s, ph_flux_sp_integrated, sp_flux_sp_integrated, proj_area_spots, proj_area_faculae
   
    pi = dacos(-1.0d0)
    to_radians = pi/180.0d0  !Conversion from degrees to radians
    
    !star parameters
    x_i = star_prop(1)       ! axis inclination
    r_s = star_prop(2)       ! star radius
    w_star = star_prop(3)    ! angular speed rotation
    !rotation drift constants [deg/day]
    B = star_prop(4)
    C = star_prop(5)
    tphot = star_prop(6)
    dtfc = star_prop(8)
    !spot parameters
    evo_rate = spot_params(1)
    diff_rotation = spot_params(2) 
    Q = spot_params(4)
    
    do 99 j = 1, n_time_points
        s = s_index_im  !immaculate S-index
        !read the time value in the input array (in days!)
        t = time_arr(j)
        spectra(j,1) = t
        !copy the immaculate photosphere vector to temporal flt_out
        do i = 1, n_wv
            flt_t(i) = flt_immaculate_ph(i)
            spectra(j,i+1) = flt_t(i)
        enddo
        !compute coordinates (ph,th) and spot radius for each spot in spot_list array for all t_j
        !structure "spot_list" = [[spot_1],[spot_2],...]
        !sum of projected area of spots&faculae (to compute filling factor)
        proj_area_spots = 0.0
        proj_area_faculae = 0.0
        
        do 98 index_spot = 1,n_spots
            !spot central coordinates
            th_0 = spot_list(index_spot,3)
            ph_0 = spot_list(index_spot,4)
            !compute the drifted angular speed for a (ph,th) position
            w_drift = w_star - diff_rotation*to_radians*( B*(dsin((90.d0-th_0)*to_radians))**2.d0 + &
                                               C*(dsin((90.d0-th_0)*to_radians))**4.d0)
            !get init time and end time of each spot
            r_spot_max = spot_list(index_spot,5)
            t_init_spot = spot_list(index_spot,1)
            t_end_spot = spot_list(index_spot,2) + t_init_spot
            t_1 = t_init_spot + r_spot_max/evo_rate
            t_2 = t_end_spot - r_spot_max/evo_rate

            !compute current spot radius
            if ((t.GT.t_init_spot).AND.(t.LT.t_1)) then
                current_spot_radius = evo_rate*(t-t_init_spot)
            else if ((t.GT.t_1).AND.(t.LT.t_2)) then
                current_spot_radius = r_spot_max
            else if ((t.GT.t_2).AND.(t.LT.t_end_spot)) then
                current_spot_radius = evo_rate*(t_end_spot-t)
            else
                current_spot_radius = 0.0
            end if

            !projection of normal-to-surface and obs direction
            amu = dsin(th_0*to_radians)*dcos(ph_0*to_radians + w_drift*t)*dsin(x_i*to_radians) + &
                    dcos(th_0*to_radians)*dcos(x_i*to_radians)
           
            !check if this pixel is visible
            if (amu.GT.0.0) then
                !LD range of the pixel
                call limb_darkening_coeff(amu,l)
                !spot surface approximation 
                !are_sp = (to_radians**2)*(pi/4.0)*(2.0*current_spot_radius)**2.0
                are_sp = 2.0d0*pi*(1.0d0-cos(current_spot_radius*to_radians))
                are_fc = Q*are_sp  !(to_radians**2)*(pi/4.0)*(Q*(ang_res**2))               
                
                ph_flux_sp_integrated = 0.0
                sp_flux_sp_integrated = 0.0
                
                proj_area_spots = proj_area_spots + are_sp*amu/pi
                proj_area_faculae = proj_area_faculae + are_fc*amu/pi
                !spectral loop; integrates flux over wavelengths
                do m = 1,n_wv
                    !compute the flux for photosphere
                    
                    dlp_ph = dlnp_ph(m,l+1)+(dlnp_ph(m,l)-dlnp_ph(m,l+1))*(amu-cth_coeffs(l+1))/ &
                                (cth_coeffs(l)-cth_coeffs(l+1))   
                    sflp_ph = flnp_ph(m)*(dlp_ph/flpk_ph(m))*are_sp*amu/(4.0d0*pi)
                    
                    ph_flux_sp_integrated = ph_flux_sp_integrated + sflp_ph 
                    !compute the flux for spots
                    dlp_sp = dlnp_sp(m,l+1)+(dlnp_sp(m,l)-dlnp_sp(m,l+1))*(amu-cth_coeffs(l+1))/ &
                                (cth_coeffs(l)-cth_coeffs(l+1))
                    sflp_sp = flnp_sp(m)*(dlp_sp/flpk_sp(m))*are_sp*amu/(4.0d0*pi)
                    sp_flux_sp_integrated = sp_flux_sp_integrated + sflp_sp 
                    !compute the flux for faculae
                    sflf = flnp_fc(m)*(dlp_ph/flpk_ph(m))*are_fc*amu/(4.0d0*pi)  !we aply same LD as ph
                    sflp = flnp_ph(m)*(dlp_ph/flpk_ph(m))*are_fc*amu/(4.0d0*pi)
                    dtfmu = 250.9d0 - 407.4d0*amu + 190.9d0*amu**2.d0
                    dflf = sflf*((tphot + dtfmu)/(tphot + dtfc))**4.d0 - sflp    !over-flux of the facula limb brightened
                    !subtract excess flux in current pixel from immaculate photosphere total flux
                    flt_excess_ph = sflp_ph - sflp_sp 
                    flt_t(m) = flt_t(m) - flt_excess_ph + dflf 
                    !save the spectrum
                    spectra(j,m+1) = flt_t(m)
                end do
               
                !subtract ph flux of photosphere and add spot + faculae contribution to the S-index
                s = s - 0.164*ph_flux_sp_integrated + 3.164*sp_flux_sp_integrated*(1.0 + Q)  
                
            endif
98      continue

    !convolve the integral flux with the specified interpolated_filter_coeffs
    call filter_convolution(flt_t, interpolated_filter_coeffs, n_wv, flt_convoluted)
    !integrate the planckian to get the bolometric flux
    call trapezoid_integration(wv, flt_convoluted, n_wv, res_out)
    !fill the output bolometric flux array
    flt_bolometric(j,1) = t
    flt_bolometric(j,2) = res_out
    s_index(j,1) = t
    s_index(j,2) = s
    filling_factor_sp(j,1) = t
    filling_factor_sp(j,2) = proj_area_spots
    filling_factor_fc(j,1) = t
    filling_factor_fc(j,2) = proj_area_faculae
    
99  continue

end subroutine compute_rotating_spotted_photosphere

subroutine compute_rotating_spotted_photosphere_planet(time_arr, n_time_points, ang_res, star_prop, spot_list, &
                                                spot_params, cth_coeffs, interpolated_filter_coeffs, wv, dlnp_ph, &
                                                flpk_ph, flnp_ph, dlnp_sp, dlnp_fc, flpk_sp, flnp_sp, flnp_fc, &
                                                flt_immaculate_ph, flt_bolometric, spectra, n_wv, n_spots)

    implicit double precision (a-h,o-z)

    !parameter (ang_res=0.25d0)           !Size of the pixels on the star (deg)
    !parameter (nph=int(360/ang_res))     !Number of pixels on longitude (360/ang_res)
    !parameter (nth=int(180/ang_res))     !Number of pixels on latitude (180/ang_res)
    !parameter (n1=nph*nth)               !Total number of pixels nph*nth
    
    integer, intent(in) :: n_wv, n_spots, n_time_points
    double precision, intent(in) :: time_arr, ang_res, star_prop(14), cth_coeffs, interpolated_filter_coeffs(n_wv)
    double precision, intent(in) :: wv(n_wv), flpk_ph(n_wv), dlnp_ph(n_wv,18), flnp_ph(n_wv)
    double precision, intent(in) :: flpk_sp(n_wv), dlnp_sp(n_wv,18), flnp_sp(n_wv), flt_immaculate_ph(n_wv)
    double precision, intent(in) :: dlnp_fc(n_wv,18), flnp_fc(n_wv)
    double precision, intent(in) :: spot_list(n_spots,5), spot_params(4)
    double precision, intent(out) :: flt_bolometric(n_time_points,2), spectra(n_time_points, n_wv+1)
    dimension cth_coeffs(*), time_arr(*)
    double precision :: flt_t(n_wv), x_i, r_s, w_star, evo_rate, diff_rotation, amu
    !declarations for external subroutines
    double precision :: flt_convoluted(n_wv), Q, dlp_ph, dlp_sp, sflp_ph, sflp_sp
    double precision :: flt_excess_ph
    double precision :: spot_coords_t(n_spots,3)
    double precision, dimension(:), allocatable :: ph, th, pht
    
    pi = dacos(-1.0d0)
    to_radians = pi/180.0d0  !Conversion from degrees to radians
    nph = int(360/ang_res)
    nth = int(180/ang_res)
    n1 = nph*nth
    
    allocate (ph(n1))
    allocate (th(n1))
    allocate (pht(n1))
    
    !star parameters
    x_i = star_prop(1)          ! axis inclination
    r_s = star_prop(2)          ! star radius
    w_star = star_prop(3)       ! angular speed rotation
    !rotation drift constants; deg/day
    B = star_prop(4)
    C = star_prop(5)
    tphot = star_prop(6)
    dtfc = star_prop(8)
    !planet parameters
    T0_planet = star_prop(9)      !planet ephermeris T0(days)
    P_planet = star_prop(10)      !period of the planet
    A_planet = star_prop(11)      !planet semiaxis in R_star
    spin_orbit = star_prop(12)
    bim = star_prop(13)           !parameter of impact
    R_planet = star_prop(14)      !planet radii
    !spot parameters
    evo_rate = spot_params(1)
    diff_rotation = spot_params(2)
    Q = spot_params(4)
 
    do 99 j = 1, n_time_points
        !read the time value in the input array (in days!)
        t = time_arr(j)
        spectra(j,1) = t
        !copy the immaculate photosphere vector to temporal flt_out
        do i = 1, n_wv
            flt_t(i) = flt_immaculate_ph(i)
            spectra(j,i+1) = flt_t(i)
        enddo
        !Position of the planet. Assume a circular orbit
        x_pla = A_planet*dsin((2.0d0*pi*(t-T0_planet)/P_planet)) !in Rstar units
        v_pla = dcos(2.0d0*pi*(t-T0_planet)/P_planet) !planet velocity to select primary transits
        !polar coordinates of the planet with respect to the center of star (rho, phi) = (rho_pl, phi_pl)
        phi_pl = datan( (x_pla*dcos(spin_orbit*to_radians))/(bim + x_pla*dsin(spin_orbit*to_radians)) )
        rho_pl = (x_pla*dcos(spin_orbit*to_radians))/dsin(phi_pl)          
!       compute array coordinates (ph,th) and spot radius for each spot in spot_list array for current time
!       spot_coords_t = (th_0(1), ph_0(1), r(1)), th_0(2), ph_0(2), r(2)), ... )     
                
        do 98 index_spot = 1,n_spots
            !spot central coordinates
            th_0 = spot_list(index_spot,3)
            ph_0 = spot_list(index_spot,4)
            !get init time and end time of each spot
            r_spot_max = spot_list(index_spot,5)
            t_init_spot = spot_list(index_spot,1)
            t_end_spot = spot_list(index_spot,2) + t_init_spot
            t_1 = t_init_spot + r_spot_max/evo_rate
            t_2 = t_end_spot - r_spot_max/evo_rate
            !compute current spot radius
            if ((t.GT.t_init_spot).AND.(t.LT.t_1)) then
                current_spot_radius = evo_rate*(t-t_init_spot)
            else if ((t.GT.t_1).AND.(t.LT.t_2)) then
                current_spot_radius = r_spot_max
            else if ((t.GT.t_2).AND.(t.LT.t_end_spot)) then
                current_spot_radius = evo_rate*(t_end_spot-t)
            else
                current_spot_radius = 0.0
            end if
            !fill the time spot array
            spot_coords_t(index_spot,1) = th_0
            spot_coords_t(index_spot,2) = ph_0
            spot_coords_t(index_spot,3) = current_spot_radius           
98      continue        
        !Division of the star in pixels
        thi = ang_res/2.0d0
        phi = ang_res/2.0d0
        thaux = ang_res/2.0d0
        phaux = -ang_res/2.0d0
        do 199 i = 1, n1
            th(i) = thaux
            ph(i) = phaux + ang_res
            phaux = ph(i)
            !jump to the next latitude step:
            if((ph(i) + ang_res).ge.360.0d0) then
                thaux = thaux + ang_res
                phaux = -ang_res/2.0d0
            endif
            !compute the drifted angular speed for a (ph,th) position
            w_drift = w_star - diff_rotation*to_radians*( B*(dsin((90.d0-th(i))*to_radians))**2.0 + &
                                            C*(dsin((90.d0-th(i))*to_radians))**4.0)
            !area of pixel element
            are = 2.0d0*ang_res*to_radians*dsin(ang_res*to_radians/2.0d0)*dsin(th(i)*to_radians)
            !projection of normal-to-surface and obs direction
            amu = dsin(th(i)*to_radians)*dcos(ph(i)*to_radians + w_drift*t)*dsin(x_i*to_radians) + &
                dcos(th(i)*to_radians)*dcos(x_i*to_radians)

            pht(i) = ph(i)*to_radians + w_drift*t
            phinc = pht(i)
            thinc = th(i)*to_radians
            adu = datan2(dsin(thinc)*dsin(phinc), dcos(x_i*to_radians)*dsin(thinc)*dcos(phinc) - &
                  dsin(x_i*to_radians)*dcos(thinc))
            adu = adu/to_radians
            adu = 180.d0 - adu
            if(adu.lt.0.0) adu = adu + 360.0d0
            if(adu.gt.360.d0) adu = adu - 360.d0
            adu = adu*to_radians
            rhoi = dsin(dacos(amu))
            !pixel-Planet distance (in Rstar units)
            dpp = dsqrt((rhoi*dcos(adu)-rho_pl*dcos(phi_pl))**2.0d0 + (rhoi*dsin(adu)-rho_pl*dsin(phi_pl))**2.0d0)
            ddif = dpp - R_planet
            dpr = 0.5d0 - 90.0d0*ddif/(dabs(ang_res*amu)) !distance-R in units of projected pixel size
            !Variable dpr is the planet coverage of the pixel (0 to 1)
            if (dpr.ge.1.0d0) then  
                dpr = 1.0d0
            elseif(dpr.le.0.0d0) then
                dpr = 0.0d0
            endif
            
            !check if this pixel is visible
            if (amu.GT.0.0) then
                if (dpr.GT.0.0.AND.v_pla.GE.0.0) then
                    call limb_darkening_coeff(amu,l)
                    
                    do m = 1,n_wv 
                        !compute the flux for photosphere cell
                        dlp_ph = dlnp_ph(m,l+1)+(dlnp_ph(m,l)-dlnp_ph(m,l+1))*(amu-cth_coeffs(l+1))/ &
                                    (cth_coeffs(l)-cth_coeffs(l+1))
                        sflp_ph = flnp_ph(m)*(dlp_ph/flpk_ph(m))*are*dpr*amu/(4.0d0*pi)
                        flt_excess_ph = sflp_ph
                        flt_t(m) = flt_t(m) - flt_excess_ph 
                        !save the spectrum
                        spectra(j,m+1) = flt_t(m)
                     enddo  
                goto 199
                endif
                !check distances between center of spot and pixel
                do 200 index_spot = 1,n_spots
                    th_0 = spot_coords_t(index_spot,1)
                    ph_0 = spot_coords_t(index_spot,2)
                    R_spot = spot_coords_t(index_spot,3)
                    R_fac = R_spot*SQRT(1.0 + Q)  !radius of faculae (if exists)
                    ! distance center spot - pixel
                    call spherical_distance(th_0, ph_0, th(i), ph(i), dist)
                    if (dist.LT.R_fac) then               
                        !LD range of the pixel
                        call limb_darkening_coeff(amu,l)
                        
                        !spectral loop
                        dflf = 0.0
                        flt_excess_ph = 0.0
                        do m = 1,n_wv
                            if ((dist.LT.R_spot).AND.((dpr.eq.0.0).OR.(dpr.eq.1.0))) then
                            !compute the flux for photosphere
                            dlp_ph = dlnp_ph(m,l+1)+(dlnp_ph(m,l)-dlnp_ph(m,l+1))*(amu-cth_coeffs(l+1))/ &
                                    (cth_coeffs(l)-cth_coeffs(l+1))
                            sflp_ph = flnp_ph(m)*(dlp_ph/flpk_ph(m))*are*amu/(4.0d0*pi)
                            !compute the flux for spots
                            dlp_sp = dlnp_sp(m,l+1)+(dlnp_sp(m,l)-dlnp_sp(m,l+1))*(amu-cth_coeffs(l+1))/ &
                                        (cth_coeffs(l)-cth_coeffs(l+1))
                            sflp_sp = flnp_sp(m)*(dlp_sp/flpk_sp(m))*are*amu/(4.0d0*pi)
                            flt_excess_ph = sflp_ph - sflp_sp
                            elseif ((dist.GT.R_spot).AND.(dist.LT.R_fac).AND.((dpr.eq.0.0).OR.(dpr.eq.1.0))) then
                                !compute the flux for faculae
                                dlp_ph = dlnp_ph(m,l+1)+(dlnp_ph(m,l)-dlnp_ph(m,l+1))*(amu-cth_coeffs(l+1))/ &
                                    (cth_coeffs(l)-cth_coeffs(l+1))
                                sflp_ph = flnp_ph(m)*(dlp_ph/flpk_ph(m))*are*amu/(4.0d0*pi)
                                dlp_fc = dlnp_fc(m,l+1)+(dlnp_fc(m,l)-dlnp_fc(m,l+1))*(amu-cth_coeffs(l+1))/ &
                                            (cth_coeffs(l)-cth_coeffs(l+1))
                                sflf = flnp_fc(m)*(dlp_ph/flpk_ph(m))*are*amu/(4.0d0*pi) !we aply same LD as ph
                                sflp = flnp_ph(m)*(dlp_ph/flpk_ph(m))*are*amu/(4.0d0*pi)
                                dtfmu = 250.9d0 - 407.4d0*amu + 190.9d0*amu**2.d0
                                dflf = sflf*((tphot + dtfmu)/(tphot + dtfc))**4.d0 - sflp   !over-flux of the facula limb brightened (Unruh et al 1999)  
                            endif
                            !subtract excess flux in current pixel from immaculate photosphere total flux 
                            flt_t(m) = flt_t(m) - flt_excess_ph + dflf
                            !save the spectrum
                            spectra(j,m+1) = flt_t(m)
                        enddo
                        endif
200             continue
             endif
                            
199      continue
     
    !convolve the integral flux with the specified interpolated_filter_coeffs
    call filter_convolution(flt_t, interpolated_filter_coeffs, n_wv, flt_convoluted)
    !integrate the planckian to get the bolometric flux
    call trapezoid_integration(wv, flt_convoluted, n_wv, res_out)
    !fill the output bolometric flux array
    flt_bolometric(j,1) = t
    flt_bolometric(j,2) = res_out

99  continue

    deallocate (ph)
    deallocate (th)
    deallocate (pht)

end subroutine compute_rotating_spotted_photosphere_planet

subroutine ccf_rotating_spotted_photosphere(time_arr, star_params, spot_params, cth_coeffs, spot_list, &
                                            dlnp_ph, flpk_ph, flnp_ph, dlnp_sp, flpk_sp, flnp_sp, &
                                            flnp_fc, rv_sym, ccf_sym_ph, ccf_sym_sp, ccf_sym_fc, ccf_im, &
                                            ccf_series, flxph, n_time_points, nccf, n_wv, n_spots)

    implicit double precision (a-h,o-z)

    integer, intent(in) :: n_time_points, nccf, n_wv, n_spots
    double precision, intent(in) :: spot_params(4), cth_coeffs(*), spot_list(n_spots,5), flxph
    double precision, intent(in) :: star_params(14), time_arr(n_time_points)
    double precision, intent(in) :: flpk_ph(n_wv), dlnp_ph(n_wv,18), flnp_ph(n_wv)
    double precision, intent(in) :: flpk_sp(n_wv), dlnp_sp(n_wv,18), flnp_sp(n_wv), flnp_fc(n_wv)
    double precision, intent(in) ::  rv_sym(nccf), ccf_sym_ph(nccf), ccf_sym_sp(nccf), ccf_sym_fc(nccf), ccf_im(nccf)
    double precision, intent(out) :: ccf_series(n_time_points, nccf)
   
    double precision :: x_i, r_s, w_star, evo_rate, diff_rotation, ang_res, amu, are_sp, are_fc, to_radians
    double precision :: Q, dlp_ph, dlp_sp, sflp_ph, sflp_sp, cxds_fc, sccfph
    double precision :: cxd(3), cxu(5), rvsh1(nccf), rvsh2(nccf), rvsh3(nccf), ccft(nccf)
    double precision :: sccf_ph_cell, sccf_sp_cell, sccf_fc_cell, w_drift, rvel, dtfmu
       
    pi = dacos(-1.0d0)
    to_radians = pi/180.0d0    !Conversion from degrees to radians
    !Spot-Photosphere bisectors (solar spectra) coefficients
    cxd(1) = -1290.79545378d0
    cxd(2) =  4529.54335802d0
    cxd(3) = -3023.08850834d0
    !rotation drift constants; deg/day
    B = star_params(4)
    C = star_params(5)
    x_i = star_params(1)       ! axis inclination
    r_s = star_params(2)       ! star radius
    w_star = star_params(3)    ! angular speed rotation (rad/day)
    tphot = star_params(6)
    dtfc = star_params(8)
    evo_rate = spot_params(1)
    diff_rotation = spot_params(2)
    Q = spot_params(4)

    do 99 j = 1,n_time_points
        !read the time value in the input array (in days!)
        t = time_arr(j)
        !refresh ccft array with the immaculate photosphere ccf vector
        do k=1,nccf
            ccft(k) = ccf_im(k)    
        enddo
        
        !compute coordinates (ph,th) and spot radius for each spot in spot_list array for all t_j
        !structure "spot_list" = [[spot_1],[spot_2],...]
        do 98 index_spot = 1,n_spots
            !spot central coordinates
            th_0 = spot_list(index_spot,3)
            ph_0 = spot_list(index_spot,4)
            !compute the drifted angular speed for a (ph,th) position
            w_drift = w_star - diff_rotation*to_radians*(B*(dsin((90.d0-th_0)*to_radians))**2.0 + &
                                                         C*(dsin((90.d0-th_0)*to_radians))**4.0)
            !get init time and end time of each spot
            r_spot_max = spot_list(index_spot,5)
            t_init_spot = spot_list(index_spot,1)
            t_end_spot = spot_list(index_spot,2) + t_init_spot
            t_1 = t_init_spot + r_spot_max/evo_rate
            t_2 = t_end_spot - r_spot_max/evo_rate

            !compute current spot radius
            if ((t.GT.t_init_spot).AND.(t.LT.t_1)) then
                current_spot_radius = evo_rate*(t-t_init_spot)
            else if ((t.GT.t_1).AND.(t.LT.t_2)) then
                current_spot_radius = r_spot_max
            else if ((t.GT.t_2).AND.(t.LT.t_end_spot)) then
                current_spot_radius = evo_rate*(t_end_spot-t)
            else
                current_spot_radius = 0.0
            end if
            !current angular resolution (spots)
            ang_res = 2.d0*current_spot_radius
            !projection of normal-to-surface and obs direction
            amu = dsin(th_0*to_radians)*dcos(ph_0*to_radians + w_drift*t)*dsin(x_i*to_radians) + &
                  dcos(th_0*to_radians)*dcos(x_i*to_radians)
           
            sccf_ph_cell = 0.0
            sccf_sp_cell = 0.0
            sccf_fc_cell = 0.0
            sccfph = 0.0
            
            !check if this spot is visible
            if (amu.GT.0.0) then
                !Radial velocity of the pixel element (km/s)
                rvel = r_s*w_drift*8.05d0*dsin(x_i*to_radians)*dsin(th_0*to_radians)*dsin(ph_0*to_radians + w_drift*t)
                !LD range of the pixel
                call limb_darkening_coeff(amu,l)
                !cifist coeffs 
                call cifist(amu, cxu)
                !area of the spot and faculae
                are_sp = (to_radians**2.0)*(pi/4.0)*ang_res**2.0
                are_fc = Q*are_sp
                
                do m = 1,n_wv
                    !compute the flux for photosphere
                    dlp_ph = dlnp_ph(m,l+1)+(dlnp_ph(m,l)-dlnp_ph(m,l+1))*(amu-cth_coeffs(l+1))/ &
                                (cth_coeffs(l)-cth_coeffs(l+1))
                    
                    sflp_ph = flnp_ph(m)*(dlp_ph/flpk_ph(m))*are_sp*amu/(4.0d0*pi)
                    
                    sccf_ph_cell = sccf_ph_cell + sflp_ph 
                    
                    !compute the flux for spots
                    dlp_sp = dlnp_sp(m,l+1)+(dlnp_sp(m,l)-dlnp_sp(m,l+1))*(amu-cth_coeffs(l+1))/ &
                                (cth_coeffs(l)-cth_coeffs(l+1))
                    sflp_sp = flnp_sp(m)*(dlp_sp/flpk_sp(m))*are_sp*amu/(4.0d0*pi)
                    
                    sccf_sp_cell = sccf_sp_cell + sflp_sp
                    
                    !compute the flux for faculae
                    sflf = flnp_fc(m)*(dlp_ph/flpk_ph(m))*are_fc*amu/(4.0d0*pi) !we aply same LD as ph
                    sflp = flnp_ph(m)*(dlp_ph/flpk_ph(m))*are_fc*amu/(4.0d0*pi)
                    sccf_fc_cell = sccf_fc_cell + sflf
                    sccfph = sccfph + sflp 
                enddo    
                    
                dtfmu = 250.9d0 - 407.4d0*amu + 190.9d0*amu**2.d0
                !dflf = sflf*((tphot + dtfmu)/(tphot + dtfc))**4.d0   !flux of the facula limb brightened
                sccf_fc_cell = sccf_fc_cell*((tphot + dtfmu)/(tphot + dtfc))**4.d0
               
                sccf_ph_cell = sccf_ph_cell/flxph
                sccf_sp_cell = sccf_sp_cell/flxph
                sccf_fc_cell = sccf_fc_cell/flxph
                sccfph = sccfph/flxph
                
                !blueshift and flux blocking on the spots    
                do m = 1,nccf
                
                    cbsu_ph = cxu(1) + cxu(2)*(1.d0-ccf_sym_ph(m)) + cxu(3)*(1.d0-ccf_sym_ph(m))**2.d0 &
                              + cxu(4)*(1.d0-ccf_sym_ph(m))**3.d0 + cxu(5)*(1.d0-ccf_sym_ph(m))**4.d0
                    
                    cbsu_sp = cxu(1) + cxu(2)*(1.d0-ccf_sym_sp(m)) + cxu(3)*(1.d0-ccf_sym_sp(m))**2.d0 &
                              + cxu(4)*(1.d0-ccf_sym_sp(m))**3.d0 + cxu(5)*(1.d0-ccf_sym_sp(m))**4.d0
                    
                    
                    cxds = cxd(1) + cxd(2)*(1.d0-ccf_sym_sp(m)) + cxd(3)*(1.d0-ccf_sym_sp(m))**2.d0 !spot blocking
                    cxds = cxds/1000.d0
                    
                    cbsuf = cxu(1) + cxu(2)*(1.d0-ccf_sym_fc(m)) + cxu(3)*(1.d0-ccf_sym_fc(m))**2.d0 + &
                            cxu(4)*(1.d0-ccf_sym_fc(m))**3.d0 + cxu(5)*(1.d0-ccf_sym_fc(m))**4.d0
                    
                    cxds_fc = cxd(1) + cxd(2)*(1.d0-ccf_sym_fc(m)) + cxd(3)*(1.d0-ccf_sym_fc(m))**2.d0
                    cxds_fc = cxds_fc/1000.d0
                    
                    rvsh1(m) = rv_sym(m) + rvel + cbsu_ph
                    
                    rvsh2(m) = rv_sym(m) + rvel + amu*cxds + cbsu_sp 
                    
                    rvsh3(m) = rv_sym(m) + rvel + amu*cxds_fc + cbsuf
                    
                enddo
                
                !subtract ccf the contribution of photosphere for current spot surface from immaculate photosphere
                do ks=1,nccf
                    do ms=1,nccf-1
                        if(rv_sym(ks).le.rvsh1(1)) then
                            ccft(ks) = ccft(ks) - ccf_sym_ph(1)*sccf_ph_cell
                            goto 22
                        endif
                        if(rv_sym(ks).ge.rvsh1(nccf)) then
                            ccft(ks) = ccft(ks) - ccf_sym_ph(nccf)*sccf_ph_cell
                            goto 22
                        endif
                        
                        if(((rv_sym(ks).ge.rvsh1(ms))).and.((rv_sym(ks).le.rvsh1(ms+1)))) then
                            ccft(ks) = ccft(ks) - ccf_sym_ph(ms)*sccf_ph_cell*( (rvsh1(ms+1) - & 
                                       rv_sym(ks))/(rvsh1(ms+1)-rvsh1(ms)) ) - & 
                                       ccf_sym_ph(ms+1)*sccf_ph_cell*((rv_sym(ks)-rvsh1(ms))/(rvsh1(ms+1)-rvsh1(ms)))
                            goto 22
                        endif
                    enddo
22                  continue
                enddo
                
               !add the spot contribution
               do ks=1,nccf
                   do ms=1,nccf-1
                       if(rv_sym(ks).le.rvsh2(1)) then
                           ccft(ks) = ccft(ks) + ccf_sym_sp(1)*sccf_sp_cell
                           goto 23
                       endif
                       if(rv_sym(ks).ge.rvsh2(nccf)) then
                           ccft(ks) = ccft(ks) + ccf_sym_sp(nccf)*sccf_sp_cell
                           goto 23
                       endif
                       if(((rv_sym(ks).ge.rvsh2(ms))).and.((rv_sym(ks).le.rvsh2(ms+1)))) then
                           ccft(ks) = ccft(ks) + ccf_sym_sp(ms)*sccf_sp_cell*(rvsh2(ms+1)- & 
                                    rv_sym(ks))/(rvsh2(ms+1)-rvsh2(ms)) + & 
                                    ccf_sym_sp(ms+1)*sccf_sp_cell*(rv_sym(ks)-rvsh2(ms))/(rvsh2(ms+1)-rvsh2(ms))
                           goto 23
                       endif
                   enddo
23                  continue
               enddo            
            
               !subtract from immaculate photosphere ccf the contribution of photosphere for current faculae surface
               do ks=1,nccf
                   do ms=1,nccf-1
                       if(rv_sym(ks).lt.rvsh1(1)) then
                           ccft(ks) = ccft(ks) - ccf_sym_ph(1)*sccfph
                           goto 26
                       endif
                       if(rv_sym(ks).gt.rvsh1(nccf)) then
                           ccft(ks) = ccft(ks) - ccf_sym_ph(nccf)*sccfph
                           goto 26
                       endif
                       if((rv_sym(ks).ge.rvsh1(ms)).and.(rv_sym(ks).lt.rvsh1(ms+1))) then
                           ccft(ks) = ccft(ks) - ccf_sym_ph(ms)*sccfph*(rvsh1(ms+1)- & 
                                    rv_sym(ks))/(rvsh1(ms+1)-rvsh1(ms)) - & 
                                    ccf_sym_ph(ms+1)*sccfph*(rv_sym(ks)-rvsh1(ms))/(rvsh1(ms+1)-rvsh1(ms))
                           goto 26
                       endif
                   enddo
26                  continue
               enddo
                
                !add the faculae contribution
               do ks=1,nccf
                   do ms=1,nccf-1
                       if(rv_sym(ks).lt.rvsh2(1)) then
                           ccft(ks) = ccft(ks) + ccf_sym_fc(1)*sccf_fc_cell
                           goto 25
                       endif
                       if(rv_sym(ks).gt.rvsh2(nccf)) then
                           ccft(ks) = ccft(ks) + ccf_sym_fc(nccf)*sccf_fc_cell
                           goto 25
                       endif
                       if((rv_sym(ks).ge.rvsh3(ms)).and.(rv_sym(ks).lt.rvsh3(ms+1))) then
                           ccft(ks) = ccft(ks) + ccf_sym_fc(ms)*sccf_fc_cell*(rvsh3(ms+1)- & 
                                    rv_sym(ks))/(rvsh3(ms+1)-rvsh3(ms)) + & 
                                    ccf_sym_fc(ms+1)*sccf_fc_cell*(rv_sym(ks)-rvsh3(ms))/(rvsh3(ms+1)-rvsh3(ms))
                           goto 25
                       endif
                   enddo
25                  continue
               enddo                           
                       
            endif
98      continue
      
      do i=1,nccf
          ccf_series(j,i) = ccft(i)
      enddo
      

99  continue


      

end subroutine ccf_rotating_spotted_photosphere

subroutine limb_darkening_coeff(amu, l)
    !returns the limb darkening coefficient, given an amu component
    double precision :: amu
    integer :: l

    if (amu.le.1.0.and.amu.gt.0.9) then
        l=1
    else if (amu.le.0.9.and.amu.gt.0.8) then
        l=2
    else if (amu.le.0.8.and.amu.gt.0.7) then
        l=3
    else if (amu.le.0.7.and.amu.gt.0.6) then
        l=4
    else if (amu.le.0.6.and.amu.gt.0.5) then
        l=5
    else if (amu.le.0.5.and.amu.gt.0.4) then
        l=6
    else if (amu.le.0.4.and.amu.gt.0.3) then
        l=7
    else if (amu.le.0.3.and.amu.gt.0.25) then
        l=8
    else if (amu.le.0.25.and.amu.gt.0.2) then
        l=9
    else if (amu.le.0.2.and.amu.gt.0.15) then
        l=10
    else if (amu.le.0.15.and.amu.gt.0.125) then
        l=11
    else if (amu.le.0.125.and.amu.gt.0.1) then
        l=12
    else if (amu.le.0.1.and.amu.gt.0.075) then
        l=13
    else if (amu.le.0.075.and.amu.gt.0.05) then
        l=14
    else if (amu.le.0.05.and.amu.gt.0.025) then
        l=15
    else if (amu.le.0.025.and.amu.gt.0.01) then
        l=16
    else if (amu.le.0.01.and.amu.ge.0.0) then
        l=17
    end if

end subroutine limb_darkening_coeff

subroutine trapezoid_integration(x, y, n_points, res_out)

    !returns the integral of y along x.
    integer, intent(in) :: n_points
    double precision, intent(in) :: x(n_points), y(n_points)
    double precision, intent(out) :: res_out
    double precision :: da

    res_out = 0.0

    do 11 i = 1, n_points-1
        da = 0.5*(x(i+1) - x(i))*(y(i+1) + y(i))
        res_out = res_out + da
11  continue
    res_out = res_out

end subroutine trapezoid_integration

subroutine rms(model, data, n_points, rms_out)

    !return root mean square between data and spot pattern
    integer, intent(in) :: n_points
    double precision, intent(in) :: model(n_points,2), data(n_points,3)
    double precision, intent(out) :: rms_out
    double precision :: summ

    summ = 0.0
    do 11 i = 1, n_points
        summ = summ + (model(i,2) - data(i,2))**2
11  continue

    rms_out = SQRT(summ/n_points)

end subroutine rms

subroutine khi2(model, data, n_points, khi2_out)

    !return root mean square between data and spot pattern
    integer, intent(in) :: n_points
    double precision, intent(in) :: model(n_points,2), data(n_points,3)
    double precision, intent(out) :: khi2_out

    khi2_out = 0.0
    do 11 i = 1, n_points
        khi2_out = khi2_out + ( (model(i,2) - data(i,2))/(data(i,3)) )**2
11  continue
        khi2_out = khi2_out/n_points

end subroutine khi2

subroutine compute_ccf_immaculate_photosphere(star_prop, cth_coeffs, flpk, dlnp, flnp, &
                                              flxph, rv_sym, ccf_sym, n_wv, nccf, rvsh, ccft)

    implicit double precision (a-h,o-z)

    parameter (ang_res=1.0d0)            !Size of the pixels on the star (deg)
    parameter (nph=int(360/ang_res))     !Number of pixels on longitude (360/ang_res)
    parameter (nth=int(180/ang_res))     !Number of pixels on latitude (180/ang_res)
    parameter (n1=nph*nth)               !Total number of pixels nph*nth
    parameter (nrp=10000)                !Maximum number of points in the spectra

    integer, intent(in) :: n_wv, nccf
    double precision, intent(in) :: star_prop(3), cth_coeffs, flpk(n_wv), dlnp(n_wv,18), flnp(n_wv), flxph
    double precision, intent(in) :: rv_sym(nccf), ccf_sym(nccf)
    double precision, intent(out) :: rvsh(nccf), ccft(nccf)
    dimension cth_coeffs(*) ! wvk(1219)
    double precision :: ph(n1), th(n1), cxu(5)
    double precision :: x_i, r_s, w, rvel, amu, cbsu, sccf=0.0d0

    character :: typ(n1)*2
    pi = dacos(-1.0d0)
    rad = pi/180.0d0                          !Conversion from degrees to radianswrite(*,*) flnp(1), flnp(2)

    x_i = star_prop(1)  ! axis inclination
    r_s = star_prop(2)  ! star radius
    w = star_prop(3)    ! angular speed rotation

    !write(*,*) "Num points:", n_wv
    do i = 1,nccf
       rvsh(i) = 0.0
       ccft(i) = 0.0
    enddo

!   Division of the star in pixels
    thi = ang_res/2.0d0
    phi = ang_res/2.0d0

    thaux = ang_res/2.0d0
    phaux = -ang_res/2.0d0

!   Loop to generate the photosphere (space)
    do 99 i = 1,n1
        th(i) = thaux
        ph(i) = phaux + ang_res
        phaux = ph(i)
!   jump to the next latitude step:
        if((ph(i) + ang_res).ge.360.0d0) then
            thaux = thaux + ang_res
            phaux = -ang_res/2.0d0
         endif
! area of pixel element:
         are = 2.0d0*ang_res*rad*dsin(ang_res*rad/2.0d0)*dsin(th(i)*rad)
! projection of normal-to-surface and obs direction:
         amu = dsin(th(i)*rad)*dcos(ph(i)*rad)*dsin(x_i*rad) + dcos(th(i)*rad)*dcos(x_i*rad)
! area of pixel element at the center of the disk:
         aredc = 2.0d0*ang_res*rad*dsin(ang_res*rad/2.0d0)*dsin(90.d0*rad)
! visibility of pixel element:
        if(amu.ge.0.0) then
            vis = 1.0
        else
            vis = 0.0
        endif
! temperature of the pixel element: all defined as photosphere (ph)
        typ(i) = 'ph'
! Radial velocity of the pixel element (km/s)
        rvel = r_s*w*8.05d0*dsin(x_i*rad)*dsin(th(i)*rad)*dsin(ph(i)*rad)
! SED of the pixel element: ph SED convolved with a LD(wv) from the Kurucz_interpol
! We generate the LD(wv) function for the pixel position
!       LD range of the pixel:
        call limb_darkening_coeff(amu,l)
        if(typ(i).eq.'ph'.and.vis.gt.0.5) then
            ! cifist bisector fit coeffs
            call cifist(amu,cxu)
            sccf = 0.0d0
            do m=1,n_wv
                dlp = dlnp(m,l+1)+(dlnp(m,l)-dlnp(m,l+1))*(amu-cth_coeffs(l+1))/(cth_coeffs(l)-cth_coeffs(l+1))
                sccf = sccf + flnp(m)*(dlp/flpk(m))*are*amu/(4.0d0*pi)
            enddo
        
        sccf = sccf/flxph

         !compute the convective blueshift
         do m=1,nccf
             cbsu = cxu(1) + cxu(2)*(1.0d0-ccf_sym(m)) + cxu(3)*(1.0d0-ccf_sym(m))**2.0 + cxu(4)*(1.0d0-ccf_sym(m))**3.0 + &
                    cxu(5)*(1.0d0-ccf_sym(m))**4.0
             rvsh(m) = rv_sym(m) + rvel + cbsu
         enddo
        !interpolate RV to find integrated CCF
        do ks = 1,nccf
            do ms = 1,nccf-1
                if (rv_sym(ks).lt.rvsh(1)) then
                    ccft(ks) = ccft(ks) + ccf_sym(1)*sccf
                    goto 19
                 endif
                 if (rv_sym(ks).gt.rvsh(nccf)) then
                    ccft(ks) = ccft(ks) + ccf_sym(nccf)*sccf
                    goto 19
                 endif
                 if ((rv_sym(ks).ge.rvsh(ms)).and.(rv_sym(ks).lt.rvsh(ms+1))) then
                     ccft(ks) = ccft(ks) + ccf_sym(ms)*sccf*(rvsh(ms+1)-rv_sym(ks))/(rvsh(ms+1)-rvsh(ms)) + &
                     ccf_sym(ms+1)*sccf*(rv_sym(ks)-rvsh(ms))/(rvsh(ms+1)-rvsh(ms))
                     goto 19
                 endif
            enddo
19      continue 
        enddo
     endif ! tancament del bucle de visibilitat
      
99   continue

end subroutine compute_ccf_immaculate_photosphere

subroutine filter_convolution(flt, interpolated_filter_coeffs, n_wv, flt_convoluted)

    !returns the integral of y along x.
    integer, intent(in) :: n_wv
    double precision, intent(in) :: flt(n_wv), interpolated_filter_coeffs(n_wv)
    double precision, intent(out) :: flt_convoluted(n_wv)

    do i = 1, n_wv
        flt_convoluted(i) = flt(i)*interpolated_filter_coeffs(i)
    end do

    return

end subroutine filter_convolution

subroutine compute_immaculate_photosphere(star_prop, cth_coeffs, flpk, dlnp, flnp, n_wv, flt, flxph, s_index_im)

    implicit double precision (a-h,o-z)

    parameter (ang_res=1.0d0)            !Size of the pixels on the star (deg)
    parameter (nph=int(360/ang_res))     !Number of pixels on longitude (360/ang_res)
    parameter (nth=int(180/ang_res))     !Number of pixels on latitude (180/ang_res)
    parameter (n1=nph*nth)               !Total number of pixels nph*nth
    parameter (nrp=10000)                !Maximum number of points in the spectra

    integer, intent(in) :: n_wv
    double precision, intent(in) :: star_prop(8), cth_coeffs, flpk(n_wv), dlnp(n_wv,18), flnp(n_wv)
    double precision, intent(out) :: flt(n_wv), flxph, s_index_im
    dimension cth_coeffs(*) ! wvk(1219)
    double precision :: ph(n1), th(n1)
    double precision :: sflp, x_i, r_s, w

    character :: typ(n1)*2
    pi = dacos(-1.0d0)
    rad = pi/180.0d0    !Conversion from degrees to radians

    x_i = star_prop(1)  ! axis inclination
    r_s = star_prop(2)  ! star radius
    w = star_prop(3)    ! angular speed rotation

    !write(*,*) "Num points:", n_wv
    do i = 1,n_wv
       flt(i) = 0.0
    enddo
    
    flxph = 0.0   !sum of integrated immaculate photosphere flux
    s_index_im = 0.0d0 
    !   Division of the star in pixels
    thi = ang_res/2.0d0
    phi = ang_res/2.0d0

    thaux = ang_res/2.0d0
    phaux = -ang_res/2.0d0

!   Loop to generate the photosphere (space)
    do 99 i = 1,n1
        th(i) = thaux
        ph(i) = phaux + ang_res
        phaux = ph(i)
!   jump to the next latitude step:
        if((ph(i) + ang_res).ge.360.0d0) then
            thaux = thaux + ang_res
            phaux = -ang_res/2.0d0
         endif
! area of pixel element:
         are = 2.0d0*ang_res*rad*dsin(ang_res*rad/2.0d0)*dsin(th(i)*rad)
! projection of normal-to-surface and obs direction:
         amu = dsin(th(i)*rad)*dcos(ph(i)*rad)*dsin(x_i*rad) + dcos(th(i)*rad)*dcos(x_i*rad)

! area of pixel element at the center of the disk:
         aredc = 2.0d0*ang_res*rad*dsin(ang_res*rad/2.0d0)*dsin(90.d0*rad)
! visibility of pixel element:
        if(amu.ge.0.0) then
            vis = 1.0
        else
            vis = 0.0
        endif
! temperature of the pixel element: all defined as photosphere (ph)
        typ(i) = 'ph'
! SED of the pixel element: ph SED convolved with a LD(wv) from the Kurucz_interpol
! We generate the LD(wv) function for the pixel position

!       LD range of the pixel:
        call limb_darkening_coeff(amu,l)
       
        if(typ(i).eq.'ph'.and.vis.gt.0.5) then
            sccf = 0.0d0
            do m=1,n_wv
                dlp = dlnp(m,l+1)+(dlnp(m,l)-dlnp(m,l+1))*(amu-cth_coeffs(l+1))/(cth_coeffs(l)-cth_coeffs(l+1))
                sflp = flnp(m)*(dlp/flpk(m))*are*amu/(4.0d0*pi)
                flt(m) = flt(m) + sflp
                sccf = sccf + sflp
            enddo
            flxph = flxph + sccf
            s_index_im = s_index_im + 0.164*sccf
        endif

99   continue

end subroutine compute_immaculate_photosphere

subroutine limb_darkening_to_btsettl(typ, nspres, wvmin, nwv, acd, wv, dlnp, flpk, flnp)
!Generates NextGen SEDs with interpolated LD factors from Kurucz_interpol
!45  format('./x99-corr-maskfit spec_ref.dat 'a11' 'f5.2)
!Declaration instructions
    implicit double precision (a-h,o-z)
    parameter (nrp=500000)  !Maximum number of points in the spectra
    integer, intent(in) :: nspres, nwv
    character*2, intent(in) :: typ
    double precision, intent(in) :: wvmin, acd
    double precision, intent(out) :: wv(nwv), dlnp(nwv,18), flpk(nwv), flnp(nwv)
    dimension acd(*), dlkp(1219,18), wvk(1219), flkp(1219)
    !dimension wv(nrp)

    if (typ.EQ.'ph') then
        open(1,file='./data/specph.tmp')
        open(2,file='./data/specph.dat')
        open(3,file='./data/spkurph.dat',status='unknown')
    else if (typ.EQ.'sp') then
        open(1,file='./data/specsp.tmp')
        open(2,file='./data/specsp.dat')
        open(3,file='./data/spkursp.dat',status='unknown')
    else if (typ.EQ.'fc') then
        open(1,file='./data/specfc.tmp')
        open(2,file='./data/specfc.dat')
        open(3,file='./data/spkurfc.dat',status='unknown')
    else
        write(*,*) "Type not known. Possible values: ph, sp, fc"
    endif

    !reads the LD factors from generated Kurucz_interpol
    read(3,*)
    do j=1,1219
       read(3,*) wvk(j),flkp(j),(dlkp(j,l), l=2,17)
    enddo
    close(3)

    do i=1,nspres
        read(1,*) wva
        if(wva.ge.wvmin) then
            do m = 1,nwv
                read(1,*) wv(m),flnp(m)
                do n = 1,1219
                    if(wv(m).gt.wvk(n).and.wv(m).le.wvk(n+1)) then
                        do l = 2,17
                            dlnp(m,l) = dlkp(n,l) + (dlkp(n+1,l) - dlkp(n,l))*(wv(m) - wvk(n))/(wvk(n+1) - wvk(n))
                        enddo
                    endif
                enddo

        dlnp(m,1) = 100000.0d0
        dlnp(m,18) = dlnp(m,17) - 0.67*(dlnp(m,16) - dlnp(m,17))
        flpk(m) = 0.0d0
        do j = 1,17
            flpk(m) = flpk(m) + ((dlnp(m,j) + dlnp(m,j+1))*0.5d0)*(sin(acd(j+1))**2.0d0 - sin(acd(j))**2.0d0)
            enddo
        write(2,*) wv(m),flnp(m),flpk(m),(dlnp(m,l), l=1,18)
            enddo
    goto 23
        endif
    enddo
23  continue
    close(2)
    close(1)

end subroutine limb_darkening_to_btsettl

subroutine limb_darkening_to_btsettl_fast(typ, n_points, wvmin, wvmax, acd, wv, dlnp, flpk, flnp)
!Generates NextGen SEDs with interpolated LD factors from Kurucz_interpol
!45  format('./x99-corr-maskfit spec_ref.dat 'a11' 'f5.2)
!Declaration instructions
    implicit double precision (a-h,o-z)
    parameter (nrp=500000)  !Maximum number of points in the spectra
    integer, intent(in) :: n_points
    character*2, intent(in) :: typ
    double precision, intent(in) :: wvmin, wvmax, acd
    double precision, intent(out) :: wv(n_points), dlnp(n_points,18), flpk(n_points), flnp(n_points)
    dimension acd(*), dlkp(1219,18), wvk(1219), flkp(1219)
    !dimension wv(nrp)
    
    if (typ.EQ.'ph') then
        open(1,file='./data/specph.tmp')
        open(2,file='./data/specph.dat')
        open(3,file='./data/spkurph.dat',status='unknown')
    
    else if (typ.EQ.'sp') then
        open(1,file='./data/specsp.tmp')
        open(2,file='./data/specsp.dat')
        open(3,file='./data/spkursp.dat',status='unknown')
    
    else if (typ.EQ.'fc') then
        open(1,file='./data/specfc.tmp')
        open(2,file='./data/specfc.dat')
        open(3,file='./data/spkurfc.dat',status='unknown')
    else
        write(*,*) "Type not known. Possible values: ph, sp, fc"
    endif

    !reads the LD factors from generated Kurucz_interpol
    read(3,*)
    do j=1,1219
       read(3,*) wvk(j),flkp(j),(dlkp(j,l), l=2,17)
    enddo
    close(3)
   
    do i=1,n_points
        read(1,*) wva
        if ((wva.ge.wvmin).AND.(wva.lt.wvmax)) then 
            do m = 1, int(n_points-i)
                read(1,*) wv(m), flnp(m)
                if (wv(m).GE.wvmax) then
                    goto 23
                endif
                do n = 1,1219
                    if ((wv(m).gt.wvk(n)).and.(wv(m).le.wvk(n+1))) then
                        do l = 2,17
                            dlnp(m,l) = dlkp(n,l) + (dlkp(n+1,l) - dlkp(n,l))*(wv(m) - wvk(n))/(wvk(n+1) - wvk(n))
                        enddo
                    endif
                enddo

        dlnp(m,1) = 100000.0d0
        dlnp(m,18) = dlnp(m,17) - 0.67*(dlnp(m,16) - dlnp(m,17))
        flpk(m) = 0.0d0
        do j = 1,17
            flpk(m) = flpk(m) + ((dlnp(m,j) + dlnp(m,j+1))*0.5d0)*(sin(acd(j+1))**2.0d0 - sin(acd(j))**2.0d0)
            enddo
        write(2,*) wv(m),flnp(m),flpk(m),(dlnp(m,l), l=1,18)
            enddo
    goto 23
        endif
    enddo
23  continue
  
    close(2)
    close(1)

end subroutine limb_darkening_to_btsettl_fast

subroutine kurucz_interpol(teff, grv, nme, spkz_out)

! ROUTINE TO INTERPOLATE KURUCZ ATMOSPHERE MODELS
! Model file: ip00k0new.pck
! Input: effective temperature and gravity:
!    1. Linear interpolation on gravity
!    2. Liner interpolation on temperature
! Output: atmosphere model, wavelength, central intensity, coefficients

!     Declaration instructions
      implicit double precision (a-h,o-z)

      parameter (nl=1221) !Wavelength values
      parameter (nm=17)   !Angle coefficients
      parameter (ni=5)    !4 models to interpolate + 1 model interpolated
      character model*29,modelo*29 ! row*115,
      dimension wv(nl),fl(nl,ni),ifl(nl,nm,ni),ang(nm),gmod(11),tmod(20)
      dimension model(ni),iok(ni)
      dimension xfl1(nm),xfl2(nm)

      double precision, intent(in) :: teff, grv
      double precision, intent(out) :: spkz_out(1221,18)
      character*2, intent(in) :: nme

!     Format instructions
1000  format('TEFF',f8.0,'  GRAVITY',f8.5)
1001  format(f9.2,1x,ES9.3,16(i6))

      ! in case of lower Teff's take the lowest possible temperature in Kurucz grid (only for LD)
      if (teff.LT.3500.0) then
          teff_c = 3500.0
      else
          teff_c = teff
      endif
       
!     Define the Teff models in the Kurucz file
      data tmod / 3500., 3750., 4000., 4250., 4500., 4750., 5000.,5250., 5500., &
                 5750., 6000., 6250., 6750., 7000.,7250., 7500., 8000., 8250., 8500., 8750. /
      data gmod / 0.0, 0.5, 1.0, 1.5, 2.0, 2.5,3.0, 3.5, 4.0, 4.5, 5.0 /
      data ang / 1.000, 0.900, 0.800, 0.700, 0.600, 0.500, 0.400,0.300, 0.250,  &
                0.200, 0.150, 0.125, 0.100, 0.075,0.050, 0.025, 0.010 /
!
!     write(*,*)'OBTAINING LIMB DARKENING COEFICIENTS'
      if(nme.eq.'ph') then
!       write(*,*) ' Photosphere'
      elseif(nme.eq.'sp') then
!       write(*,*) ' Spots'
      endif
      open(unit=1,file='./data/ip00k0new.pck',status='old',blank='zero')

!     Selection of models to interpolate according to Teff
      do 10 i=1,19
         if((teff_c.ge.tmod(i)).and.(teff_c.lt.tmod(i+1))) then
            teff1=tmod(i)
            teff2=tmod(i+1)
         endif
10    continue
      if(teff_c.lt.tmod(1)) then
         teff1=tmod(1)
         teff2=tmod(2)
         write(*,*) '   Extrapolation of models due to low Teff!!'
      endif
      if(teff_c.ge.tmod(20)) then
         teff1=tmod(19)
         teff2=tmod(20)
         write(*,*) '   Extrapolation of models due to large Teff!!'
      endif
!
!     Selection of models to interpolate according to log(g)
      do 20 i=1,10
         if((grv.ge.gmod(i)).and.(grv.lt.gmod(i+1))) then
            grv1=gmod(i)
            grv2=gmod(i+1)
         endif
20    continue
      if(grv.lt.gmod(1)) then
         grv1=gmod(1)
         grv2=gmod(2)
         write(*,*) '   Extrapolation of models due to low log(g)!!'
      endif
      if(grv.ge.gmod(11)) then
         grv1=gmod(10)
         grv2=gmod(11)
         write(*,*) '   Extrapolation of models due to large log(g)!!'
      endif
      write(model(1),1000) teff1,grv1
      write(model(2),1000) teff1,grv2
      write(model(3),1000) teff2,grv1
      write(model(4),1000) teff2,grv2

!
!     Read the corresponding models to interpolate
!     write(*,*) ' '
!     write(*,*) ' Reading the Kurucz models...'
      ioktot=0
      do k=1,4
!         write(*,*) '   ',model(k)
          iok(k)=0
      enddo
25    read(1,'(a29)',end=30) modelo
         do 40 k=1,4
            if(modelo.eq.model(k)) then
               iok(k)=1
               read(1,*)
               read(1,*)
               do 50 i=1,nl
                  read(1,1001) wv(i),fl(i,k),(ifl(i,j,k), j=2,nm)
50             continue
               ioktot=ioktot+iok(k)
               if(k.eq.4) go to 30 !!This way the code does not read till the end of the file
               go to 25
            endif
40       continue
         go to 25
30    continue
      close(unit=1)
      if(ioktot.eq.0) then
         do 60 i=1,4
            if(iok(i).eq.0) write(*,*) '   Model ',model(i),' not available!!'
60       continue
         go to 999
      endif
!      do k=1,4
!      write(2,*) model(k)
!          do i=1,nl
!             write(2,1001) wv(i),fl(i,k),(ifl(i,j,k), j=2,nm)
!          enddo
!      enddo
!
!     Interpolation of models
!     write(*,*)
!     write(*,*) ' Interpolation according to log(g) and Teff...'
!     Interpolation of first model teff1,grv
!     To round of angle values, add 0.5
      do 70 i=1,nl
!        Models Teff1, grv
         a0=fl(i,1)
         pend=(fl(i,2)-fl(i,1))/(grv2-grv1)
         fl1=a0+pend*(grv-grv1)
         do 80 j=1,nm
            a0=dfloat(ifl(i,j,1))
            pend=(dfloat(ifl(i,j,2)-ifl(i,j,1)))/(grv2-grv1)
            xfl1(j)=a0+pend*(grv-grv1)
80       continue
!        Model Teff2, grv
         a0=fl(i,3)
         pend=(fl(i,4)-fl(i,3))/(grv2-grv1)
         fl2=a0+pend*(grv-grv1)
         do 90 j=1,nm
            a0=dfloat(ifl(i,j,3))
            pend=(dfloat(ifl(i,j,4)-ifl(i,j,3)))/(grv2-grv1)
            xfl2(j)=a0+pend*(grv-grv1)
90       continue
!        Model Teff, grv
         a0=fl1
         pend=(fl2-fl1)/(teff2-teff1)
         fl(i,5)=a0+pend*(teff_c-teff1)
         do 100 j=1,nm
            a0=xfl1(j)
            pend=(xfl2(j)-xfl1(j))/(teff2-teff1)
            fl_aux=a0+pend*(teff_c-teff1)
            ifl(i,j,5)=int(fl_aux)
            fl_aux=fl_aux-dfloat(ifl(i,j,5))
            if(fl_aux.gt.0.5) ifl(i,j,5)=ifl(i,j,5)+1
100       continue
70    continue

      !write the output array
      do 120 i=1,nl
          spkz_out(i,1) = wv(i)
          spkz_out(i,2) = fl(i,5)
          do 121 j=3,nm+1
               spkz_out(i,j) = ifl(i,j-1,5)
121       continue
120   continue

!     Write output files
       if(nme.eq.'ph') then
        open(2,file='./data/spkurph.dat',status='unknown')
       else if(nme.eq.'sp') then
        open(2,file='./data/spkursp.dat',status='unknown')
       else if(nme.eq.'fc') then
        open(2,file='./data/spkurfc.dat',status='unknown')
       endif
      !open(unit=2,file='kurucz_model.tmp',status='unknown',blank='zero')
       write(2,'(a7,f8.0,a9,f8.5)') '# TEFF ',teff,' GRAVITY ',grv
       do 110 i=1,nl   !nl = 1221
          write(2,1001) wv(i),fl(i,5),(ifl(i,j,5), j=2,nm)   !nm=17
110    continue
       close(unit=2)
!        write(*,*)' '
!        write(*,*)' End'
999    continue

end subroutine kurucz_interpol

subroutine hcorr(nm, maskin, ccfrng, stp, nrv, ref_file, rv, ccf)
!     calculation of CCF (Cross-correlation Function) for photospere, spots or faculae, using a mask
!     and high-resolution Phoenix HR spectra.
!     Input: nm: type of surface ('ph', 'sp', 'fc')
!            maskin: mask (spectral type)
!            nrv: number of points in output CCF
!     Output: rv and ccf: arrays of the cross-correlation function       
      implicit double precision (a-h,o-z)
      !parameter (stp=0.25d0) !radial velocity step
      
      character*2, intent(in) :: nm, maskin, ref_file
      double precision, intent(in) :: ccfrng, stp
      integer, intent(in) :: nrv
      double precision, intent(out) :: rv(nrv), ccf(nrv) 
      
      dimension wv(2000000), wv1(20000), wv2(20000), al(20000), flx(2000000)
      dimension step(2000000), stepo(2000000)
      character*26 hmask
      character*40 str
      double precision :: ccfnor = 0.0
      integer :: nwv=0
      
      !number of CCF points
!     nrv = 2*int(ccfrng/stp)     
!     open the files to write the normalized spectra       
      if(nm.eq.'ph') then
          str = './data/tmp/phoenixHR_ph_norm.tmp'//ref_file
          open(1,file=str,status='unknown')
      elseif(nm.eq.'sp') then
          str = './data/tmp/phoenixHR_sp_norm.tmp'//ref_file
          open(1,file=str,status='unknown')
      elseif(nm.eq.'fc') then
          str = './data/tmp/phoenixHR_fc_norm.tmp'//ref_file
          open(1,file=str,status='unknown')
      endif
!     read the mask file, and count the lines (nwv).            
      i = 1      
      do while(1.EQ.1)
          read(1,*,end=34) wv(i),flx(i)
          if(i.gt.1) step(i-1)=wv(i)-wv(i-1)
          if(i.gt.1) stepo(i-1)=step(i-1)
          nwv = nwv + 1 
          i = i + 1
      enddo
34    continue
      close(1)
      call sort (nwv-1,stepo)
95    format('./data/mask'a2'.mas')
      write(hmask,95) maskin
      open(1,file=hmask, status='old')
      n=1
1     read(1,*,end=2) w1,w2,al(n)
          if (maskin.EQ.'IR') then
              wm = w1
          else 
              wm=(w1+w2)/2.d0
          endif
      wamp=wm*820.d0/3.d8
      wv1(n)=(wm-wamp/2.d0)/10.d0
      wv2(n)=(wm+wamp/2.d0)/10.d0
          
      n=n+1
      goto 1
2     nlin=n-1
      close(1)
      do k=1,nrv
          rv(k)=dfloat(k)*stp-ccfrng
          ccf(k)=0.d0
      do i=1,nlin
          wl1=wv1(i)*(1.d0+rv(k)/3.d5)
          wl2=wv2(i)*(1.d0+rv(k)/3.d5)
          jin=int((wl1-stepo(nwv-1)-wv(1))/stepo(nwv-1))

       If(jin.lt.1) jin=1
       do j=jin,nwv-1
           if(wv(j).gt.wl2) goto 3
           if((wv(j)-(wl1-step(j)))*(wv(j)-wl2).lt.0.d0) then
           if((wv(j)-(wl1-step(j)))*(wv(j)-wl1).lt.0.d0) then
               ccf(k)=ccf(k)+flx(j)*al(i)*(wv(j+1)-wl1)/step(j)
           elseif((wv(j)-wl1)*(wv(j)-(wl2-step(j))).lt.0d0) then
               ccf(k)=ccf(k)+flx(j)*al(i)
           else
               ccf(k)=ccf(k)+flx(j)*al(i)*(wl2-wv(j))/step(j)
           endif
           endif
       enddo
3      continue
       enddo
          
       if(k.eq.1) ccfnor=ccf(k)
       ccf(k) = ccf(k)/ccfnor
       enddo

       if(nm.eq.'ph') then
           open(1,file='./data/corr-specphn.dat',status='unknown')
       elseif(nm.eq.'sp') then
           open(1,file='./data/corr-specspn.dat',status='unknown')
       elseif(nm.eq.'fc') then
           open(1,file='./data/corr-specfcn.dat',status='unknown')
       endif
     
       do i=1,nrv
           ccf(i) = 1.d0 - ccf(i)
           write(1,*) rv(i), ccf(i)
       enddo
       close(1)
       
end subroutine hcorr

subroutine sort(n, arr)
!     From Numerical Recipes Software &H1216.
      INTEGER n,M,NSTACK
      DOUBLE PRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=l-1
2         arr(i+1)=a
12      continue
        if(jstack.eq.0) return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK) write(*,*) 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      
end subroutine sort

subroutine cifist(amu, cxu)

      implicit double precision (a-h,o-z)
      dimension cxu(5),amv(10),cxm(10,5)

      amv(1) = 1.d0
      amv(2) = 0.9d0
      amv(3) = 0.8d0
      amv(4) = 0.7d0
      amv(4) = 0.6d0
      amv(4) = 0.5d0
      amv(4) = 0.4d0
      amv(4) = 0.3d0
      amv(4) = 0.2d0
      amv(4) = 0.1d0
 
! mu1
 
      cxm(1,1) = 6.26663288564337506d-2
      cxm(1,2) = -2.0670933431270311d0  
      cxm(1,3) = 4.8373645518385970d0  
      cxm(1,4) = -6.4516296955583172d0    
      cxm(1,5) = 3.2655301356607209d0
 
! mu2
 
      cxm(2,1) = -0.69431946130967059d0     
      cxm(2,2) = 2.4039636042874215d0 
      cxm(2,3) = -5.2709508982761601d0    
      cxm(2,4) = 3.9784546433290275d0  
      cxm(2,5) = -0.77400176062747106d0
 
! mu3
 
      cxm(3,1) = -0.93259382612796937d0   
      cxm(3,2) = 3.7052371760766065d0 
      cxm(3,3) = -8.0873003920809197d0    
      cxm(3,4) = 6.8106821000360460d0  
      cxm(3,5) = -1.9069955429552707d0
 
! mu4
 
      cxm(4,1) = -1.0719225011468299d0   
      cxm(4,2) = 4.7355761410870096d0  
      cxm(4,3) = -10.463490540258894d0 
      cxm(4,4) = 9.2119708126641431d0  
      cxm(4,5) = -2.8367029084003463d0
 
! mu5
 
      cxm(5,1) = -1.9847283247315139d0   
      cxm(5,2) = 10.373894688794481d0  
      cxm(5,3) = -22.838311272740228d0    
      cxm(5,4) = 21.099573534520708d0  
      cxm(5,5) = -7.1080317369820998d0
 
! mu6
 
      cxm(6,1) = -2.7222970291487565d0   
      cxm(6,2) = 15.147888265002793d0  
      cxm(6,3) = -33.419644715367561d0    
      cxm(6,4) = 31.339604644844020d0  
      cxm(6,5) = -10.818449896888437d0
 
! mu7
 
      cxm(7,1) = -3.8435431779824496d0   
      cxm(7,2) = 22.381494580870587d0  
      cxm(7,3) = -49.497213554810209d0    
      cxm(7,4) = 46.974167758989005d0  
      cxm(7,5) = -16.493747061693675d0
 
! mu8
 
      cxm(8,1) = -5.3053173888148306d0   
      cxm(8,2) = 31.798214713677620d0  
      cxm(8,3) = -70.258482051995955d0    
      cxm(8,4) = 67.049595174910948d0  
      cxm(8,5) = -23.742151830133839d0
 
! mu9
 
      cxm(9,1) = -8.0208534457096885d0   
      cxm(9,2) = 48.285687363847394d0  
      cxm(9,3) = -105.34777535451389d0    
      cxm(9,4) = 100.04480989529668d0  
      cxm(9,5) = -35.341406918690780d0
 
! mu10
 
      cxm(10,1) = -13.178925404208321d0    
      cxm(10,2) = 76.993453025476114d0  
      cxm(10,3) = -163.43475518573922d0    
      cxm(10,4) = 152.54368824442119d0  
      cxm(10,5) = -53.071375270718036d0
 

      do k=1,9
          if(amu.le.amv(k).and.amu.gt.amv(k+1)) then
              do i=1,5
                  cxu(i)=cxm(k,i)+((amu-amv(k+1))/(amv(k+1)-amv(k)))*(cxm(k+1,i)-cxm(k,i))
              enddo
          endif
      enddo
      
      if(amu.le.amv(10)) then
          do i=1,5
              cxu(i)=cxm(10,i)+((amu-amv(10))/(amv(10)-amv(9)))*(cxm(10,i)-cxm(9,i))
          enddo
      endif
 
      
end subroutine cifist

subroutine cifist2(amu, cxu)

      implicit double precision (a-h,o-z)
      double precision cxu(5),amv(10),cxm(10,5)

      amv(1)=1.d0
      amv(2)=0.9d0
      amv(3)=0.8d0
      amv(4)=0.7d0
      amv(4)=0.6d0
      amv(4)=0.5d0
      amv(4)=0.4d0
      amv(4)=0.3d0
      amv(4)=0.2d0
      amv(4)=0.1d0

! mu1

      cxm(1,1)=6.26663288564337506d-2
      cxm(1,2)=-2.0670933431270311d0
      cxm(1,3)=4.8373645518385970d0
      cxm(1,4)=-6.4516296955583172d0
      cxm(1,5)=3.2655301356607209d0

! mu2

      cxm(2,1)=-0.69431946130967059d0
      cxm(2,2)=2.4039636042874215d0
      cxm(2,3)=-5.2709508982761601d0
      cxm(2,4)=3.9784546433290275d0
      cxm(2,5)=-0.77400176062747106d0

! mu3

      cxm(3,1)=-0.93259382612796937d0
      cxm(3,2)=3.7052371760766065d0
      cxm(3,3)=-8.0873003920809197d0
      cxm(3,4)=6.8106821000360460d0
      cxm(3,5)=-1.9069955429552707d0

! mu4

      cxm(4,1)=-1.0719225011468299d0
      cxm(4,2)=4.7355761410870096d0
      cxm(4,3)=-10.463490540258894d0
      cxm(4,4)=9.2119708126641431d0
      cxm(4,5)=-2.8367029084003463d0

! mu5

      cxm(5,1)=-1.9847283247315139d0
      cxm(5,2)=10.373894688794481d0
      cxm(5,3)=-22.838311272740228d0
      cxm(5,4)=21.099573534520708d0
      cxm(5,5)=-7.1080317369820998d0

! mu6

      cxm(6,1)=-2.7222970291487565d0
      cxm(6,2)=15.147888265002793d0
      cxm(6,3)=-33.419644715367561d0
      cxm(6,4)=31.339604644844020d0
      cxm(6,5)=-10.818449896888437d0

! mu7

      cxm(7,1)=-3.8435431779824496d0
      cxm(7,2)=22.381494580870587d0
      cxm(7,3)=-49.497213554810209d0
      cxm(7,4)=46.974167758989005d0
      cxm(7,5)=-16.493747061693675d0

! mu8

      cxm(8,1)=-5.3053173888148306d0
      cxm(8,2)=31.798214713677620d0
      cxm(8,3)=-70.258482051995955d0
      cxm(8,4)=67.049595174910948d0
      cxm(8,5)=-23.742151830133839d0

! mu9

      cxm(9,1)=-8.0208534457096885d0
      cxm(9,2)=48.285687363847394d0
      cxm(9,3)=-105.34777535451389d0
      cxm(9,4)=100.04480989529668d0
      cxm(9,5)=-35.341406918690780d0

! mu10

      cxm(10,1)=-13.178925404208321d0
      cxm(10,2)=76.993453025476114d0
      cxm(10,3)=-163.43475518573922d0
      cxm(10,4)=152.54368824442119d0
      cxm(10,5)=-53.071375270718036d0
      
      do k=1,9
        if(amu.le.amv(k).and.amu.gt.amv(k+1)) then
        do i=1,5
            cxu(i)=cxm(k,i)+((amu-amv(k+1))/(amv(k+1)-amv(k)))*(cxm(k+1,i)-cxm(k,i))
        enddo
        endif
      enddo
      
      if(amu.le.amv(10)) then
      do i=1,5
        cxu(i)=cxm(10,i)+((amu-amv(10))/(amv(10)-amv(9)))*(cxm(10,i)-cxm(9,i))
      enddo
      endif

      
end subroutine cifist2

subroutine spherical_distance(th_0, ph_0, th, ph, dist)
     
     double precision :: th_0, ph_0, th, ph, dist, pi, to_radians, cos_D
     
     pi = dacos(-1.0d0)
     to_radians = pi/180.0d0  !Conversion from degrees to radians

     cos_D = dsin(pi/2 - to_radians*th_0)*dsin(pi/2 - to_radians*th) + &
             dcos(pi/2 - to_radians*th_0)*dcos(pi/2 - to_radians*th)*dcos(to_radians*abs(ph_0 - ph))
             
     dist = dacos(cos_D)/to_radians
     
end subroutine spherical_distance

subroutine check_spot_config(spot_map, spot_params, n_spots, bool_out)

    implicit double precision (a-h,o-z)
    
    !integer, intent(in) :: n_spots
    double precision, intent(in) :: spot_map(n_spots,5), spot_params(4)
    logical, intent(out) :: bool_out
    
    bool_out = .True.
    evo_rate = spot_params(1)
    diff_rotation = spot_params(2)
    Q = spot_params(4)

    do i=1,n_spots
        do j=i+1,n_spots
            t_ini_0 = spot_map(i,1)
            t_ini = spot_map(j,1)
            t_life_0 = spot_map(i,2)
            t_life = spot_map(j,2)
            r_0 = spot_map(i,5)
            r = spot_map(j,5)
            th_0 = spot_map(i,3)
            ph_0 = spot_map(i,4)
            th = spot_map(j,3)
            ph = spot_map(j,4)
            
            call spherical_distance(th_0,ph_0,th,ph,dist)
            
            if ((dist.LT.(sqrt(Q+1)*(r_0 + r))).AND.(.not.((t_ini.GT.(t_ini_0 + &
                    t_life_0).OR.(t_ini + t_life).LT.t_ini_0)))) then
                    
                    bool_out = .False.
            endif
        enddo
    enddo
            
end subroutine check_spot_config

subroutine check_spot_config_w(spot_map, spot_params, n_spots, bool_out, sp1, sp2)

    implicit double precision (a-h,o-z)
    
    !integer, intent(in) :: n_spots
    double precision, intent(in) :: spot_map(n_spots,5), spot_params(4)
    logical, intent(out) :: bool_out
    integer, intent(out) :: sp1, sp2
    
    sp1 = 0
    sp2 = 0
    bool_out = .True.
    evo_rate = spot_params(1)
    diff_rotation = spot_params(2)
    Q = spot_params(4)

    do i=1,n_spots
         
        do j=i+1,n_spots
            t_ini_0 = spot_map(i,1)
            t_ini = spot_map(j,1)
            t_life_0 = spot_map(i,2)
            t_life = spot_map(j,2)
            r_0 = spot_map(i,5)
            r = spot_map(j,5)
            th_0 = spot_map(i,3)
            ph_0 = spot_map(i,4)
            th = spot_map(j,3)
            ph = spot_map(j,4)
            
            call spherical_distance(th_0,ph_0,th,ph,dist)
            
            if ((dist.LT.(sqrt(Q+1)*(r_0 + r))).AND.(.not.((t_ini.GT.(t_ini_0 + &
                    t_life_0).OR.(t_ini + t_life).LT.t_ini_0)))) then
                    sp1 = i
                    sp2 = j
                    bool_out = .False.
            endif
        enddo
    enddo
            
end subroutine check_spot_config_w

subroutine y_gauss(x, params, n_points, y)

    integer, intent(in) :: n_points
    double precision, intent(in) :: x(n_points), params(4)
    double precision, intent(out) :: y(n_points)
   
    do 11 i = 1, n_points
        y(i) = params(1)*exp(-((x(i)-params(2))**2)/(2.0*(params(3)**2))) + params(4)
11  continue


end subroutine y_gauss





