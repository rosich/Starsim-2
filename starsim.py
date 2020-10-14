#!/usr/bin/python
#-----------------------------------------------------------------------------------------------------------------------------------
#Name:      starsim_2.py
#
#Purpose:   Python implementation of E.Herrero's starsim_v11 photosphere simulator to allow parallelism and inverse problem.
#           StarSim-2 provides simulation of spots and faculae effects on both photometric and RV time series due to star rotation.
#
#Author:    Albert Rosich (rosich@ice.cat).
#
#First Version:   2014/05/04
#
#Last Update: 2020/05/31

#-----------------------------------------------------------------------------------------------------------------------------------
"""
The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""
#-----------------------------------------------------------------------------------------------------------------------------------

import getopt
import copy
import os
import numpy as np
import sys
import glob
import time
from math import sin, cos, acos, sqrt, pi, log, log10, exp, asin
import random
import multiprocessing as mp
import bz2
sys.path.append('./bin')
import mstarsim
sys.path.append('./src')
import starsim_generic
from EnvGlobals import Glob
from CStarsim import Star, Planet, StarsimSimulation, Spectra, Starsim_IO 
from CStarsim import Spots, CPhotometry, CRadialVelocity, CBIS, CFWHM, CCONTRAST

       
if __name__ == '__main__':
    """
    """
#---------------------------------------------------------------------------#
# COMMON PART: ALL VARIABLES DECLARED HERE ARE ACCESSIBLE IN ALL SCOPES     #
#---------------------------------------------------------------------------#
    
    io = Starsim_IO() #makes available IO functions
    #option capture
    options, remainder = getopt.gnu_getopt(sys.argv[1:], 'c:i:m:s:l:xcfv', ['planet','ncpus=','inhibit_msg','mode=','spotmap=',
                                                                            'label=', 'indices', 'bw', 'save_ccfs', 'save_spectra'])
    Glob.ncpus = mp.cpu_count()
    mode = None
    num_spots = 0
    label = ''
  
    #argument parsing
    for opt, arg in options:
        """
        """
        if opt in ('-n', '--ncpus'):
            Glob.ncpus = int(arg)
        elif opt in ('-i', '--inhibit_msg'):
            Glob.inhibit_msg = True
        elif opt in ('-m', '--mode'):
            mode = str(arg)
        elif opt in ('-s', '--spotmap'):
            num_spots = int(arg)
        elif opt in ('-p', '--planet'):
            Glob.planet = True
        elif opt in ('-f', '--save_ccfs'):
            Glob.save_ccfs = True
        elif opt in ('-x', '--indices'):    
            Glob.indices = True
        elif opt in ('-c', '--bw'):
            Glob.bw = True
        elif opt in ('-l', '--label'):
            label = str(arg)
        elif opt in ('-v', '--save_spectra'):
            Glob.save_spectra = True

    
    # print some init info
    io.print_start_info()
#---------------------------------------------------------------------------------------------------------------#        
#                                            SELECT OPTION                                                      #
#                                            -------------                                                      #                                                                                                             #   
#---------------------------------------------------------------------------------------------------------------#

    if mode == 'ph':
        """compute a single band lightcurve. If file passed, its time-vector is captured to generate the curve,
           otherwise an equallyspaced time array is generated using t1, t2 and t_exp params of starsim.conf file.
        """ 
        #generate a simulation object
        Simulation = StarsimSimulation()
        
        OB = []
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
        except:
            print("Could not be read the multiobject configuration or no data selected")
        
        if len(OB) > 0 and not Glob.planet:
            for ph in OB:
                #is photometry object?
                if ph.data_type == 'ph':     
                    #ph.set_ph_zero_point(ph.z_ph)
                    #write lightcurve file
                    logL_ph = ph.compute_lc_series(write_lc=True)
                    try:
                        print("\nStatistics:")
                        print("\tNull-model logL:", ph.logL_null) 
                        print("\tModel logL:", logL_ph[0])
                        #print("\tDlogL(model-data):", logL_ph - ph.logL_null)
                        #print("")
                        #generate data curve with z_ph implemented + jitter
                        ph.write_lc_data()
                        #write data series + residuals (last column)
                        ph.compute_lc_series(write_lc=True)
                    except:
                        pass
                    ph.obs_time = np.linspace(ph.dt1, ph.dt2, num=int(1440.0*(ph.dt2-ph.dt1)/ph.t_exp))
                    ph.obs_data = []
                    
                    #compute the HR lc curve
                    ph.compute_lc_series(write_lc=True, mark='HR')
                  
        
        else:
            #no data files in multifit.conf
            print("No photometry objects in ./multifit.conf")
            Ph = CPhotometry(Simulation)
        
            if not Glob.planet:
                #compute lc
                Ph.compute_lc_series(write_lc=True)
                
            else:
                #compute lc
                #Ph.compute_lc_series(write_lc=True)
                
                #transit duration (days)
                T_dur = (Simulation.P_planet/np.pi)*asin(sqrt((Simulation.R_star+Simulation.R_planet)**2 - Simulation.bim**2)/Simulation.A_planet)
                Ph.obs_time = []
                N_transits = (Simulation.dt2 - Simulation.dt1) / Simulation.P_planet
                t_init_transit = (int(((Ph.dt1 - Simulation.T0_planet)/ Simulation.P_planet)) + 1)*Simulation.P_planet + Simulation.T0_planet 
                
                for i in range(int(N_transits)):
                     Ph.obs_time.append(t_init_transit + i*Simulation.P_planet)
                print("Number of centers:", len(Ph.obs_time))
                #compute lc centers
                Ph.compute_lc_series(write_lc=True, mark='centers')
                
                #compute all the transit (remember time_cadence parameter should be enough short!)
                for i in range(int(N_transits)):
                    Ph.obs_time.extend(np.linspace(t_init_transit-2.*T_dur+i*Simulation.P_planet, 
                                       t_init_transit+2.*T_dur+i*Simulation.P_planet, num=int(1440.0*(T_dur)/Ph.t_exp)))
                
                #compute lc
                Ph.compute_lc_series(write_lc=True)
                
    elif mode == 'ph_grid':
        """computes photometry (full series) by integrating user-defined size grid (StarSim-1)
        """
        
        Simulation = StarsimSimulation()
        #photometry object
        Ph = CPhotometry(Simulation)
        
        Glob.planet = True
        Simulation.R_planet = 0.0
        #compute lc
        Ph.compute_lc_series(write_lc=True)
        
    elif mode == 'transit_centers':
        """
        """
        #generate a simulation object
        Simulation = StarsimSimulation()

        Ph = CPhotometry(Simulation)

        #Ph.set_wv_range(1881.0, 1975.0)     
        Glob.planet = True
  
        T_dur = (Simulation.P_planet/np.pi)*asin(sqrt((Simulation.R_star+Simulation.R_planet)**2 - Simulation.bim**2)/Simulation.A_planet)
        Ph.obs_time = []
        N_transits = (Simulation.dt2 - Simulation.dt1) / Simulation.P_planet
        t_init_transit = (int(((Ph.dt1 - Simulation.T0_planet)/ Simulation.P_planet)) + 1)*Simulation.P_planet + Simulation.T0_planet
                
        for i in range(int(N_transits)):
            Ph.obs_time.append(t_init_transit + i*Simulation.P_planet)
        print("Number of centers:", len(Ph.obs_time))
        #compute lc centers
        Ph.compute_lc_series(write_lc=True, mark='centers')
        
        #compute all the transit (remember time_cadence parameter should be enough short!)
        for i in range(int(N_transits)):
            Ph.obs_time.extend(np.linspace(t_init_transit-2.*T_dur+i*Simulation.P_planet,
                                            t_init_transit+2.*T_dur+i*Simulation.P_planet, num=int(1440.0*(T_dur)/Ph.t_exp)))
                
        #compute lc
        Ph.compute_lc_series(write_lc=True, mark='transits')

        #compute activity signal without any planet
        Glob.planet = False
        Ph.compute_lc_series(write_lc=True) 
        
    elif mode == 'transit_centers_mod':
        """
        """
        def find_nearest(array, value):
            idx = abs(np.abs(array)-value).argmin()
            return array[idx]
        
        lock = mp.Lock()
        
        N_iterations = 10
        
        for iter in range(N_iterations):
            # generate a simulation object
            Simulation = StarsimSimulation()
            
            Ph = CPhotometry(Simulation)
            
            spot_temp = random.gauss(500.0, 100.0) 
            # take the nearest temperature. Carefully set the maximum range allowed, ie "... in range(70)"
            delta_t_sp = find_nearest([100.0 + 20.0*j for j in range(70)], spot_temp)   
            # load spectra...
            lock.acquire()
            Simulation.delta_t_sp = delta_t_sp
            Ph.set_delta_t()
            #Simulation.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=Simulation.delta_t_fc)
            lock.release()
            
            print("sampled T_sp:", spot_temp, "approximated by", delta_t_sp)
            
            
            Glob.planet = True
    
            T_dur = (Simulation.P_planet/np.pi)*asin(sqrt((Simulation.R_star+Simulation.R_planet)**2 - Simulation.bim**2)/Simulation.A_planet)
            Ph.obs_time = []
            N_transits = (Simulation.dt2 - Simulation.dt1) / Simulation.P_planet
            t_init_transit = (int(((Ph.dt1 - Simulation.T0_planet)/ Simulation.P_planet)) + 1)*Simulation.P_planet + Simulation.T0_planet
                    
            for i in range(int(N_transits)):
                Ph.obs_time.append(t_init_transit + i*Simulation.P_planet)
            print("Number of centers:", len(Ph.obs_time))
            #compute lc centers
            Ph.compute_lc_series(write_lc=True, mark='centers_' + str(delta_t_sp))
            
            #compute all the transit (remember time_cadence parameter should be enough short!)
            for i in range(int(N_transits)):
                Ph.obs_time.extend(np.linspace(t_init_transit-2.*T_dur+i*Simulation.P_planet,
                                                t_init_transit+2.*T_dur+i*Simulation.P_planet, num=int(1440.0*(T_dur)/Ph.t_exp)))
                    
            #compute lc
            Ph.compute_lc_series(write_lc=True, mark='transits_' + str(delta_t_sp))

            #compute activity signal without any planet
            Glob.planet = False
            Ph.compute_lc_series(write_lc=True) 
            
            # delete objects defined at the begining
            Simulation.__del__
            Ph.__del__
    
    elif mode == 'transit_centers2':
        """
        """
        def single_sim(band_lo, band_hi, tspot):
          """
          """
          lock = mp.Lock()
          #generate a simulation object
          Simulation = StarsimSimulation()
          Ph = CPhotometry(Simulation)
         
          delta_t_sp = max(np.round(tspot /20)*20, 20.0)

          lock.acquire()
          Simulation.delta_t_sp = delta_t_sp
          Ph.set_delta_t()
          #Simulation.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=Simulation.delta_t_fc)
          lock.release()
         
          Glob.planet = True

          LB, UB = band_lo, band_hi  
          Ph.set_wv_range(LB,UB)  
     
          print Simulation.R_star

          T_dur = (Simulation.P_planet/np.pi)*asin(sqrt((Simulation.R_star+Simulation.R_planet)**2 - Simulation.bim**2)/Simulation.A_planet)
         
          Ph.obs_time = []
          N_transits = (Simulation.dt2 - Simulation.dt1) / Simulation.P_planet
          t_init_transit = (int(((Ph.dt1 - Simulation.T0_planet)/ Simulation.P_planet)) + 1)*Simulation.P_planet + Simulation.T0_planet
             
          for i in range(int(N_transits)):
                  Ph.obs_time.append(t_init_transit + i*Simulation.P_planet)
          print "Number of centers:", len(Ph.obs_time)
          #compute lc centers
          Ph.compute_lc_series(write_lc=True, mark='centers_{:0.1f}'.format(delta_t_sp))
     
          #compute all the transit (remember time_cadence parameter should be enough short!)
          for i in range(int(N_transits)):
                  Ph.obs_time.extend(np.linspace(t_init_transit-2.*T_dur+i*Simulation.P_planet,
                                              t_init_transit+2.*T_dur+i*Simulation.P_planet, num=int(1440.0*(T_dur)/Ph.t_exp)))
             
          #compute lc
          Ph.compute_lc_series(write_lc=True, mark='transits_{:0.1f}'.format(delta_t_sp))

          #compute activity signal without any planet
          Glob.planet = False

          Ph.compute_lc_series(write_lc=True)
          del(Simulation)
          del(Ph)
          #Simulation.__del__
          #Ph.__del__
       
        """
        #import copy
        #Tspot_min, Tspot_max = 200.0, 600.0
        Tspot_min, Tspot_max = 260.0, 320.0
        Glob.band_lo, Glob.band_up = 1100.0, 1155.0
        for Tspot in np.arange(Tspot_min, Tspot_max, 20.0):
            #Glob_temp = copy.deepcopy(Glob)
            #Glob = Glob_temp
            print('Processing band {:0.1f}-{:0.1f}, Tspot:{:0.1f}'.format(Glob.band_lo, Glob.band_up, Tspot))
            print (Tspot)
            single_sim(Glob.band_lo, Glob.band_up, Tspot)
        """
        single_sim(1881.0, 1975.0, 300.0)
        
    elif mode == 'parallel_transit_centers':
        """
        """
        def run(params):
             """
             """
             band_lo, band_hi, tspot = params[:]
             single_sim(band_lo, band_hi, tspot)

        def single_sim(band_lo, band_hi, tspot):
            """
            """
            lock = mp.Lock()
            
            #generate a simulation object
            Simulation = StarsimSimulation()
            Ph = CPhotometry(Simulation)

            delta_t_sp = max(np.round(tspot /20)*20, 20.0)

            lock.acquire()
            Simulation.delta_t_sp = delta_t_sp
            Ph.set_delta_t()
            #Simulation.update_spectral_data(delta_t_sp=Simulation.delta_t_sp,delta_t_fc=Simulation.delta_t_fc)
            lock.release()

            Glob.planet = True

            LB, UB = band_lo, band_hi
            Ph.set_wv_range(LB,UB)

            #print Simulation.R_star

            T_dur = (Simulation.P_planet/np.pi)*asin(sqrt((Simulation.R_star+Simulation.R_planet)**2- Simulation.bim**2)/Simulation.A_planet)

            Ph.obs_time = []
            N_transits = (Simulation.dt2 - Simulation.dt1) /Simulation.P_planet
            t_init_transit = (int(((Ph.dt1 - Simulation.T0_planet)/Simulation.P_planet)) + 1)*Simulation.P_planet + Simulation.T0_planet

            for i in range(int(N_transits)):
                    Ph.obs_time.append(t_init_transit +i*Simulation.P_planet)
            print "Number of centers:", len(Ph.obs_time)
            #compute lc centers
            Ph.compute_lc_series(write_lc=True,mark='centers_{:0.1f}'.format(delta_t_sp))

            #compute all the transit (remember time_cadence parametershould be enough short!)
            for i in range(int(N_transits)):
                Ph.obs_time.extend(np.linspace(t_init_transit-2.*T_dur+i*Simulation.P_planet,t_init_transit+2.*T_dur+i*Simulation.P_planet,
                                   num=int(1440.0*(T_dur)/Ph.t_exp)))

            #compute lc
            Ph.compute_lc_series(write_lc=True, mark='transits_{:0.1f}'.format(delta_t_sp))

            #compute activity signal without any planet
            Glob.planet = False

            Ph.compute_lc_series(write_lc=True)
            del(Simulation)
            del(Ph)
     
     
        ncpus = mp.cpu_count()
        pool = mp.Pool(ncpus)

        T_spot = np.arange(100.0, 460.0,20)
        filters = [500.0, 600.0, 810.0, 1100.0, 1155.0, 1212.75, 1273.3875, 1337.056875, 
                   1403.90971875, 1474.10520469, 1547.81046492, 1625.20098817, 1706.46103758, 
                   1791.78408946,1881.37329393,1975.44195862]
        
        params = []
        for t in T_spot:
            for i in range(len(filters)-1):
                config = [round(filters[i],1), round(filters[i+1],1), t]
                params.append(config)
        
        print params
        
        try:
            parallel_runs = pool.map_async(run, params).get(9999999)

        except KeyboardInterrupt:
            print ''
            print_message('Parallel pool killed. Bye!', index=3, color=31)
            pool.terminate()
            sys.exit()

        except:
            print "Unexpected error in pool.map_async(parallel_optimal):", sys.exc_info()[0]
            raise

        pool.terminate()
        
    elif mode == 'rv':
        """run RV forward problem for a given spotmap (./data/spot_list.dat)
        """
        #generate a simulation object
        Simulation = StarsimSimulation()
        #list of objects
        OB = []
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
   
        except:
            print("Could not be read the multiobject configuration or no data selected")
  
        if len(OB) > 0:
            for rv in OB:
                if rv.data_type == 'rv':
                    print("Going to use the multifit.conf data: ", rv.science_data) 
                    #compute rv
                    logL_rv = rv.compute_rv_series(write_rv=True)
                    print("\nStatistics:")
                    print("\tNull-model logL:", rv.logL_null) 
                    print("\tModel logL:", logL_rv)
                    print("\tDlogL(model-data):", logL_rv - rv.logL_null)
                    print("")
                    #write data series with added jitter
                    rv.write_rv_data()
                    #modify time array to generate HR curves
                    rv.obs_time = np.linspace(rv.dt1-10.0, rv.dt2+10.0, num=int(1440.0*(rv.dt2-rv.dt1)/rv.t_exp))
                    rv.obs_data = []
                    #compute the HR rv curve
                    rv.compute_rv_series(write_rv=True, mark='HR')
                    
        else:
            #no data files in multifit.conf
            RV = CRadialVelocity(Simulation)
            #compute rv
            RV.compute_rv_series(write_rv=True) 
                                                                            
    elif mode == 'bis':
        """run RV forward problem for a given spotmap (./data/spot_list.dat)
        """
        #generate a simulation object
        Simulation = StarsimSimulation()
        #list of objects
        OB = []
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
   
        except:
            print("Could not be read the multiobject configuration or no data selected")
  
        if len(OB) > 0:
            for bis in OB:
                if bis.data_type == 'bis':
                    bis.write_bis_data()
                    print("Going to use the multifit.conf data: ", bis.science_data) 
                    #compute rv
                    logL_bis = bis.compute_BIS_series(write_rv=True)
                    print("\nStatistics:")
                    print("\tNull-model logL:", bis.logL_null) 
                    print("\tModel logL:", logL_bis)
                    print("\tDlogL(model-data):", logL_bis - bis.logL_null)
                    print("")
                    #modify time array to generate HR curves
                    bis.obs_time = np.linspace(bis.dt1-10.0, bis.dt2+10.0, num=int(1440.0*(bis.dt2-bis.dt1)/bis.t_exp))
                    bis.obs_data = []
                    
                    #compute the HR rv curve
                    bis.compute_BIS_series(write_rv=True)
        else:
            #no data files in multifit.conf
            bis = CBIS(Simulation)
            #compute rv
            bis.compute_BIS_series(write_rv=True) 
        
    elif mode == 'fwhm':
        """run RV forward problem for a given spotmap (./data/spot_list.dat)
        """
        #generate a simulation object
        Simulation = StarsimSimulation()
        #list of objects
        OB = []
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
   
        except:
            print("Could not be read the multiobject configuration or no data selected")
  
        if len(OB) > 0:
            for fwhm in OB:
                if fwhm.data_type == 'fwhm':
                    print("Going to use the multifit.conf data: ", fwhm.science_data) 
                    #compute rv
                    logL_fwhm = fwhm.compute_FWHM_series(write_rv=True)
                    print("\nStatistics:")
                    print("\tNull-model logL:", fwhm.logL_null) 
                    print("\tModel logL:", logL_fwhm)
                    print("\tDlogL(model-data):", logL_fwhm - fwhm.logL_null)
                    print("")
        else:
            #no data files in multifit.conf
            fwhm = CFWHM(Simulation)
            #compute rv
            fwhm.compute_FWHM_series(write_rv=True) 
            
    elif mode == 'contrast':
        """run RV forward problem for a given spotmap (./data/spot_list.dat)
        """
        #generate a simulation object
        Simulation = StarsimSimulation()
        #list of objects
        OB = []
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
   
        except:
            print("Could not be read the multiobject configuration or no data selected")
  
        if len(OB) > 0:
            for contrast in OB:
                if contrast.data_type == 'contrast':
                    print("Going to use the multifit.conf data: ", contrast.science_data) 
                    #compute rv
                    logL_contrast = contrast.compute_CONTRAST_series(write_rv=True)
                    print("\nStatistics:")
                    print("\tNull-model logL:", contrast.logL_null) 
                    print("\tModel logL:", logL_contrast)
                    print("\tDlogL(model-data):", logL_contrast - contrast.logL_null)
                    print("")
        else:
            #no data files in multifit.conf
            contrast = CCONTRAST(Simulation)
            #compute rv
            contrast.compute_CONTRAST_series(write_rv=True) 
    
    elif mode == 'FOM':
        """
        """
        def inversion(parallel_bunch):
            """receives a list of Simulation objects, inverts observational data 
               and computes forward problem
            """        
            bunch_data = []
            for S in parallel_bunch:
                #init random number generator
                np.random.seed()
                #random jitter
                s = s_old + random.gauss(0.0,0.001)
                while not (0.0 < s < 0.005):
                    s = s_old + random.gauss(0.0,0.001)
                    
                for k in range(len(OB)):
                    if S.Observables[k].data_type == 'ph':
                        # assign a new jitter to reference series
                        S.Observables[k].jitter = s
                
                #model inversion
                Glob.inhibit_msg = False
                joint_logL, id, info, RV_WEIGHT = S.simulated_annealing_local(S.Observables)
                
                bunch_data.append([s, joint_logL])
                    
            return bunch_data
            
        def parallel_inversion(N_inv=1):
            """
            """
            Glob.inhibit_msg = True
            Simulations = [StarsimSimulation() for i in range(N_inv)]
            #prepare the objects
            for Simulation in Simulations:
                #load the data and specific parameters for this data, in each object of the population
                OB = Simulation.read_multifit_config('./multifit.conf')
                #reset new map
                Simulation.spot_map[0] = new_spot_map[:]
                S.spot_map[0] = new_spot_map[:]
                
            pool = mp.Pool(Glob.ncpus)
            parallel_bunches = S.split_list(Simulations, Glob.ncpus)
            
            try:
                parallel_data_blocks = pool.map_async(inversion, parallel_bunches).get(9999999)
            except KeyboardInterrupt:
                io.print_message('Parallel pool killed. Bye!', index=3, color=31)
                pool.terminate()
                sys.exit()
                
            #join parallel_data_blocks
            SimulationsOut = []
            for parallel_data_block in parallel_data_blocks:
                for parallel_data_piece in parallel_data_block: 
                    SimulationsOut.append(parallel_data_piece)
            
            pool.terminate()
            
            return SimulationsOut
        
        
        #generate object simualtion
        S = StarsimSimulation()
        OB = S.read_multifit_config('./multifit.conf')
        
        global s_old
        s_old = 0.0025
        
        results = []
        for i in range( int( len(S.spot_map[0])/10-1)):
            for ii in range(5):
                #compute the FOM for each spot
                spot_map = S.spot_map[0]
            
                FOMs = []
                for j in range(len(spot_map)):
                    FOMs.append([j, spot_map[j][1]*sin((pi/180.0)*spot_map[j][2])*spot_map[j][4]**2])
            
            
                #sort the list of FOMs
                FOMs_sorted = sorted(FOMs, key=lambda x : x[1], reverse=False)
                less_significant_spot = FOMs_sorted[0]  #first in sorted list
            
                #remove the least significant spot
                new_spot_map = np.delete(spot_map, less_significant_spot[0], 0)
                S.spot_map[0] = new_spot_map[:]
        
            Glob.inhibit_msg = True
            #optimize again
            #compute the joint inverse problem (local) with all observables
            stats = parallel_inversion(N_inv=12)
            print (stats)
            #fit a polinomial
            stats = np.array(stats)
            #z = np.polyfit(stats[:,0], stats[:,1], 2)
            #s_max = -z[1]/(2.*z[0])
            #logL_max = z[0]*s_max**2 + z[1]*s_max + z[2]
            stats_sorted = sorted(stats, key=lambda x: x[1], reverse=True)[0]
            s_max, logL_max = stats_sorted[0], stats_sorted[1]
            s_old = s_max
            print ('s_max:', s_max, 'logL_max:', logL_max)
                
            #mean_joint_logL = np.mean(stats)
            #sigma = np.std(stats)
            
            results.append([len(new_spot_map), logL_max, s_max])
            
            print (len(new_spot_map), logL_max, s_max)
        
        io.write_file('FOM_jitter_150' + '.dat', results, ' ', '')
      
    elif mode == 'inversion':
        """multifit inversion
        """
        #generate object simualtion
        Simulation = StarsimSimulation()
        
        t_min, t_max = Simulation.dt1 - 2.*Simulation.spots_lifetime, Simulation.dt2
        
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
        except:
            io.print_message("Could not be read the multiobject configuration or no data selected", 33, 5)
            sys.exit()
        #modify stellar params and lc offset
        for object in OB:
            if object.data_type == 'ph':
                pass
                #object.set_ph_zero_point(object.z_ph)
        
        
        #fit observables
        res = Simulation.simulated_annealing(Simulation.Observables)
        print("")

        #copy the inverted map to the simulation environment 
        Simulation.spot_map[0] = Simulation.InvertedSpotMap 
                
        #compute the forward problem (HR)
        #generate curves
        if len(OB) > 0:
            for ob in OB:
                #iterate over the collection of Observables
                if ob.data_type == 'rv':
                    #RV object
                    print("\nRV object")
                    logL_rv = ob.compute_rv_series(write_rv=False)
                    print("\nStatistics:")
                    print("\tNull-model logL:", ob.logL_null) 
                    print("\tModel logL:", logL_rv)
                    print("\tDlogL(model-data):", logL_rv - ob.logL_null)
                    print("")
                    #write data series with added jitter
                    ob.write_rv_data()
                    #write data series + residuals (last col)
                    ob.compute_rv_series(write_rv=True)
                    #modify time array to generate HR curves
                    ob.obs_time = np.linspace(ob.dt1-10.0, ob.dt2+10.0, num=int(1440.0*(ob.dt2-ob.dt1)/ob.t_exp))
                    ob.obs_data = []
                    #compute the HR rv curve
                    ob.compute_rv_series(write_rv=True, mark='HR')
                
                elif ob.data_type == 'ph':
                    #Photometry object
                    print("\nPhotometry object:", ob.name)
                    logL_ph, spectra = ob.compute_lc_series(write_lc=False)
                    print("\nStatistics:")
                    print("\tNull-model logL:", ob.logL_null)
                    print("\tModel logL:", logL_ph)
                    print("\tDlogL(model-data):", logL_ph - ob.logL_null)
                    print("")
                    #generate data curve with z_ph implemented + jitter
                    #ob.write_lc_data()
                    # write data series + residuals (last col)
                    ob.compute_lc_series(write_lc=True)
                    #residuals
                    ph_residuals = []
                    for i in range(len(ob.obs_data)):
                        ph_residuals.append([ob.obs_data[i][0], ob.obs_data[i][1] - ob.mean_normalized_series[i][1],
                                             (ob.obs_data[i][1] - ob.mean_normalized_series[i][1])/ob.obs_data[i][2], 
                                             ob.obs_data[i][2], ob.name ])
                    
                    io.write_file("./output/PH_MODEL_" + str(ob.name) + ".residuals", ph_residuals, ' ', '')
                    ob.obs_time = np.linspace(ob.dt1-10.0, ob.dt2+10.0, num=int(1440.0*(ob.dt2-ob.dt1)/ob.t_exp))
                    ob.obs_data = []
                    #compute the HR lc curve
                    ob.compute_lc_series(write_lc=True, mark='HR_' + str(ob.name))
    
    elif mode == 'local_inversion':
        """multifit inversion
        """
        #generate object simualtion
        Simulation = StarsimSimulation()
        
        t_min, t_max = Simulation.dt1 - 2.*Simulation.spots_lifetime, Simulation.dt2
        
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
        except:
            print("Could not be read the multiobject configuration or no data selected")
            raise
        #modify stellar params and lc offset
        for object in OB:
            if object.data_type == 'ph':
                pass
                #object.set_ph_zero_point(object.z_ph)
        
        
        #fit observables
        res = Simulation.simulated_annealing_local(Simulation.Observables)
        print("")

        #copy the inverted map to the simulation environment 
        Simulation.spot_map[0] = Simulation.InvertedSpotMap 
                
        #compute the forward problem (HR)
        #generate curves
        if len(OB) > 0:
            for ob in OB:
                #iterate over the collection of Observables
                if ob.data_type == 'rv':
                    #RV object
                    print("\nRV object")
                    logL_rv = ob.compute_rv_series(write_rv=False)
                    print("\nStatistics:")
                    print("\tNull-model logL:", ob.logL_null) 
                    print("\tModel logL:", logL_rv)
                    print("\tDlogL(model-data):", logL_rv - ob.logL_null)
                    print("")
                    #write data series with added jitter
                    ob.write_rv_data()
                    #write data series + residuals (last col)
                    ob.compute_rv_series(write_rv=True)
                    #modify time array to generate HR curves
                    ob.obs_time = np.linspace(ob.dt1-10.0, ob.dt2+10.0, num=int(1440.0*(ob.dt2-ob.dt1)/ob.t_exp))
                    ob.obs_data = []
                    #compute the HR rv curve
                    ob.compute_rv_series(write_rv=True, mark='HR')
                
                elif ob.data_type == 'ph':
                    #Photometry object
                    print("\nPhotometry object")
                    logL_ph, spectra = ob.compute_lc_series(write_lc=False)
                    print("\nStatistics:")
                    print("\tNull-model logL:", ob.logL_null)
                    print("\tModel logL:", logL_ph)
                    print("\tDlogL(model-data):", logL_ph - ob.logL_null)
                    print("")
                    #generate data curve with z_ph implemented + jitter
                    ob.write_lc_data()
                    # write data series + residuals (last col)
                    ob.compute_lc_series(write_lc=True)
                    #residuals
                    ph_residuals = []
                    for i in range(len(ob.obs_data)):
                        ph_residuals.append([ob.obs_data[i][0], ob.obs_data[i][1] - ob.mean_normalized_series[i][1],
                                             (ob.obs_data[i][1] - ob.mean_normalized_series[i][1])/ob.obs_data[i][2], 
                                             ob.obs_data[i][2], ob.name ])
                    
                    io.write_file("./output/PH_MODEL_" + str(ob.name) + ".residuals", ph_residuals, ' ', '')
                    print (ob.dt1, ob.dt2)
                    ob.obs_time = np.linspace(ob.dt1-10.0, ob.dt2+10.0, num=int(1440.0*(ob.dt2-ob.dt1)/ob.t_exp))
                    ob.obs_data = []
                    #compute the HR lc curve
                    ob.compute_lc_series(write_lc=True, mark='HR_' + str(ob.name))
    
    elif mode == 'N-inversion':
        """fit multiple identic objects
        """        
        def split_list(alist, wanted_parts):
            """
            """
            length = len(alist)
            if length < wanted_parts:
                wanted_parts = length
                
            return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
                    for i in range(wanted_parts) ]
        
        def model_interpolation(time_data, y_data, time_model, y_model):
            """model,data interploation of arbitrarily sampled model
            """
            y_interpolated_model = []
            for j in range(len(y_data)):
                for i in range(len(y_model)-1):
                    if (time_model[i] <= time_data[j]) and (time_model[i+1] > time_data[j]):
                        interp_y = y_model[i] + ((y_model[i+1]-y_model[i])/(time_model[i+1]-time_model[i]))*(time_data[j]-time_model[i])
                        y_interpolated_model.append(interp_y)

            return np.array(y_interpolated_model)
                    
        def compute_joint_fit(parallel_bunch):
            """receives a list of Simulation objects, inverts observational data 
                and computes the model (forward problem)
            """        
            def find_nearest(array,value):
                idx = abs(np.abs(array)-value).argmin()
                return array[idx]
        
            bunch_data = []
            for S in parallel_bunch:
                """
                for object in S.Observables:
                    if object.data_type == 'ph':
                        object.set_ph_zero_point(1.0)
                    #recompute CCF for the new rotation period
                    elif object.data_type == 'rv':
                        object.init_ccfs()
                """    
                Glob.inhibit_msg = False   
                
                #init random number generator
                np.random.seed()
                
                #gen new random map
                #S.spot_map[0] = S.gen_rnd_map(S.spot_map[3])
                
                #MCSA model inversion
                joint_logL, id, info, RV_WEIGHT = S.simulated_annealing(S.Observables)
                
                #finally compute the log likelihood and time series of forward problem
                logL_ph = 0.0
                logL_rv = 0.0
                logL_bis = 0.0
            
                for i in range(len(S.Observables)):
                    if S.Observables[i].data_type == 'ph':
                            
                        logL_ph += S.Observables[i].compute_lc_series(write_lc=False)[0]
                        #compute HR curve
                        S.Observables[i].obs_time = np.linspace(S.Observables[i].dt1 - 10.0, S.Observables[i].dt2 + 10.0, 
                                                    num=int(1440.0*(S.Observables[i].dt2 - S.Observables[i].dt1)/S.Observables[i].t_exp))
                        
                        S.Observables[i].obs_data = []
                        #forward problem
                        S.Observables[i].compute_lc_series(write_lc=True, mark='')
                        
                    elif S.Observables[i].data_type == 'rv':
                        
                        logL_rv += S.Observables[i].compute_rv_series(write_rv=False)
                        #compute HR curve
                        S.Observables[i].obs_time = np.linspace(S.Observables[i].dt1 - 10.0, S.Observables[i].dt2 + 10.0, 
                                                    num=int(1440.0*(S.Observables[i].dt2 - S.Observables[i].dt1)/S.Observables[i].t_exp))
                        S.Observables[i].obs_data = []
                        #forward problem
                        S.Observables[i].compute_rv_series(write_rv=True)
                    
                    elif S.Observables[i].data_type == 'bis':
                        
                        logL_bis += S.Observables[i].compute_BIS_series(write_rv=False)
                        #compute HR curve
                        S.Observables[i].obs_time = np.linspace(S.Observables[i].dt1 - 10.0, S.Observables[i].dt2 + 10.0, 
                                                    num=int(1440.0*(S.Observables[i].dt2 - S.Observables[i].dt1)/S.Observables[i].t_exp))
                        S.Observables[i].obs_data = []
                        #forward problem
                        S.Observables[i].compute_BIS_series(write_rv=True)
                
                bunch_data.append([S, S.id, logL_ph, logL_rv, joint_logL])
                io.write_file('./maps' + 
                              '/spot_map_' + str(S.id) + '.dat', S.spot_map[0], ' ', '')
            
            return bunch_data
        
        
        #lock object
        lock = mp.Lock()
        #instantiate a Simulation object
        Sim = StarsimSimulation()
        
        #read data and configuration
        try:
            #list of objects to joint minimize 
            OB = Sim.read_multifit_config('./multifit.conf')
        except:
            print("Could not be read the multiobject configuration or no data selected")
            sys.exit()

        t_min, t_max = Sim.dt1 - 2.*Sim.spots_lifetime, Sim.dt2 
        
        #generate an initial population of StarsimSimulation objects
        INIT_POPULATION = 32
        io.print_message("Inversion of Observational Data", 0, 33)
        io.print_message("Generating initial population: " + str(INIT_POPULATION) + " simulations", 2, 94)
        
        #inhibit the screen messages, for cleanness
        Glob.inhibit_msg = True
        
        Simulations = [StarsimSimulation() for iter in range(INIT_POPULATION)]
        #prepare the objects
        for Simulation in Simulations:
            #inhibit the screen messages, for cleanness
            #Glob.inhibit_msg = True
            #load the data and specific parameters for this data, in each object of the population
            OB = Simulation.read_multifit_config('./multifit.conf')
            #Simulation.update_spectral_data(delta_t_sp=Delta_SP, delta_t_fc=Simulation.delta_t_fc)
            #check if photmetry object
            phot_objects = []
            for object in OB:
                if object.data_type == 'ph':
                    #select a photometric zero point
                    #object.set_ph_zero_point(object.z_ph)
                    phot_objects.append(object.info)
                    
        pool = mp.Pool(Glob.ncpus)
        parallel_bunches = split_list(Simulations, Glob.ncpus)
        
        try:
            parallel_data_blocks = pool.map_async(compute_joint_fit, parallel_bunches).get(9999999)
        except KeyboardInterrupt:
            io.print_message('Parallel pool killed. Bye!', index=3, color=31)
            pool.terminate()
            sys.exit()
        # free memory
        pool.terminate()
        
        print ("Joining all pieces from CPUs...")
        #join parallel_data_blocks+
        SimulationsOut = []
        for parallel_data_block in parallel_data_blocks:
            for parallel_data_piece in parallel_data_block: 
                SimulationsOut.append(parallel_data_piece)
        
        stats_ = np.array(SimulationsOut)
        
        #computing statisitcs
        mean = stats_[:,-1].astype(float).mean()
        std = stats_[:,-1].astype(float).std()
        p99 = np.percentile(stats_[:,-1], 99.9)
        maxLogL = max(stats_[:,-1])
        minLogL = min(stats_[:,-1])
        
        print ("Statistical logL:", mean + 3.*std)
        print ("Percentile .99:", p99)
        print ("Min.", minLogL)
        print ("Max.", maxLogL)
        print ("Mean:", mean)
        print ("St.dev:", std)
        print ("mean/std:", mean/std)
        
        #fitting statistics file
        try:
            io.write_file("./stats.dat", stats_, ' ', '')
        except:
            raise

        #compute (mean,sigma) series
        mean_models_ph = []
        mean_imm_models_ph = []
        mean_models_rv = []
        concat_models_rv = []
        mean_models_ff_sp = []
        z_stats = []
        
        for j in range(len(OB)):
            if OB[j].data_type == 'ph': 
                band_normean = []  #normalized w.r.t. its mean 
                band_normim = [] #normalized w.r.t. immaculate photosphere
                band_ff = []
                band_z = []
                for simulation in SimulationsOut:
                    times = np.array(simulation[0].Observables[j].mean_normalized_series)[:,0]
                    band_normean.append(np.array(simulation[0].Observables[j].mean_normalized_series)[:,1])
                    band_normim.append(np.array(simulation[0].Observables[j].imm_normalized_series)[:,1])
                    band_ff.append(np.array(simulation[0].Observables[j].ff_sp[:,1]))
                    #z's
                    band_z.append(simulation[0].Observables[j].z)
                                    
                band_normean, band_norim = np.array(band_normean), np.array(band_normim)
                band_ff = np.array(band_ff)
                normean, normean_std = np.mean(band_normean, axis=0), np.std(band_normean, axis=0)
                normim, normim_std = np.mean(band_norim, axis=0), np.std(band_normim, axis=0)
                #normalization over the mean and immaculate photosphere
                model_series_normean = list(zip(times, normean, normean_std))
                model_series_normim = list(zip(times, normim, normim_std))
                #z stats
                z_stats.append([str(OB[j].name), np.mean(band_z), np.std(band_z)])
                
                mean_models_ph.append(model_series_normean)
                mean_imm_models_ph.append(model_series_normim)
                
                #std_model = list(zip(times, std))
                io.write_file('./output/PH_MODEL_norm_imm_' + str(OB[j].name), model_series_normim, ' ', '')
                io.write_file('./output/PH_MODEL_norm_mean_' + str(OB[j].name), model_series_normean, ' ', '')
                
                #filling factor
                mean_ff, std_ff = np.mean(band_ff, axis=0), np.std(band_ff, axis=0)
                model_series_ff = list(zip(times, mean_ff, std_ff))
                
                mean_models_ff_sp.append(model_series_ff)
                
                io.write_file('./output/MODEL_FF_SP_MODEL_' + str(OB[j].name), model_series_ff, ' ', '')
            
            
            elif OB[j].data_type == 'rv': 
                band_rv = []
                for simulation in SimulationsOut:
                    times = np.array(simulation[0].Observables[j].normalized_series)[:,0]
                    band_rv.append(np.array(simulation[0].Observables[j].normalized_series)[:,1])
                
                band_rv = np.array(band_rv)
                mean, std = np.mean(band_rv, axis=0), np.std(band_rv, axis=0)
                model_series = list(zip(times, mean, std))
                mean_model = list(zip(times, mean))
                mean_models_rv.append(mean_model)
                concat_models_rv.extend(mean_model) #series of concatenated models
                std_model = list(zip(times, std))
                io.write_file('./output/MEAN_RV_MODEL_' + str(OB[j].name), mean_model, ' ', '')
                io.write_file('./output/STD_RV_MODEL_' + str(OB[j].name), std_model, ' ', '')    
                io.write_file('./output/MODEL_SERIES_FULL_' + str(OB[j].name), model_series, ' ', '') 
            
                MODEL_SERIES_REDUCED = []    
                for i in range(int(len(model_series)/100)):
                    MODEL_SERIES_REDUCED.append(model_series[100*i])
                
                io.write_file('./output/MODEL_SERIES_REDUCED_' + str(OB[j].name), MODEL_SERIES_REDUCED, ' ', '')
        
        io.write_file('z.dat', z_stats, ' ', '')
        #gen data with implemented z's
        for j in range(len(OB)):
            for i in range(len(z_stats)):
                if OB[j].name == str(z_stats[i][0]):
                    OB[j].set_ph_zero_point(float(z_stats[i][1]))
                    OB[j].write_lc_data()
        
        #compute logL stastistic (MEAN MODEL - DATA)
        ph_residuals = []
        rv_residuals = []
        bis_residuals = []
        
        ii = 0
        jj = 0
        #HR MEAN / RESIDUALS LOGL
        for j in range(len(OB)):
                if OB[j].data_type == 'ph': 
                    
                    mean_model = np.array(mean_models_ph[jj])
                    jj += 1
                    y_model, y_data_0 = mean_model[:,1].astype(float), np.array(OB[j].y_data_0)
                    y_data = np.array(OB[j].y_data)
                    time_model = mean_model[:,0].astype(float)
                    
                    inv_sigma2 = 1.0 / (OB[j].err**2 + OB[j].jitter**2)
                    inv_sigma2_0 = 1.0 / (OB[j].err_0**2 + OB[j].jitter**2)
                    
                    #interpolate densely sampled model
                    y_interpolated_model = model_interpolation(OB[j].obs_time_0, OB[j].obs_data_0, 
                                                                time_model, y_model)
                    
                    logL_0_ph = -0.5*(np.sum(((y_data - np.mean(y_data))**2)*inv_sigma2_0 + np.log(2.*np.pi) - np.log(inv_sigma2_0)) )
                    
                    logL_ph = -0.5*(np.sum(((y_data-y_interpolated_model)**2)*inv_sigma2 + np.log(2.*np.pi) - np.log(inv_sigma2)) )
        
                    for i in range(len(y_data)):
                        ph_residuals.append([OB[j].obs_time_0[i], y_data[i] - y_interpolated_model[i], OB[j].err[i] ])
                    
                    io.write_file("./output/MEAN_PH_MODEL_" + str(OB[j].name) + ".residuals", ph_residuals, ' ', '')
                    #write data series + offset + jitter
                    OB[j].write_lc_data()
                    print("\nPhotometry model", OB[j].name, ":")
                    print("\n\tNull model logL:", OB[j].logL_null)
                    print("\tlogL average photometric curve:", logL_ph)
                    print("\tDlogL:", logL_ph - OB[j].logL_null)
                    
                            
                elif OB[j].data_type == 'rv':
                    """
                    """
                    mean_model = np.array(mean_models_rv[ii])
                    ii += 1
                    y_model, y_data_0 = mean_model[:,1].astype(float), np.array(OB[j].y_data_0)
                    y_data = np.array(OB[j].y_data)
                    time_model = mean_model[:,0].astype(float)
                    
                    inv_sigma2 = 1.0 / (OB[j].err**2 + OB[j].jitter**2)
                    inv_sigma2_0 = 1.0 / (OB[j].err_0**2 + OB[j].jitter**2)
                    
                    #interpolate densely sampled model
                    y_interpolated_model = model_interpolation(OB[j].obs_time_0, OB[j].obs_data_0, 
                                                            time_model, y_model)
        
                    logL_0_rv = -0.5*(np.sum(((y_data - np.mean(y_data))**2)*inv_sigma2_0 + np.log(2.*np.pi) - np.log(inv_sigma2_0)) )
                    
                    logL_rv = -0.5*(np.sum(((y_data-y_interpolated_model)**2)*inv_sigma2 + np.log(2.*np.pi) - np.log(inv_sigma2)) ) 

                    for i in range(len(y_data)):
                        rv_residuals.append([OB[j].obs_time_0[i], y_data[i] - y_interpolated_model[i], OB[j].err[i] ])
                    
                    io.write_file("./output/MEAN_RV_MODEL_" + str(OB[j].name) + ".residuals", rv_residuals, ' ', '')
                    #write data series + jitter
                    OB[j].write_rv_data()
                    print("Radial Velocity:")
                    print("\n\tNull model logL:", OB[j].logL_null)
                    print("\tlogL average RV curve:", logL_rv)
                    print("\tDlogL:", logL_rv - OB[j].logL_null)
                    
                elif Simulation.Observables[i].data_type == 'bis':
                    """
                    """
                    y_model, y_data = mean_bis_model[:,1].astype(float), Simulation.Observables[i].y_data_0
                    time_model = mean_bis_model[:,0].astype(float)
                    inv_sigma2 = 1.0/(Simulation.Observables[i].err**2 + Simulation.Observables[i].jitter**2)
                    
                    #interpolate highly sampled model
                    y_interpolated_model = model_interpolation(Simulation.Observables[i].obs_time_0, Simulation.Observables[i].obs_data_0, 
                                                                time_model, y_model)
            
                    logL_0_bis = -0.5*(np.sum(((y_data)**2)*inv_sigma2 + np.log(2.*np.pi) - np.log(inv_sigma2)) )
                    logL_bis = -0.5*(np.sum(((y_data-y_interpolated_model)**2)*inv_sigma2 + np.log(2.*np.pi) - np.log(inv_sigma2)) ) 

                    for j in range(len(y_data)):
                        bis_residuals.append([Simulation.Observables[i].obs_time_0[j], y_data[j] - y_interpolated_model[j], Simulation.Observables[i].err[j] ])
                    
                    io.write_file("./output/MEAN_BIS_MODEL.residuals", rv_residuals, ' ', '')
                    #write data series + jitter
                    Simulation.Observables[i].write_bis_data()
                    print("Bisector span:")
                    print("\n\tNull model logL:", Simulation.Observables[0].logL_null)
                    print("\tlogL average BIS curve:", logL_bis)
                    print("\tDlogL:", logL_bis - Simulation.Observables[0].logL_null)
        
    elif mode == 'forward':
        """compute the observables from a list of spot maps and a list of params.
        (forward problem)
        """
        histo = []
        S = StarsimSimulation()
        
        try:
            #list of objects to joint minimize 
            OB = S.read_multifit_config('./multifit.conf')
        except:
            print("Could not be read the multiobject configuration or no data selected")
            raise
        
        maps = []
        maps.extend(glob.glob('./W52_reloaded/maps/spot_map*') )    
        
        #read the param. file
        param_data = io.read_file('./W52_reloaded/param_good.dat', ' ')
        
        #check lengths
        if len(maps) != len(param_data):
            print ("Warning: the number of maps is not equal to the parameter list")
            
        print(len(maps), "maps read...")
        print(len(param_data), "parameter lines read...")
        
        
        MODELS = [[] for j in range(len(OB))]
        FF = [[] for j in range(len(OB))]
        for map_ in maps:
            try:
                #list of objects to joint minimize 
                OB = S.read_multifit_config('./multifit.conf')
            except:
                print("Could not be read the multiobject configuration or no data selected")
                raise
            #read spot map
            spot_lst = S.read_file(map_, ' ')
            #copy the read map to the current simulation object
            S.spot_map[0] = spot_lst
            #look for the assoc. parameters to this map
            ref_code = str((map_.split('/')[-1]).split('_')[2].split('.')[0])
            
            for i in range(len(param_data)):
                if ref_code == str(param_data[i][0]):
                    param_lst = param_data[i]
            
            S.id = ref_code
            #load the parameters linked with the map
            S.P_rot = float(param_lst[2])
            # set Q ratio
            S.set_Q(float(param_lst[3]))
            # delta T spots
            delta_t_sp = float(param_lst[4])
            S.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=S.delta_t_fc)
      
            joint_logL = 0
            for i in range(len(OB)):
                if S.Observables[i].data_type == 'ph':
                    # assign a new jitter to reference series
                    S.Observables[i].jitter = float(param_lst[5])
           
                    joint_logL += S.Observables[i].compute_lc_series(write_lc=False)[0]
                    
                    try:
                        #S.Observables[i].write_lc_data()
                        #write data series + residuals (last column)
                        S.Observables[i].compute_lc_series(write_lc=False)
                    
                    except:
                        pass
                    
                    
                    #HR time sampling
                    S.Observables[i].obs_time = np.linspace(S.Observables[i].dt1-0.0, 
                                                            S.Observables[i].dt2+0.0, 
                                                            num=int(1440.0*(S.Observables[i].dt2-S.Observables[i].dt1)/S.Observables[i].t_exp))
                    S.Observables[i].obs_data = []
                    #compute the HR lc curve
                    S.Observables[i].compute_lc_series(write_lc=False, mark='HR_' + S.Observables[i].name)
                    
                    MODELS[i].append(S.Observables[i].mean_normalized_series)
                    FF[i].append(S.Observables[i].ff_sp)
                    
            #histo.append([S.id, joint_logL, S.P_rot, S.Q, float(param_lst[4]), float(param_lst[5])])    
            #print ("Joint logL:", S.id, joint_logL, S.P_rot, S.Q, float(param_lst[4]), float(param_lst[5]))
        
        #io.write_file('./histogram.dat', histo, ' ', '')
            
        
        MODELS = np.array(MODELS)
        FF = np.array(FF)
        
        for i in range(len(OB)):
            data_model = []
            data_ff = []
            mean_model, std_model = MODELS[i].mean(axis=0), MODELS[i].std(axis=0)  
            mean_ff, std_ff = FF[i].mean(axis=0), FF[i].std(axis=0)
            for j in range(len(mean_model)):
                data_model.append([mean_model[j][0], mean_model[j][1], std_model[j][1]])
                data_ff.append([mean_ff[j][0], mean_ff[j][1], std_ff[j][1]])
                
            io.write_file('./output/MEAN_PH_' + str(OB[i].name), data_model, ' ', '')
            io.write_file('./output/FF_' + str(OB[i].name), data_ff, ' ', '')
        
    elif mode == 'bayesian-UNDERTESTING':
        """perform bayesian optimization for the parameters
        """     
        
        def find_nearest(array,value):
            idx = abs(np.abs(array)-value).argmin()
            return array[idx]
        
        def compute_joint_fit(parallel_bunch):
            """receives a list of Simulation objects, inverts observational data 
               and computes forward problem
            """        

            bunch_data = []
            for S in parallel_bunch:
                for i in range(len(S.Observables)):
                    if S.Observables[i].data_type == 'ph':
                        S.Observables[i].set_ph_zero_point(S.Observables[i].z_ph)
                    #recompute CCF for the new rotation period
                    elif S.Observables[i].data_type == 'rv':
                        S.Observables[i].init_ccfs()
                
                #gen intit random spot map
                S.spot_map[0] = S.gen_rnd_map(S.n_spots)
                
                #model inversion   
                joint_logL = S.simulated_annealing(S.Observables)[0]
                     
                bunch_data.append(joint_logL)
              
            return bunch_data
            
        def multi_inversion(P_rot, Q, delta_T_sp):
            """
            """
            for Simulation in Simulations:
                #set new parameters
                delta_t_sp = find_nearest(temps, delta_T_sp)
                Simulation.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=Simulation.delta_t_fc)
                Simulation.P_rot = P_rot
                Simulation.set_Q(Q)
                
                #update objects
                for i in range(len(Simulation.Observables)):
                    if Simulation.Observables[i].data_type == 'ph':
                        pass
                        #Simulation.Observables[i].set_ph_zero_point(Simulation.Observables[i].z_ph)
                    #recompute CCF for the new rotation period
                    elif Simulation.Observables[i].data_type == 'rv':
                        Simulation.Observables[i].init_ccfs()
            
            pool = mp.Pool(Glob.ncpus)
            parallel_bunches = Simulation.split_list(Simulations, Glob.ncpus)
            
            try:
                parallel_data_blocks = pool.map_async(compute_joint_fit, parallel_bunches).get(9999999)
            except KeyboardInterrupt:
                io.print_message('Parallel pool killed. Bye!', index=3, color=31)
                pool.terminate()
                sys.exit()
                
            #join parallel_data_blocks
            SimulationsOut = []
            for parallel_data_block in parallel_data_blocks:
                for parallel_data_piece in parallel_data_block: 
                    SimulationsOut.append(parallel_data_piece)
            
            
            raw_points.append( [P_rot, Q, delta_T_sp, np.mean(SimulationsOut), np.std(SimulationsOut)] )
            
            return max(SimulationsOut)
        
        def inversion(P_rot, Q, delta_T_sp):
            """
            """
            #set new parameters
            delta_t_sp = find_nearest(temps, delta_T_sp)
            Simulation.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=Simulation.delta_t_fc)
            Simulation.P_rot = P_rot
            Simulation.set_Q(Q)
            
            #model inversion   
            joint_logL = S.simulated_annealing(S.Observables)[0]
            
            return joint_logL
   
        
        #generate object simualtion
        S = StarsimSimulation()
        
        try:
            #list of objects to joint minimize 
            OB = S.read_multifit_config('./multifit.conf')
        except:
            print("Could not be read the multiobject configuration or no data selected")
            raise
        
        io.print_message("Inversion of Observational Data", 2, 91)
        Glob.inhibit_msg = True
        
        Simulations = [StarsimSimulation() for i in range(Glob.ncpus)]
        #prepare the objects
        for Simulation in Simulations:
            #load the data and specific parameters for this data, in each object of the population
            OB = Simulation.read_multifit_config('./multifit.conf')
        
        
        from bayes_opt import BayesianOptimization
        import warnings
        warnings.filterwarnings("ignore")
        
        raw_points = []
        
        temps = [100.0 + 25.0*i for i in range(60)]
        
        params =  {'P_rot': (16.0, 20.0), 'Q': (0.0, 3.0), 'delta_T_sp': (50.0, 2000.0)}
        
        """
        for j in range(10):
        
            print('Iteration:', j)
            randP, randQ, randT = random.uniform(params['P_rot'][0], params['P_rot'][1]), \
                                  random.uniform(params['Q'][0], params['Q'][1]), \
                                  random.uniform(params['delta_T_sp'][0], params['delta_T_sp'][1])
     
            inversion(randP, randQ, randT)
        
        
        io.write_file('raw_points_' + str(int(random.uniform(0,1000))) + '.dat', raw_points, ' ', '')
        """
        
        
        bo = BayesianOptimization(inversion, params)
    
        bo.maximize(init_points=15, n_iter=5, acq='ei', xi=1e-5)
        
        bo.fit(bo.X, bo.Y)
        
        
        net_to_predict = []
        for i in range(0,50):
            x = 10.0 + i/3.0
            for j in range(0,50):
                y = j/10.0
                net_to_predict.append([x,y])
    
        mu, sigma = bo.gp.predict(net_to_predict, return_std=True)
       
        results = []
        for i in range(len(mu)):
            results.append([net_to_predict[i][0], net_to_predict[i][1],  mu[i]])

        SortedCandidates = sorted(results, key=lambda x: x[2], reverse=True)
        print "Max.", SortedCandidates[0]
    
        io.write_file('GaussianSurface.dat', results, ' ', '')
        
    elif mode == 'param_fit':
        """fit theta parameters of the inverse problem
        """
        def flatten(lis):
            """return a nested list, flattened
            """
            f_lis = []
            for item in lis:
                if type(item) == type([]):
                    f_lis.extend(flatten(item))
                else:
                    f_lis.append(item)
            
            return f_lis
        
        def find_nearest(array, value):
            idx = abs(np.abs(array)-value).argmin()
            return array[idx]
        
        def inversion(parallel_bunch):
            """
            """ 
            bunch_data = []
            
            for S in parallel_bunch: 
                
                #reset the random generator
                np.random.seed()
                # set uniformly distributed parameters on its defined limits
                S.P_rot = random.uniform(L_P_rot[0], L_P_rot[1])
                # set Q ratio
                S.set_Q(random.uniform(L_Q[0], L_Q[1]))
                # delta T spots
                delta_t_sp_ = random.uniform(L_T_sp[0], L_T_sp[1])
                delta_t_sp = find_nearest([50.0 + 25.0*j for j in range(int((L_T_sp[1]-50.0)/25.0))], delta_t_sp_)
                #lock
                l.acquire()
                S.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=Simulation.delta_t_fc)
                l.release()
                #z_ph_REF = random.uniform(L_z_ph[0], L_z_ph[1])
                
                
                """
                for i in range(len(PH_OBJS)):
                    if S.Observables[i].name == BEST_SERIES:
                        # set a new zero-point to reference series
                        S.Observables[i].set_ph_zero_point(z_ph_REF)
                        # assign a new jitter to reference series
                        #S.Observables[i].jitter = ph_jitter_VAR
                        # model inversion (with only the reference photometric series)
                        joint_logL_band, SA_id, info, RV_WEIGHT = S.simulated_annealing( [S.Observables[i]] )
                    
                #compute the forward problem in ph and set the new computed z_ph's
                zphs = []
                ph_jitters = []
                for Obs in S.Observables:
                    if Obs.data_type == 'ph' and Obs.name != BEST_SERIES:
                        Obs.compute_lc_series(write_lc=False)
                        # compute the other zero-points using 1/<f(t)>
                        z_ph = 1.0 / Obs.mean_series
                        zphs.append(z_ph)
                        # set new zero-point to current Observable
                        Obs.set_ph_zero_point(z_ph)
                        # assign new jitter to current Observable 
                        #Obs.jitter = ph_jitter_VAR #random.uniform(L_jitter_ph[0], L_jitter_ph[1])
                        #ph_jitters.append(Obs.jitter)
                """
            
                #compute the joint inverse problem with all observables
                joint_logL, SA_id, info, RV_WEIGHT = S.simulated_annealing(S.Observables)
                
                out = flatten([SA_id, joint_logL, S.P_rot, S.Q, S.delta_t_sp])
                bunch_data.append(out)
                
                #save the spotmap to disk
                io.write_file('./maps' + 
                              '/spot_map_' + str(SA_id) + '.dat', S.spot_map[0], ' ', '')
                
                
            return bunch_data
        
        
        l = mp.Lock()
        #instantiate a Simulation object
        Simulation = StarsimSimulation()
        
        t_min, t_max = Simulation.dt1 - 2.*Simulation.spots_lifetime, Simulation.dt2 
       
        #generate an initial population of StarsimSimulation objects
        INIT_POPULATION = 16
        io.print_message("Inversion of Observational Data", 2, 91)
        io.print_message("Generating initial population: " + str(INIT_POPULATION) + " simulations", 2, 94)
        
        #inhibit the screen messages, for cleanness
        Glob.inhibit_msg = True
        
        Simulations = [StarsimSimulation() for i in range(INIT_POPULATION)]
        #prepare the objects
        for Simulation in Simulations:
            #load the data and specific parameters for this data, in each object of the population
            OB = Simulation.read_multifit_config('./multifit.conf')
        
        #group photometric objects
        PH_OBJS = []
        for Obj in OB:
            if Obj.data_type == 'ph': 
                PH_OBJS.append(Obj)
        
        
        #-------------------------------------------------
        #PARAMs LIMITS
        #L_z_ph = [1.0, 1.04]
        L_T_sp = [50.0, 2000.0]
        L_P_rot = [14.99, 15.01]
        L_Q = [0.0, 1.0]
        #L_jitter_ph = [0.0005, 0.0075]
        #L_inc = [20.0, 90.0]
        #-------------------------------------------------
        
        print("Initial population ready. Computing inverse problem...")
        
        for i in range(1000):
            pool = mp.Pool(Glob.ncpus)
            parallel_bunches = Simulation.split_list(Simulations, Glob.ncpus)
     
            try:
                parallel_data_blocks = pool.map_async(inversion, parallel_bunches).get(9999999)
            except KeyboardInterrupt:
                io.print_message('Parallel pool killed. Bye!', index=3, color=31)
                pool.terminate()
                sys.exit()
            #free memory
            pool.terminate()
           
            # join parallel_data_blocks
            SimulationsOut = []
            for parallel_data_block in parallel_data_blocks:
                for parallel_data_piece in parallel_data_block: 
                    SimulationsOut.append(parallel_data_piece)
    
            io.write_file('./param_fit/Theta_' + str( int(random.uniform(0,10000)) ) + '.dat', SimulationsOut, ' ', '')
    
    elif mode == 'N-depth':
        """
        """
        def get_ARIEL_resolution(spectra):
            """
            """
            ARIEL_SPECTRA = []

            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(len(spectra)):
                if spectra[i][0] > 500.0 and spectra[i][0] < 600.0:
                    #vis channel
                    f_avg += spectra[i][1]
                    lam_avg += spectra[i][0]
                    n += 1
            
            ARIEL_SPECTRA.append([0.55, f_avg/n, 0.05])
            
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(len(spectra)):
                if spectra[i][0] > 600.0 and spectra[i][0] < 800.0:
                    #FGS1 channel
                    f_avg += spectra[i][1]
                    lam_avg += spectra[i][0]
                    n += 1
            
            ARIEL_SPECTRA.append([0.7, f_avg/n, 0.1])
            
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(len(spectra)):
                if spectra[i][0] > 850.0 and spectra[i][0] < 1100.0:
                    #FGS2 channel
                    f_avg += spectra[i][1]
                    lam_avg += spectra[i][0]
                    n += 1
            
            ARIEL_SPECTRA.append([0.975, f_avg/n, 0.125])
            
            #NIR Channel
            NIR = []
            for i in range(len(spectra)):
                if spectra[i][0] > 1100.0 and spectra[i][0] < 1950.0:
                    #NIR channel
                    NIR.append([spectra[i][0], spectra[i][1]])
            NIR_channel = []
            i_init = 0
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(i_init, len(NIR)):
                if NIR[i][0] - NIR[i_init][0] <= 100.0:
                    lam_avg += NIR[i][0]
                    f_avg += NIR[i][1]
                    n += 1
                else:    
                    NIR_channel.append([0.001*lam_avg/n, f_avg/n, 0.05])
                    i_init = i
                    n = 0    
                    f_avg = 0.0
                    lam_avg = 0.0
            
            ARIEL_SPECTRA.extend(NIR_channel)
            
            #Ch0 Channel
            Ch0 = []
            for i in range(len(spectra)):
                if spectra[i][0] > 1950.0 and spectra[i][0] < 3950.0:
                    Ch0.append([spectra[i][0], spectra[i][1]])
            Ch0_channel = []
            i_init = 0
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(i_init, len(Ch0)):
                if Ch0[i][0] - Ch0[i_init][0] <= 30.0:
                    lam_avg += Ch0[i][0]
                    f_avg += Ch0[i][1]
                    n += 1
                else:    
                    Ch0_channel.append([0.001*lam_avg/n, f_avg/n, 0.015])
                    i_init = i
                    n = 0    
                    f_avg = 0.0
                    lam_avg = 0.0
            
            ARIEL_SPECTRA.extend(Ch0_channel)
            
            #Ch1 Channel
            Ch1 = []
            for i in range(len(spectra)):
                if spectra[i][0] > 3950.0 and spectra[i][0] < 7850.0:
                    Ch1.append([spectra[i][0], spectra[i][1]])
            Ch1_channel = []
            i_init = 0
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(i_init, len(Ch1)):
                if Ch1[i][0] - Ch1[i_init][0] <= 180.0:
                    lam_avg += Ch1[i][0]
                    f_avg += Ch1[i][1]
                    n += 1
                else:    
                    Ch1_channel.append([0.001*lam_avg/n, f_avg/n, 0.091])
                    i_init = i
                    n = 0    
                    f_avg = 0.0
                    lam_avg = 0.0
            
            ARIEL_SPECTRA.extend(Ch1_channel)
            
                    
            return ARIEL_SPECTRA        

        def parallel_fit(parallel_bunch):
            """
            """
            for S in parallel_bunch:
                
                print (S.id, S.P_rot)
                Glob.planet = False
                #fit observables
                while (1):
                    res = S.simulated_annealing(S.Observables)
                    logL_map = res[0]
                    #print("")
                    #print ("Map generated with logL:", logL_map)
                    if logL_map > 1630.6:
                        #print ("Map ACCEPTED.")
                        break
                
                #SPOTTED PHOTOSPHERE
                #create a photometry object
                Ph = CPhotometry(S)
                #copy the inverted map to the simulation environment 
                S.spot_map[0] = S.InvertedSpotMap
        
                #transit duration (days)
                T_dur = (S.P_planet/np.pi)*asin(sqrt((S.R_star+S.R_planet)**2 - S.bim**2)/S.A_planet)
                N_transits = int((Ph.dt2 - Ph.dt1) / S.P_planet )
                t_init_transit = (int(((Ph.dt1 - S.T0_planet)/ S.P_planet)) + 1)*S.P_planet + S.T0_planet 

                #print("Number of transit centers:", N_transits)
                #print("Computing the transit events considering activity")
                #print("Number of transits:", int(N_transits))
            
                #iterate over all transit epochs: for each transit, compute a normalized depth curve
                for i in range(int(N_transits)):
                    #mid-transit time
                    time_stamp = t_init_transit + i*S.P_planet
                    #print("Time of transit:", time_stamp)    
                    Ph.obs_time = [time_stamp]
        
                    #no planet
                    Glob.planet = True
                    S.R_planet = 0.0
                    #print("Computing activity signal:")
                    logL, spectra_act = Ph.compute_lc_series(write_lc=False, mark='activity')
                    
                    #transits
                    Glob.planet = True
                    S.R_planet = R_planet
                    #print("Computing transits + activity:")
                    logL, spectra_act_planet = Ph.compute_lc_series(write_lc=False, mark='activity+transits_' + str(i) )
                    
                    #the first element of spectra arrays is the epoch (JD). The rest, fluxes.               
                    Depth_series = []
                    for j in range(len(Ph.wv)):
                        Depth_series.append( [float(Ph.wv[j]), 
                                            (float(spectra_act[0][j+1])-float(spectra_act_planet[0][j+1]))/float(spectra_act[0][j+1])] )         
                    
                    Depth_series_ARIEL = get_ARIEL_resolution(Depth_series)
                    io.write_file('./output/SPOTTED+PL_' + str(time_stamp) + '_' + str(logL_map) + '.dat', Depth_series_ARIEL, ' ', '')
                    
                    #compute the normalized transit depth
                    Depth_activity_corrected = []
                    for j in range(len(Depth_series)):
                        Depth_activity_corrected.append( [Depth_series[j][0], Depth_series[j][1] - Depth_imm_series[j][1] ] )
                    
                    #write Depth-time array on a file
                    Depth_activity_corrected_ARIEL = get_ARIEL_resolution(Depth_activity_corrected)
                    io.write_file('./output/norm_depth_' + str(time_stamp) + '_' + str(logL_map) + '_' + str(S.id) +'.dat', Depth_activity_corrected_ARIEL, ' ', '')
                    
    
        l = mp.Lock()    
        #generate object simulation
        Simulation = StarsimSimulation()
        
        R_planet = Simulation.R_planet
        
        T0_pl = Simulation.T0_planet
        
        t_min, t_max = Simulation.dt1 - 2.*Simulation.spots_lifetime, Simulation.dt2
        
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
        except:
            print("Could not be read the multiobject configuration or no data selected")
            raise
        
        #create a photometry object
        Ph = CPhotometry(Simulation)

        #IMMACULATE PHOTOSPHERE
        print("Computing the transit events considering immaculate photosphere")
        #reset the epochs to avoid the spot map
        #Ph.dt2 = Ph.dt2 - Ph.dt1
        dt1, dt2 = Ph.dt1, Ph.dt2
        
        Ph.dt1 = 0.0 
        Simulation.T0_planet = Ph.dt1
    
        Ph.obs_time = []
        for i in range(1,2):
            Ph.obs_time = [Ph.dt1 + i*Simulation.P_planet]
    
        #compute lcs
        Glob.planet = True  #use grid integration method
        #immaculate photosphere without transits
        Simulation.R_planet = 0.0
        print("Computing immaculate signal:")
        logL, spectra_imm = Ph.compute_lc_series(write_lc=False, mark='immaculate')
        
        #include transiting planet
        Glob.planet = True
        Simulation.R_planet = R_planet
        #immaculate photosphere + planet
        print("Computing transits with immaculate photosphere:")
        logL, spectra_imm_planet = Ph.compute_lc_series(write_lc=False, mark='immaculate+transits')
    
        Depth_imm_series = []
        for j in range(len(Ph.wv)):
            Depth_imm_series.append( [float(Ph.wv[j]), 
                                     (float(spectra_imm[0][j+1])-float(spectra_imm_planet[0][j+1]))/float(spectra_imm[0][j+1])] )
        
        Depth_imm_series_ARIEL = get_ARIEL_resolution(Depth_imm_series)
        io.write_file('./IMM.dat', Depth_imm_series_ARIEL, ' ', '')
        
        Simulation.T0_planet = T0_pl
        
        print("Immaculate depth curve computed")
        
        
        #read the param. file
        param_data = io.read_file('./W52_reloaded/param_good.dat', ' ')
        
        #generate an initial population of StarsimSimulation objects
        INIT_POPULATION = len(param_data)
        io.print_message("Inversion of Observational Data", 2, 91)
        io.print_message("Generating initial population: " + str(INIT_POPULATION) + " simulations", 2, 94)
        
        #inhibit the screen messages, for cleanness
        Glob.inhibit_msg = True
        
        Simulations = [StarsimSimulation() for i in range(INIT_POPULATION)]
        
        #prepare the objects
        k = 0
        for S in Simulations:
            #load the data and specific parameters for this data, in each object of the population
            OB = S.read_multifit_config('./multifit.conf')
            S.id = int(param_data[k][0])
            #load the parameters linked with the map
            S.P_rot = float(param_data[k][2])
            # set Q ratio
            S.set_Q(float(param_data[k][3]))
            # delta T spots
            delta_t_sp = float(param_data[k][4])
            l.acquire()
            S.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=S.delta_t_fc)
            l.release()
            
            for i in range(len(OB)):
                if S.Observables[i].data_type == 'ph':
                    # assign a new jitter to reference series
                    S.Observables[i].jitter = float(param_data[k][5])
            
            k += 1
        
        pool = mp.Pool(Glob.ncpus)
        parallel_bunches = Simulation.split_list(Simulations, Glob.ncpus)
    
        try:
            parallel_data_blocks = pool.map_async(parallel_fit, parallel_bunches).get(9999999)
        except KeyboardInterrupt:
            io.print_message('Parallel pool killed. Bye!', index=3, color=31)
            pool.terminate()
            sys.exit()
        #free memory
        pool.terminate()
    
    elif mode == 'multiple_inversion':
        """
        """
        def parallel_fit(parallel_bunch):
            """
            """
            for S in parallel_bunch:
                print (S.id, S.P_rot)
                #fit observables
                while (1):
                    #Glob.inhibit_msg = False
                    res = S.simulated_annealing(S.Observables)
                    logL_map = res[0]
                    #print("")
                    #print ("Map generated with logL:", logL_map)
                    if logL_map > 1630.6:
                        #print ("Map ACCEPTED.")
                        break
                
                io.write_file('./maps' + 
                            '/spot_map_' + str(S.id) + '.dat', S.spot_map[0], ' ', '')
                    
        
        Simulation = StarsimSimulation()
        
        #gen Lock object
        l = mp.Lock()
        #read the param. file
        param_data = io.read_file('./W52_reloaded/param_good.dat', ' ')
        
        #generate an initial population of StarsimSimulation objects
        INIT_POPULATION = len(param_data)
        io.print_message("Inversion of Observational Data", 2, 91)
        io.print_message("Generating initial population: " + str(INIT_POPULATION) + " simulations", 2, 94)
        
        #inhibit the screen messages, for cleanness
        Glob.inhibit_msg = True
        
        Simulations = [StarsimSimulation() for i in range(INIT_POPULATION)]
        
        #prepare the objects
        k = 0
        for S in Simulations:
            #load the data and specific parameters for this data, in each object of the population
            OB = S.read_multifit_config('./multifit.conf')
            S.id = int(param_data[k][0])
            #load the parameters linked with the map
            S.P_rot = float(param_data[k][2])
            # set Q ratio
            S.set_Q(float(param_data[k][3]))
            # delta T spots
            delta_t_sp = float(param_data[k][4])
            l.acquire()
            S.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=S.delta_t_fc)
            l.release()
            
            for i in range(len(OB)):
                if S.Observables[i].data_type == 'ph':
                    # assign a new jitter to reference series
                    S.Observables[i].jitter = float(param_data[k][5])
            
            k += 1
        
        pool = mp.Pool(Glob.ncpus)
        parallel_bunches = Simulation.split_list(Simulations, Glob.ncpus)
    
        try:
            parallel_data_blocks = pool.map_async(parallel_fit, parallel_bunches).get(9999999)
        except KeyboardInterrupt:
            io.print_message('Parallel pool killed. Bye!', index=3, color=31)
            pool.terminate()
            sys.exit()
        #free memory
        pool.terminate()
    
    elif mode == 'N-depth_deprecated':
        """
        """
        def get_ARIEL_resolution(spectra):
            """
            """
            ARIEL_SPECTRA = []

            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(len(spectra)):
                if spectra[i][0] > 500.0 and spectra[i][0] < 600.0:
                    #vis channel
                    f_avg += spectra[i][1]
                    lam_avg += spectra[i][0]
                    n += 1
            
            ARIEL_SPECTRA.append([0.55, f_avg/n, 0.05])
            
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(len(spectra)):
                if spectra[i][0] > 600.0 and spectra[i][0] < 800.0:
                    #FGS1 channel
                    f_avg += spectra[i][1]
                    lam_avg += spectra[i][0]
                    n += 1
            
            ARIEL_SPECTRA.append([0.7, f_avg/n, 0.1])
            
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(len(spectra)):
                if spectra[i][0] > 850.0 and spectra[i][0] < 1100.0:
                    #FGS2 channel
                    f_avg += spectra[i][1]
                    lam_avg += spectra[i][0]
                    n += 1
            
            ARIEL_SPECTRA.append([0.975, f_avg/n, 0.125])
            
            #NIR Channel
            NIR = []
            for i in range(len(spectra)):
                if spectra[i][0] > 1100.0 and spectra[i][0] < 1950.0:
                    #NIR channel
                    NIR.append([spectra[i][0], spectra[i][1]])
            NIR_channel = []
            i_init = 0
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(i_init, len(NIR)):
                if NIR[i][0] - NIR[i_init][0] <= 100.0:
                    lam_avg += NIR[i][0]
                    f_avg += NIR[i][1]
                    n += 1
                else:    
                    NIR_channel.append([0.001*lam_avg/n, f_avg/n, 0.05])
                    i_init = i
                    n = 0    
                    f_avg = 0.0
                    lam_avg = 0.0
            
            ARIEL_SPECTRA.extend(NIR_channel)
            
            #Ch0 Channel
            Ch0 = []
            for i in range(len(spectra)):
                if spectra[i][0] > 1950.0 and spectra[i][0] < 3950.0:
                    Ch0.append([spectra[i][0], spectra[i][1]])
            Ch0_channel = []
            i_init = 0
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(i_init, len(Ch0)):
                if Ch0[i][0] - Ch0[i_init][0] <= 30.0:
                    lam_avg += Ch0[i][0]
                    f_avg += Ch0[i][1]
                    n += 1
                else:    
                    Ch0_channel.append([0.001*lam_avg/n, f_avg/n, 0.015])
                    i_init = i
                    n = 0    
                    f_avg = 0.0
                    lam_avg = 0.0
            
            ARIEL_SPECTRA.extend(Ch0_channel)
            
            #Ch1 Channel
            Ch1 = []
            for i in range(len(spectra)):
                if spectra[i][0] > 3950.0 and spectra[i][0] < 7850.0:
                    Ch1.append([spectra[i][0], spectra[i][1]])
            Ch1_channel = []
            i_init = 0
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(i_init, len(Ch1)):
                if Ch1[i][0] - Ch1[i_init][0] <= 180.0:
                    lam_avg += Ch1[i][0]
                    f_avg += Ch1[i][1]
                    n += 1
                else:    
                    Ch1_channel.append([0.001*lam_avg/n, f_avg/n, 0.091])
                    i_init = i
                    n = 0    
                    f_avg = 0.0
                    lam_avg = 0.0
            
            ARIEL_SPECTRA.extend(Ch1_channel)
            
                    
            return ARIEL_SPECTRA        

        #generate object simulation
        Simulation = StarsimSimulation()
        
        R_planet = Simulation.R_planet
        
        T0_pl = Simulation.T0_planet
        
        t_min, t_max = Simulation.dt1 - 2.*Simulation.spots_lifetime, Simulation.dt2
        
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
        except:
            print("Could not be read the multiobject configuration or no data selected")
            raise
        
        #modify stellar params and lc offset
        for object in OB:
            if object.data_type == 'ph':
                object.set_ph_zero_point(object.z_ph)
        
        #create a photometry object
        Ph = CPhotometry(Simulation)

        #IMMACULATE PHOTOSPHERE
        print("Computing the transit events considering immaculate photosphere")
        #reset the epochs to avoid the spot map
        #Ph.dt2 = Ph.dt2 - Ph.dt1
        dt1, dt2 = Ph.dt1, Ph.dt2
        
        Ph.dt1 = 0.0 
        Simulation.T0_planet = Ph.dt1
    
        Ph.obs_time = []
        for i in range(1,2):
            Ph.obs_time = [Ph.dt1 + i*Simulation.P_planet]
    
        #compute lcs
        Glob.planet = True  #use grid integration method
        #immaculate photosphere without transits
        Simulation.R_planet = 0.0
        print("Computing immaculate signal:")
        logL, spectra_imm = Ph.compute_lc_series(write_lc=False, mark='immaculate')
        
        #include transiting planet
        Glob.planet = True
        Simulation.R_planet = R_planet
        #immaculate photosphere + planet
        print("Computing transits with immaculate photosphere:")
        logL, spectra_imm_planet = Ph.compute_lc_series(write_lc=False, mark='immaculate+transits')
    
        Depth_imm_series = []
        for j in range(len(Ph.wv)):
            Depth_imm_series.append( [float(Ph.wv[j]), 
                                     (float(spectra_imm[0][j+1])-float(spectra_imm_planet[0][j+1]))/float(spectra_imm[0][j+1])] )
        
        Depth_imm_series_ARIEL = get_ARIEL_resolution(Depth_imm_series)
        io.write_file('./IMM.dat', Depth_imm_series_ARIEL, ' ', '')
        
        Simulation.T0_planet = T0_pl
        
        print("Immaculate depth curve computed")
        
        #run a number of independent fits (maps) and compute depth curves with activity
        for k in range(50):
            
            #Ph.dt1, Ph.dt2 = dt1, dt2
            Glob.planet = False
            #fit observables
            res = Simulation.simulated_annealing(Simulation.Observables)
            logL_map = res[0]
            print("")
           
            #SPOTTED PHOTOSPHERE
            #create a photometry object
            Ph = CPhotometry(Simulation)
            #copy the inverted map to the simulation environment 
            Simulation.spot_map[0] = Simulation.InvertedSpotMap
    
            #transit duration (days)
            T_dur = (Simulation.P_planet/np.pi)*asin(sqrt((Simulation.R_star+Simulation.R_planet)**2 - Simulation.bim**2)/Simulation.A_planet)
            N_transits = int((Ph.dt2 - Ph.dt1) / Simulation.P_planet )
            t_init_transit = (int(((Ph.dt1 - Simulation.T0_planet)/ Simulation.P_planet)) + 1)*Simulation.P_planet + Simulation.T0_planet 

            print("Number of transit centers:", N_transits)
            
            print("Computing the transit events considering activity")
        
            print("Number of transits:", int(N_transits))
         
            #iterate over all transit epochs: for each transit, compute a normalized depth curve
            for i in range(int(N_transits)):
                #mid-transit time
                time_stamp = t_init_transit + i*Simulation.P_planet
                print("Time of transit:", time_stamp)    
                Ph.obs_time = [time_stamp]
    
                #no planet
                Glob.planet = True
                Simulation.R_planet = 0.0
                print("Computing activity signal:")
                logL, spectra_act = Ph.compute_lc_series(write_lc=False, mark='activity')
                 
                #transits
                Glob.planet = True
                Simulation.R_planet = R_planet
                print("Computing transits + activity:")
                logL, spectra_act_planet = Ph.compute_lc_series(write_lc=False, mark='activity+transits_' + str(i) )
                
                #the first element of spectra arrays is the epoch (JD). The rest, fluxes.               
                Depth_series = []
                for j in range(len(Ph.wv)):
                    Depth_series.append( [float(Ph.wv[j]), 
                                         (float(spectra_act[0][j+1])-float(spectra_act_planet[0][j+1]))/float(spectra_act[0][j+1])] )         
                
                Depth_series_ARIEL = get_ARIEL_resolution(Depth_series)
                io.write_file('./output_W52/SPOTTED+PL_' + str(time_stamp) + '_' + str(logL_map) + '.dat', Depth_series_ARIEL, ' ', '')
                
                #compute the normalized transit depth
                Depth_activity_corrected = []
                for j in range(len(Depth_series)):
                    Depth_activity_corrected.append( [Depth_series[j][0], Depth_series[j][1] - Depth_imm_series[j][1] ] )
                
                #write Depth-time array on a file
                Depth_activity_corrected_ARIEL = get_ARIEL_resolution(Depth_activity_corrected)
                io.write_file('./output_W52/norm_depth_' + str(time_stamp) + '_' + str(logL_map) + '.dat', Depth_activity_corrected_ARIEL, ' ', '')

    elif mode == 'static_N-depth':
        """do not invert the light curves. 
           Spectra resolution: high (BT_Settl_Ariel)
           Filter: flat.dat
        """
        def get_ARIEL_resolution(spectra):
            """
            """
            ARIEL_SPECTRA = []

            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(len(spectra)):
                if spectra[i][0] > 500.0 and spectra[i][0] < 600.0:
                    #vis channel
                    f_avg += spectra[i][1]
                    lam_avg += spectra[i][0]
                    n += 1
            
            ARIEL_SPECTRA.append([0.55, f_avg/n, 0.05])
            
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(len(spectra)):
                if spectra[i][0] > 600.0 and spectra[i][0] < 800.0:
                    #FGS1 channel
                    f_avg += spectra[i][1]
                    lam_avg += spectra[i][0]
                    n += 1
            
            ARIEL_SPECTRA.append([0.7, f_avg/n, 0.1])
            
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(len(spectra)):
                if spectra[i][0] > 850.0 and spectra[i][0] < 1100.0:
                    #FGS2 channel
                    f_avg += spectra[i][1]
                    lam_avg += spectra[i][0]
                    n += 1
            
            ARIEL_SPECTRA.append([0.975, f_avg/n, 0.125])
            
            #NIR Channel
            NIR = []
            for i in range(len(spectra)):
                if spectra[i][0] > 1100.0 and spectra[i][0] < 1950.0:
                    #NIR channel
                    NIR.append([spectra[i][0], spectra[i][1]])
            NIR_channel = []
            i_init = 0
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(i_init, len(NIR)):
                if NIR[i][0] - NIR[i_init][0] <= 100.0:
                    lam_avg += NIR[i][0]
                    f_avg += NIR[i][1]
                    n += 1
                else:    
                    NIR_channel.append([0.001*lam_avg/n, f_avg/n, 0.05])
                    i_init = i
                    n = 0    
                    f_avg = 0.0
                    lam_avg = 0.0
            
            ARIEL_SPECTRA.extend(NIR_channel)
            
            #Ch0 Channel
            Ch0 = []
            for i in range(len(spectra)):
                if spectra[i][0] > 1950.0 and spectra[i][0] < 3950.0:
                    Ch0.append([spectra[i][0], spectra[i][1]])
            Ch0_channel = []
            i_init = 0
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(i_init, len(Ch0)):
                if Ch0[i][0] - Ch0[i_init][0] <= 30.0:
                    lam_avg += Ch0[i][0]
                    f_avg += Ch0[i][1]
                    n += 1
                else:    
                    Ch0_channel.append([0.001*lam_avg/n, f_avg/n, 0.015])
                    i_init = i
                    n = 0    
                    f_avg = 0.0
                    lam_avg = 0.0
            
            ARIEL_SPECTRA.extend(Ch0_channel)
            
            #Ch1 Channel
            Ch1 = []
            for i in range(len(spectra)):
                if spectra[i][0] > 3950.0 and spectra[i][0] < 7850.0:
                    Ch1.append([spectra[i][0], spectra[i][1]])
            Ch1_channel = []
            i_init = 0
            f_avg = 0.0
            lam_avg = 0.0
            n = 0
            for i in range(i_init, len(Ch1)):
                if Ch1[i][0] - Ch1[i_init][0] <= 180.0:
                    lam_avg += Ch1[i][0]
                    f_avg += Ch1[i][1]
                    n += 1
                else:    
                    Ch1_channel.append([0.001*lam_avg/n, f_avg/n, 0.091])
                    i_init = i
                    n = 0    
                    f_avg = 0.0
                    lam_avg = 0.0
            
            ARIEL_SPECTRA.extend(Ch1_channel)
            
                    
            return ARIEL_SPECTRA        

        def parallel_forward(parallel_bunch):
            """
            """
            for S in parallel_bunch:
                #reinit a photometry object
                Ph = CPhotometry(S)
                #transit duration (days)
                T_dur = (S.P_planet/np.pi)*asin(sqrt((S.R_star+S.R_planet)**2 - S.bim**2)/S.A_planet)
                N_transits = int((Ph.dt2 - Ph.dt1) / S.P_planet )
                t_init_transit = (int(((Ph.dt1 - S.T0_planet)/ S.P_planet)) + 1)*S.P_planet + S.T0_planet 
                
                #iterate over all transit epochs: for each transit, compute a normalized depth curve
                for i in range(int(N_transits)):
                    #mid-transit time
                    time_stamp = t_init_transit + i*S.P_planet
                    #print("Time of transit:", time_stamp)    
                    Ph.obs_time = [time_stamp]
        
                    #no planet
                    Glob.planet = True
                    S.R_planet = 0.0
                    #print("Computing activity signal:")
                    logL, spectra_act = Ph.compute_lc_series(write_lc=False, mark='activity')
                    
                    #transits
                    Glob.planet = True
                    S.R_planet = R_planet
                    #print("Computing transits + activity:")
                    logL, spectra_act_planet = Ph.compute_lc_series(write_lc=False, mark='activity+transits_' + str(i) )
                    
                    #the first element of spectra arrays is the epoch (JD). The rest, fluxes.               
                    Depth_series = []
                    for j in range(len(Ph.wv)):
                        Depth_series.append( [float(Ph.wv[j]), 
                                            (float(spectra_act[0][j+1])-float(spectra_act_planet[0][j+1]))/float(spectra_act[0][j+1])] )         
                    
                    Depth_series_ARIEL = get_ARIEL_resolution(Depth_series)
                    io.write_file('./output/SPOTTED+PL_' + str(time_stamp) + '_' + str(S.id) + '.dat', Depth_series_ARIEL, ' ', '')
                    
                    #compute the normalized transit depth
                    Depth_activity_corrected = []
                    for j in range(len(Depth_series)):
                        Depth_activity_corrected.append( [Depth_series[j][0], Depth_series[j][1] - Depth_imm_series[j][1] ] )
                    
                    #write Depth-time array on a file
                    Depth_activity_corrected_ARIEL = get_ARIEL_resolution(Depth_activity_corrected)
                    io.write_file('./output/norm_depth_' + str(time_stamp) + '_' + str(S.id) +'.dat', Depth_activity_corrected_ARIEL, ' ', '')
                    
                
        l = mp.Lock()    
        #generate object simulation
        S = StarsimSimulation()
        
        R_planet = S.R_planet
        
        T0_pl = S.T0_planet
    
        #create a photometry object
        Ph = CPhotometry(S)
        
        #IMMACULATE PHOTOSPHERE
        print("Computing the transit events considering immaculate photosphere")
        #reset the epochs to avoid the spot map
        #Ph.dt2 = Ph.dt2 - Ph.dt1
        dt1, dt2 = Ph.dt1, Ph.dt2
        
        Ph.dt1 = 0.0 
        S.T0_planet = Ph.dt1
    
        Ph.obs_time = []
        for i in range(1,2):
            Ph.obs_time = [Ph.dt1 + i*S.P_planet]
    
        #compute lcs
        Glob.planet = True  #use grid integration method
        #immaculate photosphere without transits
        S.R_planet = 0.0
        print("Computing immaculate signal:")
        logL, spectra_imm = Ph.compute_lc_series(write_lc=False, mark='immaculate')
        
        #include transiting planet
        Glob.planet = True
        S.R_planet = R_planet
        #immaculate photosphere + planet
        print("Computing transits with immaculate photosphere:")
        logL, spectra_imm_planet = Ph.compute_lc_series(write_lc=False, mark='immaculate+transits')
    
        Depth_imm_series = []
        for j in range(len(Ph.wv)):
            Depth_imm_series.append( [float(Ph.wv[j]), 
                                     (float(spectra_imm[0][j+1])-float(spectra_imm_planet[0][j+1]))/float(spectra_imm[0][j+1])] )
        
        Depth_imm_series_ARIEL = get_ARIEL_resolution(Depth_imm_series)
        io.write_file('./IMM.dat', Depth_imm_series_ARIEL, ' ', '')
        
        S.T0_planet = T0_pl
        
        print("Immaculate depth curve computed")
        
        #reinit a photometry object
        Ph = CPhotometry(S)
        #transit duration (days)
        T_dur = (S.P_planet/np.pi)*asin(sqrt((S.R_star+S.R_planet)**2 - S.bim**2)/S.A_planet)
        N_transits = int((Ph.dt2 - Ph.dt1) / S.P_planet )
        t_init_transit = (int(((Ph.dt1 - S.T0_planet)/ S.P_planet)) + 1)*S.P_planet + S.T0_planet 
        
        #iterate over all transit epochs: for each transit, compute a normalized depth curve
        for i in range(int(N_transits)):
            #mid-transit time
            time_stamp = t_init_transit + i*S.P_planet
            #print("Time of transit:", time_stamp)    
            Ph.obs_time = [time_stamp]

            #no planet
            Glob.planet = True
            S.R_planet = 0.0
            #print("Computing activity signal:")
            logL, spectra_act = Ph.compute_lc_series(write_lc=False, mark='activity')
            
            #transits
            Glob.planet = True
            S.R_planet = R_planet
            #print("Computing transits + activity:")
            logL, spectra_act_planet = Ph.compute_lc_series(write_lc=False, mark='activity+transits_' + str(i) )
            
            #the first element of spectra arrays is the epoch (JD). The rest, fluxes.               
            Depth_series = []
            for j in range(len(Ph.wv)):
                Depth_series.append( [float(Ph.wv[j]), 
                                    (float(spectra_act[0][j+1])-float(spectra_act_planet[0][j+1]))/float(spectra_act[0][j+1])] )         
            
            Depth_series_ARIEL = get_ARIEL_resolution(Depth_series)
            io.write_file('./output/SPOTTED+PL_' + str(time_stamp) + '_' + str(S.id) + '.dat', Depth_series_ARIEL, ' ', '')
            
            #compute the normalized transit depth
            Depth_activity_corrected = []
            for j in range(len(Depth_series)):
                Depth_activity_corrected.append( [Depth_series[j][0], Depth_series[j][1] - Depth_imm_series[j][1] ] )
            
            #write Depth-time array on a file
            Depth_activity_corrected_ARIEL = get_ARIEL_resolution(Depth_activity_corrected)
            io.write_file('./output/norm_depth_' + str(time_stamp) + '_' + str(S.id) +'.dat', Depth_activity_corrected_ARIEL, ' ', '')
                
        
        
        
        
        
        
        
        
        
        
        
        sys.exit()
        
        maps = []
        maps_all = []
        maps_all.extend(glob.glob('./good_5000/spot*') )
        
        #read the param. file
        param_data = io.read_file('./good_5000/good_sols_5000.dat', ' ')
        
        #apply constraints in logL
        #read the list with logL
        results = io.read_file('./good_5000/good_sols_5000.dat', ' ')
        maps_refs = []
        for j in range(len(results)):
            if float(results[j][1]) >= 0.0:
                maps_refs.append(results[j][0])
                for map_all in maps_all:
                    ref_code = str((map_all.split('/')[-1]).split('_')[2].split('.')[0])
                    if ref_code == results[j][0]:
                        maps.append(map_all)
        
        print ("Total maps:", len(maps))
        
        #check lengths
        if len(maps) != len(param_data):
            print ("Warning: the number of maps is not equal to the parameter list")
            
        print(len(maps), "maps read...")
        print(len(param_data), "parameter lines read...")
        
        #generate an initial population of StarsimSimulation objects
        INIT_POPULATION = len(maps)
        io.print_message("Generating initial population: " + str(INIT_POPULATION) + " simulations", 2, 94)
        
        #inhibit the screen messages, for cleanness
        Glob.inhibit_msg = True
        
        Simulations = [StarsimSimulation() for i in range(INIT_POPULATION)]
    
        #prepare the objects
        for j,S in enumerate(Simulations):
            #load the data and specific parameters for this data, in each object of the population
            OB = S.read_multifit_config('./multifit.conf')
            
            #read spot map
            spot_lst = S.read_file(maps[j], ' ')
            #copy the read map to the current simulation object
            S.spot_map[0] = spot_lst
            
            #look for the assoc. parameters to this map
            ref_code = str((maps[j].split('/')[-1]).split('_')[2].split('.')[0])
            
            for i in range(len(param_data)):
                if ref_code == str(param_data[i][0]):
                    param_lst = param_data[i]
            
            S.id = ref_code
            #load the parameters linked with the map
            S.P_rot = float(param_lst[2])
            # set Q ratio
            S.set_Q(float(param_lst[3]))
            #delta T spots
            delta_t_sp = float(param_lst[4])
            l.acquire()
            S.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=S.delta_t_fc)
            l.release()
            """
            for i in range(len(OB)):
                if S.Observables[i].data_type == 'ph':
                    # assign a new jitter to reference series
                    S.Observables[i].jitter = float(param_lst[5])
            """
        
        #run parallel
        pool = mp.Pool(Glob.ncpus)
        parallel_bunches = S.split_list(Simulations, Glob.ncpus)
        
        print("Running parallel instances")
        
        try:
            parallel_data_blocks = pool.map_async(parallel_forward, parallel_bunches).get(9999999)
        except KeyboardInterrupt:
            io.print_message('Parallel pool killed. Bye!', index=3, color=31)
            pool.terminate()
            sys.exit()
        #free memory
        pool.terminate()
        
    elif mode == 'testing':
        """
        """
        Simulation = StarsimSimulation()
       
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
        except:
            print("Could not be read the multiobject configuration or no data selected")
            pass
        try:
            #modify stellar params and lc offset
            for iter in range(10000):
                for object in OB:
                    if object.data_type == 'ph':
                        logL_ph = object.compute_lc_series(write_lc=False)
        except:
            Ph = CPhotometry(Simulation)
            for iter in range(100000):
                Ph.compute_lc_series(write_lc=True)
        
    elif mode == 'mcmc':
        """
        """
        def find_nearest(array, value):
            idx = abs(np.abs(array)-value).argmin()
            return array[idx]
    
        def forward_photometric_model(theta):
            """
            """
            p_rot, Q, DT, s = theta
            
            joint_logL = 0.0
            
            Simulation.P_rot = p_rot
            Simulation.Q = Q
            delta_t_sp_ = DT
            delta_t_sp = find_nearest([50.0 + 25.0*j for j in range(int((2000.0-50.0)/25.0))], delta_t_sp_)
            Simulation.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=Simulation.delta_t_fc)
            jitter = s
            
            for ph in OB:
                ph.jitter = jitter
                logL_ph = ph.compute_lc_series(write_lc=False)
                joint_logL += logL_ph[0]
            
            return joint_logL
        
        def log_prior(theta):
            """
            """
            p_rot, Q, DT, s = theta
            
            if 17.0 < p_rot < 18.0 and 0.0 < Q < 3.0 and 350.0 < DT < 1000.0 and 0.0 < s < 0.01:
                
                return 0.0
            
            return -np.inf

        def log_probability(theta):
            """
            """
            lp = log_prior(theta)
            
            if not np.isfinite(lp):
                return -np.inf
        
            return lp + forward_photometric_model(theta)

        def MH_simple(prob, x0, n):
            """
            """
            x_prev = x0
            x = []
            i = 0
            k = 0
            while i < n:
                x_star = [random.uniform(15.0, 19.0),random.uniform(0.0,3.0),random.uniform(50.0,2000.0),random.uniform(0.0,0.01)]
                p_star = prob(x_star)
                p_prev = prob(x_prev)
                
                U = np.random.uniform()
                
                #acceptance ratio
                A = p_star / p_prev
            
                if U < A:
                    x.append(x_star)
                    x_prev = x_star
                    i += 1
                    k += 1
                else:
                    #x_prev = x[i]
                    i += 1
                    
            return x
        
        def metropolis(prob, n):
            """
            """
            x_0 = [[random.uniform(15.0, 19.0),random.uniform(0.0,3.0),
                    random.uniform(50.0,2000.0),random.uniform(0.0,0.01)] for i in range(n) ]
            
            for i in range(10):
                for j in range(len(x_0)):
                    x_prev = x_0[j]
                    x_star = [random.uniform(15.0, 19.0),random.uniform(0.0,3.0),random.uniform(50.0,2000.0),random.uniform(0.0,0.01)]
                    
                    p_star = prob(x_star)
                    p_prev = prob(x_prev)
                    
                    U = np.random.uniform()
                    
                    #acceptance ratio
                    A = p_star - p_prev
                    if U < A:
                        x_0[j] = x_star
                    
                
            return x_0
        
        Simulation = StarsimSimulation()
      
        #read data and configuration
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
        except:
            print("Could not be read the multiobject configuration or no data selected")
            raise
            sys.exit()
        
        Glob.inhibit_msg = True
        
        
        #samples = metropolis(forward_photometric_model, 50)
        #io.write_file('samples.dat', samples, ' ', '')
        
        import emcee
        from multiprocessing import Pool

        nwalkers = 64
        
        p_0 = [17.75, 0.0, 825.0, 0.0056]
        ndim = len(p_0)
        
        pos = p_0 + [1.0, 1.0, 150.0, 0.0008]* np.random.randn(nwalkers, ndim)
        nwalkers, ndim = pos.shape
  
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)
        sampler.run_mcmc(pos, 2000, progress=True)
      
        flat_samples = sampler.get_chain(discard=150, thin=2, flat=True)
        print len(flat_samples)
        
        io.write_file('samples.dat', flat_samples, ' ', '')
    
    elif mode == 'inverse_nn':
        """
        """
        def find_nearest(array, value):
            idx = abs(np.abs(array)-value).argmin()
            return array[idx]
        
        def flatten(lis):
            """return a nested list, flattened
            """
            f_lis = []
            for item in lis:
                if type(item) == type([]):
                    f_lis.extend(flatten(item))
                else:
                    f_lis.append(item)
            
            return f_lis
        
        def spotmap_reshape(serial_spotmap):
            """
            """
            spot = []
            spotmap_reshaped = []
            for i in range(len(serial_spotmap)):
                spot.append(serial_spotmap[i] * norm[i % 5])
                if len(spot) == 5:
                    spotmap_reshaped.append(spot)
                    spot = []
                
            return spotmap_reshaped
        
        def gen_dataset(n):
            """
            """
            def randomize_spotmap(spotmap):
                """
                """
                spotmap_ = copy.deepcopy(spotmap)  #tmp map
                
                while(1):
                    for j in range(1):
                        spot = int(random.uniform(0,len(spotmap_)))
                        parameter = int(random.uniform(2,5))
                        #modify coordinates and sizes
                        if parameter == 2: #colatitude
                            spotmap_[0][spot][parameter] = random.uniform(0.0,180.0)
                        elif parameter == 3: #longitude
                            spotmap_[0][spot][parameter] = random.uniform(0.0,360.0)
                        elif parameter == 4: #radii
                            spotmap_[0][spot][parameter] = random.uniform(0.0,10.0)
                            #spotmap_[0][spot][parameter] = abs(spotmap_[0][spot][parameter])
                    if Simulation.check_params_compatibility(spotmap_):
                        Simulation.spot_map = spotmap_
                        break
                    
                #serialize map
                serial_spotmap = []
                
                for i in range(len(spotmap_[0])):
                    for j in range(0,5):
                        serial_spotmap.append(spotmap_[0][i][j] / norm[j])
         
                return serial_spotmap 
                    
            Simulation = StarsimSimulation()
        
            OB = []
            try:
                #list of objects to joint minimize 
                OB = Simulation.read_multifit_config('./multifit.conf')
            except:
                print("Could not be read the multiobject configuration or no data selected")
            
            LCS = []
            params = []
            
            Glob.inhibit_msg = True
            for i in range(int(n)):
                #randomize spotmap
                serial_spotmap = randomize_spotmap(Simulation.spot_map)
                #change spot temperature
                delta_t_sp_ = random.uniform(100.0, 2000.0)
                delta_t_sp = find_nearest([50.0 + 25.0*j for j in range(int((2000.0-50.0)/25.0))], delta_t_sp_)
                #lock
                l.acquire()
                Simulation.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=Simulation.delta_t_fc)
                l.release()
                period = random.uniform(12.5,17.5)
                Simulation.P_rot = period
                lcs = []
                if len(OB) > 0:
                    for ph in OB:
                        #is photometry object?
                        if ph.data_type == 'ph':     
                            logL_ph = ph.compute_lc_series(write_lc=False)
                            lcs.extend(ph.mean_normalized_series[:,1] - 1.0)
                            
                p = []
                
                #mean = np.mean(ph.ff_sp[:,1])
                #for j in range(len(ph.ff_sp[:,1])):
                #    p.append(ph.ff_sp[j,1] - mean)
                
                #for j in range(len(serial_spotmap)):
                #    p.append(serial_spotmap[j])
            
                #for j in range(len(lcs)):
                #    p.append(lcs[j])
                
                p.append((delta_t_sp-1000.)/1000.0)
                p.append((period-14.5)/14.5)
                params.append(p)
                LCS.append(lcs)
                
            t = ph.mean_normalized_series[:,0]    
        
            X = np.array(LCS)
            Y = np.array(params)
           
            
            return X, Y, t
        
        def bootstrapped_set(x, sigma, N):
            """
            """
            L = []
            for iter in range(N):
                X_ = []
                for i in range(len(x)):
                    X_.extend(x[i] + random.gauss(0.0, sigma))
                L.append(np.array(X_))
                
            return np.array(L)
            
        def gen_single_set(p_rot, temp):
            """
            """
            Simulation = StarsimSimulation()
        
            OB = []
            try:
                #list of objects to joint minimize 
                OB = Simulation.read_multifit_config('./multifit.conf')
            except:
                print("Could not be read the multiobject configuration or no data selected")
            
            
            LCS = []
            params = []
            
            Glob.inhibit_msg = True
            for i in range(1):
                delta_t_sp_ = temp
                delta_t_sp = find_nearest([50.0 + 25.0*j for j in range(int((2000.0-50.0)/25.0))], delta_t_sp_)
                #lock
                l.acquire()
                Simulation.update_spectral_data(delta_t_sp=delta_t_sp, delta_t_fc=Simulation.delta_t_fc)
                l.release()
                Simulation.P_rot = p_rot
                lcs = []
                if len(OB) > 0:
                    for ph in OB:
                        #is photometry object?
                        if ph.data_type == 'ph':     
                            logL_ph = ph.compute_lc_series(write_lc=False)
                            lcs.extend(ph.mean_normalized_series[:,1] - 1.0)
                            
                LCS.append(lcs)
                p = []
                mean = np.mean(ph.ff_sp[:,1])
                for j in range(len(ph.ff_sp[:,1])):
                    p.append(ph.ff_sp[j,1] - mean)
                p.append((delta_t_sp-1000.0)/1000.0)
                p.append((p_rot-14.5)/14.5)
                params.append(p)
                
            t = ph.mean_normalized_series[:,0]    
        
            X = np.array(LCS)
            Y = np.array(params)
            
            return X, Y, t
        
        
        Simulation = StarsimSimulation()
        #spot parameters normalization
        norm = [abs(Simulation.dt1 - 2.*Simulation.spots_lifetime), 3.*Simulation.spots_lifetime, 180.0, 360.0, 10.0]
        
        OB = []
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
        except:
            print("Could not be read the multiobject configuration or no data selected")
            
        
        import tensorflow as tf
        from tensorflow import keras
        from tensorflow.keras.layers import LeakyReLU

        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' #info messages are not printed

        l = mp.Lock()
    
        train = True
    
        if train:
            
            X, Y, t = gen_dataset(100)
                
            x_train, y_train = np.array(X[:int(0.9*len(X))]), np.array(Y[:int(0.9*len(Y))])
            x_test, y_test = np.array(X[int(0.9*len(X)):]), np.array(Y[int(0.9*len(Y)):])
            
            x_train_1 = []
            x_train_2 = []
            x_train_3 = []
            for i in range(len(x_train)):
                x_train_1.append(x_train[i][:len(x_train[0])/3])
                x_train_2.append(x_train[i][len(x_train[0])/3:2*len(x_train[0])/3])
                x_train_3.append(x_train[i][2*len(x_train[0])/3:])
            
            x_test_1 = []
            x_test_2 = []
            x_test_3 = []    
            for i in range(len(x_test)):
                x_test_1.append(x_test[i][:len(x_test[0])/3])
                x_test_2.append(x_test[i][len(x_test[0])/3:2*len(x_test[0])/3])
                x_test_3.append(x_test[i][2*len(x_test[0])/3:])
            
            x_train_1 = np.array(x_train_1)
            x_train_2 = np.array(x_train_2)
            x_train_3 = np.array(x_train_3)
            x_test_1 = np.array(x_test_1)
            x_test_2 = np.array(x_test_2)
            x_test_3 = np.array(x_test_3)
            
            print x_train.shape
            print y_train.shape
            print x_train_1.shape
            print x_test_1.shape
            #encoding_dim = y_train.shape[1]   # n floats 
            
            # this is our input placeholder
            input_1 = tf.keras.layers.Input(shape=(x_train.shape[0]/3,x_train.shape[1],1))
            input_2 = tf.keras.layers.Input(shape=(x_train.shape[0]/3,x_train.shape[1],1))
            input_3 = tf.keras.layers.Input(shape=(x_train.shape[0]/3,x_train.shape[1],1))
            
            x_train_1 = x_train_1.reshape(x_train_1.shape[1], x_train_1.shape[0]).astype('float32')
            x_test_1 = x_test_1.reshape(x_test_1.shape[1], x_test_1.shape[0]).astype('float32')
            
            x = tf.keras.layers.Conv1D(filters=64, kernel_size=3, activation='relu')(input_1)
            x = tf.keras.models.Model(inputs=input_1, outputs=x)
            y = tf.keras.layers.Conv1D(filters=64, kernel_size=3, activation='relu',input_shape=(x_train_2.shape[1],x_train_2.shape[0]))(input_2)
            y = tf.keras.models.Model(inputs=input_2, outputs=y)
            z = tf.keras.layers.Conv1D(filters=64, kernel_size=3, activation='relu',input_shape=(x_train_3.shape[1],x_train_3.shape[0]))(input_3)
            z = tf.keras.models.Model(inputs=input_3, outputs=z)
            """
            # the first branch operates on the first input
            x = tf.keras.layers.Dense(64, activation='tanh')(input_1)
            x = tf.keras.layers.Dense(32, activation='tanh')(x)
            x = tf.keras.layers.Dense(2, activation='tanh')(x)
            x = tf.keras.models.Model(inputs=input_1, outputs=x)
            # the second branch opreates on the second input
            y = tf.keras.layers.Dense(64, activation='tanh')(input_2)
            y = tf.keras.layers.Dense(32, activation='tanh')(y)
            y = tf.keras.layers.Dense(2, activation='tanh')(y)
            y = tf.keras.models.Model(inputs=input_2, outputs=y)
            # the third branch opreates on the third input
            z = tf.keras.layers.Dense(64, activation='tanh')(input_3)
            z = tf.keras.layers.Dense(32, activation='tanh')(z)
            z = tf.keras.layers.Dense(2, activation='tanh')(z)
            z = tf.keras.models.Model(inputs=input_3, outputs=z)
            """
            # combine the output of the two branches
            combined = tf.keras.layers.concatenate([x.output, y.output, z.output])
            
            # apply a FC layer and then a regression prediction on the
            # combined outputs
            c = tf.keras.layers.Dense(128, activation='tanh')(combined)
            c = tf.keras.layers.Dense(y_train.shape[1], activation='tanh')(c)
            
            # our model will accept the inputs of the two branches and
            # then output a single value
            model = tf.keras.models.Model(inputs=[x.input, y.input, z.input], outputs=c)
            
            model.compile(optimizer='adam', loss='mse')
        
            
            model.fit([x_train_1, x_train_2, x_train_3], y_train,   
                        epochs=1500,
                        batch_size=128,
                        shuffle=True,
                        validation_data=([x_test_1, x_test_2, x_test_3], y_test))
            
            model_data = model.predict([x_test_1,x_test_2,x_test_3])
            
            periods = []
            temperatures = []
            for j in range(len(model_data)):
                temperatures.append([1000. + 1000.*y_test[j][-2], 1000. + 1000.*model_data[j][-2]])
                periods.append([14.5+14.5*y_test[j][-1], 14.5 + 14.5*model_data[j][-1]])
                
            io.write_file('periods.dat', periods, ' ', '')
            io.write_file('temperatures.dat', temperatures, ' ', '')
            
            model.save('encoder_1K.h5', overwrite=True, include_optimizer=True)
            
        if not train:
            model = keras.models.load_model('encoder_1K.h5', compile=True)
            
        """
        X_obs_data, Y_obs_data, t = gen_single_set(p_rot=13.0, temp=700.0)
        
        bootstrappend_set_ = bootstrapped_set(X_obs_data, 0.001, 10)
        predicted = model.predict(bootstrappend_set_) 
        
        for i in range(len(predicted)):
            print 1000. + 1000.*predicted[i][-2], 14.5 + 14.5*predicted[i][-1]
        
        decoded_data = decoder.predict(y_test)
        
        series = []
        for j in range(len(decoded_data[0])/3):
            series.append([t[j], x_test[0][j], x_test[1][j], x_test[2][j], x_test[3][j],
                           decoded_data[0][j], decoded_data[1][j], decoded_data[2][j], decoded_data[3][j]])
           
            
        io.write_file('reconstructed.dat', series, ' ', '') 
        
        """
        #print spotmap_reshape(x_test[0])
        #print spotmap_reshape(decoded_imgs[0])
        """
        for i in range(len(x_test)):
            true_map = spotmap_reshape(x_test[i][0:25])
            Simulation.spot_map[0] = true_map
            Simulation.update_spectral_data(delta_t_sp=2000.*x_test[i][-2], delta_t_fc=Simulation.delta_t_fc)
            Simulation.P_rot = 17.5*x_test[i][-1]
            logL_true = -Simulation.compute_joint_statistic(OB)
            print "True:", Simulation.delta_t_sp, Simulation.P_rot, logL_true
            decoded_map = spotmap_reshape(decoded_maps[i][0:25])
            Simulation.spot_map[0] = decoded_map
            Simulation.update_spectral_data(delta_t_sp=2000.*decoded_maps[i][-2], delta_t_fc=Simulation.delta_t_fc)
            Simulation.P_rot = 17.5*decoded_maps[i][-1]
            logL_decoded = -Simulation.compute_joint_statistic(OB)
            print "Rec.:", Simulation.delta_t_sp, Simulation.P_rot, logL_decoded
            
        
        #print encoded_imgs[0]
        """
        """

        out_dim = len(Y[0])
        #NN encoded --> params
        encoded_input = tf.keras.layers.Input(shape=(encoding_dim,))
        layers = tf.keras.layers.Dense(encoding_dim, activation='tanh')(encoded_input)
        layers = tf.keras.layers.Dense(50, activation='tanh')(layers)
        layers = tf.keras.layers.Dense(50, activation='tanh')(layers)
        layers = tf.keras.layers.Dense(50, activation='tanh')(layers)
        layers = tf.keras.layers.Dense(50, activation='tanh')(layers)
        layers = tf.keras.layers.Dense(50, activation='tanh')(layers)
        layers = tf.keras.layers.Dense(50, activation='tanh')(layers)
        #layers = tf.keras.layers.Dense(500, activation='tanh')(layers)
        #layers = tf.keras.layers.Dense(56, activation='tanh')(layers)
        layers = tf.keras.layers.Dense(out_dim, activation='relu')(layers)
        
        #model
        model = tf.keras.models.Model(encoded_input, layers)
        
        #compile all
        model.compile(optimizer='adam', loss='mse')
        
        encoded_train = encoder.predict(x_train)
        encoded_test = encoded_imgs
        
        print encoded_train
        #fit training data
        model.fit(encoded_train, y_train,
                    epochs=20,
                    batch_size=32,
                    shuffle=True,
                    validation_data=(encoded_test, y_test))
        
        y_predicted = model.predict(encoded_imgs)
        
        y_pred_reshaped = encoded_imgs.reshape(5,3)
        
        print y_pred_reshaped 
        
        
            
        series_ = []
        for j in range(len(y_predicted[0])-2):
            series_.append([t[j], y_test[0][j], y_test[1][j], y_test[2][j], y_test[3][j],
                           y_predicted[0][j], y_predicted[1][j], y_predicted[2][j], y_predicted[3][j]])
           
            
        io.write_file('ff.dat', series_, ' ', '') 
            
        
        import matplotlib.pyplot as plt

        n = 10  # how many digits we will display
        plt.figure(figsize=(20, 5))
        
        for i in range(n):
            # display reconstruction
            ax = plt.subplot(2, n, i + 1 + n)
            plt.imshow(encoded_imgs[i].reshape(9, 9))
            plt.gray()
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        
            
        plt.show()
        """
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        """
        train = True

        if train:
            
            X, Y, t = gen_dataset(15000)
            
            x_train, y_train = np.array(X[:int(0.9*len(X))]), np.array(Y[:int(0.9*len(Y))])
            x_test, y_test = np.array(X[int(0.9*len(X)):]), np.array(Y[int(0.9*len(Y)):])
            
            print x_train.shape
            print y_train.shape
            #print x_train[0]
            print len(y_train[0])
            
            io.write_file('X.dat', X, ' ', '')
            io.write_file('Y.dat', Y, ' ', '')
            
            out_dim = len(Y[0])
        
            #this is our input placeholder
            input_signal = tf.keras.layers.Input(shape=(x_train.shape[1],))
            #layer block
            layers = tf.keras.layers.Dense(150, activation='tanh')(input_signal)
            layers = tf.keras.layers.Dense(150, activation='tanh')(layers)
            layers = tf.keras.layers.Dense(150, activation='tanh')(layers)
            layers = tf.keras.layers.Dense(150, activation='tanh')(layers)
            layers = tf.keras.layers.Dense(150, activation='tanh')(layers)
            layers = tf.keras.layers.Dense(150, activation='tanh')(layers)
            layers = tf.keras.layers.Dense(150, activation='tanh')(layers)
            layers = tf.keras.layers.Dense(150, activation='tanh')(layers)
            layers = tf.keras.layers.Dense(150, activation='tanh')(layers)
            layers = tf.keras.layers.Dense(150, activation='tanh')(layers)
            layers = tf.keras.layers.Dense(out_dim, activation='linear')(layers)
            #model
            model = tf.keras.models.Model(input_signal, layers)
            
            #compile all
            model.compile(optimizer='adam', loss='mse')
            
            #fit training data
            model.fit(x_train, y_train,
                        epochs=1000,
                        batch_size=32,
                        shuffle=True,
                        validation_data=(x_test, y_test))
        
            #save model training and config. 
            model.save('UJ_map.h5', overwrite=True, include_optimizer=True)
        
        
        if not train:
            model = keras.models.load_model('UJ_map.h5', compile=True)
        
        
        X, Y, t = gen_dataset(150)
        x_train, y_train = np.array(X[:int(0.9*len(X))]), np.array(Y[:int(0.9*len(Y))])
        x_test, y_test = np.array(X[int(0.9*len(X)):]), np.array(Y[int(0.9*len(Y)):])
        
        #refit
        #fit training data
        model.fit(x_train, y_train,
                    epochs=125,
                    batch_size=32,
                    shuffle=True,
                    validation_data=(x_test, y_test))
        
        #X, Y, t = gen_dataset(100)
        #x_test, y_test = np.array(X[:]), np.array(Y[:])
            
        #y_predicted = model.predict(x_test)
        
        #bootstrapping
        x, y, t = gen_single_set(15.0, 700.0)
        X_boot = bootstrapped_set(x, 0.05, 50000)

        
        y_boot = model.predict(X_boot)
        
        mcmc_2 = []
        for j in range(len(y_boot)):
            #print y_boot[j][-1], 2000.0*y_boot[j][-2]
            mcmc_2.append([y_boot[j][-1], 2000.0*y_boot[j][-2]])
        
        io.write_file('mcmc_2.dat', mcmc_2, ' ', '') 
        
        y_predicted = model.predict(x_test)
        
        for j in range(len(y_predicted)):
            print y_test[j][-2], y_predicted[j][-2]
            print y_test[j][-1], y_predicted[j][-1]
        
        series = []
        for j in range(len(y_predicted[0])-2):
            series.append([t[j], y_test[0][j], y_test[1][j], y_predicted[0][j], y_predicted[1][j]])
           
            
        io.write_file('ff.dat', series, ' ', '') 
        """
        
    elif num_spots > 0:
        """generates a random spot map with a given number of spots.
        """
        #USAGE: random spot map ( --spotmap=## )
        #read science data, in case file passed in command line 
        #create an object simulation that encloses all the simulation parameters
        simulation = StarsimSimulation()
        try:
            #photometry
            science_data_ph = simulation.get_science_data(sys.argv[1])  
            time_array_lc, lc_data = science_data_ph[0], science_data_ph[1]
            #time limits
            simulation.dt1 = float(time_array_lc[0])
            simulation.dt2 = float(time_array_lc[-1]) 
            
        except:
            #if no lc, use the time limits defined in starsim.conf
            lc_data = list()
            time_array_lc = time_array_lc = np.linspace(simulation.dt1, simulation.dt2,  
                            num=int(1440.0*(simulation.dt2-simulation.dt1)/simulation.t_exp))

        spot_lst = simulation.gen_rnd_map(num_spots)
        simulation.write_file('./data/spot_list.dat', spot_lst, ' ', '')
        #check no collisions among spots
        simulation.init_spotmap()
        while not simulation.check_params_compatibility(simulation.spot_map):
            spot_lst = simulation.gen_rnd_map(num_spots)
            simulation.write_file('./data/spot_list.dat', spot_lst, ' ', '')
            #check no collisions among spots
            Glob.inhibit_msg = True
            simulation.init_spotmap()
    
    else:
        print("Enter a valid mode --mode=[ph/rv/forward/inversion/N-inversion/param_fit/N-depth/static_N-depth...]")
        print('') 
                
        








