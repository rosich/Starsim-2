#!/usr/bin/python
#-----------------------------------------------------------------------------------------------------------------------------------
#Name:      starsim_2.py
#
#Purpose:   Python implementation of E.Herrero's starsim_v11 photosphere simulator to allow parallelism and inverse problem.
#           StarSim-2 provides simulation of spots and faculae effects on both photometric and RV time series due to star rotation.
#
#Author:    Albert Rosich (rosich@ice.cat).
#
#Created:   2014/05/04
#
#Last Update: 2018/08/27

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
import time
from math import sin, cos, acos, sqrt, pi, log, log10, exp
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


def print_start_info():
    """
    """
    io.print_message('', 7,94)
    io.print_message('                     STARSIM 2.0                       ', 7, 90)
    #self.print_message('                                                       ', 7, 94)
    io.print_message('     SIMULATION OF STELLAR PHOTOMETRIC / RV ACTIVITY   ', 7, 94)                                     
    #self.print_message('                                                       ', 7, 94)
    io.print_message('                                                       ', 7, 90)
    io.print_message('', 7,90)
    #io.print_message("Read configuration file " + str(self.conf_file_path), index=2, color=31)
    #self.print_message( 'Filter TF data file: ' + str(conf_file.get('files','filter_path')), 3, 37)
    #self.print_message( '', 7,37)
    #self.print_message( 'Number of pixels on the star: ' + str(n1), 3, 37)
    #self.print_message( 'Size of the pixel: ' + str(da1) + 'deg', 3, 37)
    #self.print_message( 'Radius of the star: ' + str(star.R_star) + 'Rsun', 3, 37)
    if Glob.ncpus != 0:
        io.print_message( 'Detected CPUs / using CPUs: ' + str(mp.cpu_count()) + "/" + str(Glob.ncpus), 5, 95)
        io.print_message('', 5,95)
        io.print_message('INIT. STARSIM...', 6,90)
        
def split_list(alist, wanted_parts):
    """
    """
    length = len(alist)
    if length < wanted_parts:
        wanted_parts = length
        
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
            for i in range(wanted_parts) ]
  
        
if __name__ == '__main__':
    """
    """
#---------------------------------------------------------------------------#
# COMMON PART: ALL VARIABLES DECLARED HERE ARE ACCESSIBLE IN ALL SCOPES     #
#---------------------------------------------------------------------------#
    
    io = Starsim_IO() #makes available IO functions
    #option capture
    options, remainder = getopt.gnu_getopt(sys.argv[1:], 'c:i:m:s:l:xcf', ['planet','ncpus=','inhibit_msg','mode=','spotmap=', 'label=', 'indices', 'bw', 'save_ccfs'])
    Glob.ncpus = mp.cpu_count()
    mode = None
    label = ''
  
    #argument parsing
    for opt, arg in options:
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
    
    #print some init info
    print_start_info()
#---------------------------------------------------------------------------------------------------------------#        
#                                            SELECT OPTION                                                      #
#                                            -------------                                                      #
#                                                                                                               #   
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
            print "Could not be read the multiobject configuration or no data selected"
        
        if len(OB) > 0 and not Glob.planet:
            for ph in OB:
                #is photometry object?
                if ph.data_type == 'ph':     
                    ph.set_ph_zero_point(ph.z_ph)
                    #write lightcurve file
                    logL_ph = ph.compute_lc_series(write_lc=True)
                    print "\nStatistics:"
                    print "\tNull-model logL:", ph.logL_null 
                    print "\tModel logL:", logL_ph
                    print "\tDlogL(model-data):", logL_ph - ph.logL_null
                    print ""
                    #generate data curve with z_ph implemented + jitter
                    ph.write_lc_data()
                    #write data series + residuals (last col)
                    ph.compute_lc_series(write_lc=True)
                    ph.obs_time = np.linspace(ph.dt1-10.0, ph.dt2+10.0, num=int(1440.0*(ph.dt2-ph.dt1)/ph.t_exp))
                    ph.obs_data = []
                    #compute the HR lc curve
                    ph.compute_lc_series(write_lc=True, mark='HR')
            
        else:
            #no data files in multifit.conf
            print "No photometry objects in ./multifit.conf"
            Ph = CPhotometry(Simulation)
            #compute lc
            Ph.compute_lc_series(write_lc=True)
                                                 
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
            print "Could not be read the multiobject configuration or no data selected"
  
        if len(OB) > 0:
            for rv in OB:
                if rv.data_type == 'rv':
                    print "Going to use the multifit.conf data: ", rv.science_data 
                    #compute rv
                    logL_rv = rv.compute_rv_series(write_rv=True)
                    print "\nStatistics:"
                    print "\tNull-model logL:", rv.logL_null 
                    print "\tModel logL:", logL_rv
                    print "\tDlogL(model-data):", logL_rv - rv.logL_null
                    print ""
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
            print "Could not be read the multiobject configuration or no data selected"
  
        if len(OB) > 0:
            for bis in OB:
                if bis.data_type == 'bis':
                    print "Going to use the multifit.conf data: ", bis.science_data 
                    #compute rv
                    logL_bis = bis.compute_BIS_series(write_rv=True)
                    print "\nStatistics:"
                    #print "\tNull-model logL:", bis.logL_null 
                    #print "\tModel logL:", logL_bis
                    #print "\tDlogL(model-data):", logL_bis - bis.logL_null
                    print ""
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
            print "Could not be read the multiobject configuration or no data selected"
  
        if len(OB) > 0:
            for fwhm in OB:
                if fwhm.data_type == 'fwhm':
                    print "Going to use the multifit.conf data: ", fwhm.science_data 
                    #compute rv
                    logL_fwhm = fwhm.compute_FWHM_series(write_rv=True)
                    print "\nStatistics:"
                    print "\tNull-model logL:", fwhm.logL_null 
                    print "\tModel logL:", logL_fwhm
                    print "\tDlogL(model-data):", logL_fwhm - fwhm.logL_null
                    print ""
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
            print "Could not be read the multiobject configuration or no data selected"
  
        if len(OB) > 0:
            for contrast in OB:
                if contrast.data_type == 'contrast':
                    print "Going to use the multifit.conf data: ", contrast.science_data 
                    #compute rv
                    logL_contrast = contrast.compute_CONTRAST_series(write_rv=True)
                    print "\nStatistics:"
                    print "\tNull-model logL:", contrast.logL_null 
                    print "\tModel logL:", logL_contrast
                    print "\tDlogL(model-data):", logL_contrast - contrast.logL_null
                    print ""
        else:
            #no data files in multifit.conf
            contrast = CCONTRAST(Simulation)
            #compute rv
            contrast.compute_CONTRAST_series(write_rv=True) 
    
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
            print "Could not be read the multiobject configuration or no data selected"
            raise
        #modify stellar params and lc offset
        for object in OB:
            if object.data_type == 'ph':
                object.set_ph_zero_point(object.z_ph)
                
        #fit observables
        res = Simulation.simulated_annealing(Simulation.Observables)
        print ""

        #copy the inverted map to the simulation environment 
        Simulation.spot_map[0] = Simulation.InvertedSpotMap 
                
        #compute the forward problem (HR)
        #generate curves
        if len(OB) > 0:
            for ob in OB:
                if ob.data_type == 'rv':
                    print "\nRV object"
                    logL_rv = ob.compute_rv_series(write_rv=False)
                    print "\nStatistics:"
                    print "\tNull-model logL:", ob.logL_null 
                    print "\tModel logL:", logL_rv
                    print "\tDlogL(model-data):", logL_rv - ob.logL_null
                    print ""
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
                    print "\nPhotometry object"
                    logL_ph = ob.compute_lc_series(write_lc=False)
                    print "\nStatistics:"
                    print "\tNull-model logL:", ob.logL_null
                    print "\tModel logL:", logL_ph
                    print "\tDlogL(model-data):", logL_ph - ob.logL_null
                    print ""
                    #generate data curve with z_ph implemented + jitter
                    ob.write_lc_data()
                    #write data series + residuals (last col)
                    ob.compute_lc_series(write_lc=True)
                    ob.obs_time = np.linspace(ob.dt1-10.0, ob.dt2+10.0, num=int(1440.0*(ob.dt2-ob.dt1)/ob.t_exp))
                    ob.obs_data = []
                    #compute the HR lc curve
                    ob.compute_lc_series(write_lc=True, mark='HR')
                                     
    elif mode == 'N-inversion':
        """fit multiple identic objects
        """        
        def sigma_clip(a_array):
            """
            """
            from scipy.stats import sigmaclip
            
            for i in range(a_array):
                points = []
                for j in range(a_array[0]):
                    points.append(a_array[i][j])
                #sigma clipping
                c, low, upp = sigmaclip(points, 2.0, 2.0)
                print c
        
        def model_interpolation(time_data, y_data, time_model, y_model):
            """model,data interploation of arbitrarily sampled model
            """
            y_interpolated_model = []
            for j in range(len(y_data)):
                for i in range(len(y_model)-1):
                    if (time_model[i] <= time_data[j]) and (time_model[i+1] >= time_data[j]):
                        interp_y = y_model[i] + ((y_model[i+1]-y_model[i])/(time_model[i+1]-time_model[i]))*(time_data[j]-time_model[i])
                        y_interpolated_model.append(interp_y)

            return np.array(y_interpolated_model)
                   
        def compute_joint_fit(parallel_bunch):
            """receives a list of Simulation objects, inverts observational data 
               and computes forward problem
            """        
            def find_nearest(array,value):
                idx = abs(np.abs(array)-value).argmin()
                return array[idx]
              
            Glob.inhibit_msg = True
         
            bunch_data = []
            for S in parallel_bunch:
                for i in range(len(S.Observables)):
                    if S.Observables[i].data_type == 'ph':
                        S.Observables[i].set_ph_zero_point(S.Observables[i].z_ph)
                    #recompute CCF for the new rotation period
                    elif S.Observables[i].data_type == 'rv':
                        S.Observables[i].init_ccfs()
                
                #model inversion   
                joint_logL, id, info, RV_WEIGHT = S.simulated_annealing(S.Observables)
                print "process:", S.id, "joint_logL:", joint_logL
                #finally compute the log likelihood and time series of forward problem
                logL_ph = 0.0
                logL_rv = 0.0
                for i in range(len(S.Observables)):
                    if S.Observables[i].data_type == 'ph':
                          
                        logL_ph += S.Observables[i].compute_lc_series(write_lc=False)
                        #compute HR curve
                        S.Observables[i].obs_time = np.linspace(S.Observables[i].dt1-20.0, S.Observables[i].dt2+20.0, 
                                                    num=int(1440.0*(S.Observables[i].dt2-S.Observables[i].dt1)/S.Observables[i].t_exp))
                        
                        S.Observables[i].obs_data = []
                        #forward problem
                        S.Observables[i].compute_lc_series(write_lc=True)
                     
                    elif S.Observables[i].data_type == 'rv':
                        
                        logL_rv += S.Observables[i].compute_rv_series(write_rv=False)
                        #compute HR curve
                        S.Observables[i].obs_time = np.linspace(S.Observables[i].dt1-10.0, S.Observables[i].dt2+10.0, 
                                                    num=int(1440.0*(S.Observables[i].dt2-S.Observables[i].dt1)/S.Observables[i].t_exp))
                        S.Observables[i].obs_data = []
                        #forward problem
                        S.Observables[i].compute_rv_series(write_rv=True)
                        
                bunch_data.append([S.id, logL_ph, logL_rv, joint_logL])
                io.write_file('./maps' + 
                              '/spot_map_' + str(S.id) + '.dat', S.spot_map[0], ' ', '')
                
            return bunch_data
        
        #lock object
        lock = mp.Lock()
        #instantiate a Simulation object
        Simulation = StarsimSimulation()
        
        #read data and configuration
        try:
            #list of objects to joint minimize 
            OB = Simulation.read_multifit_config('./multifit.conf')
        except:
            print "Could not be read the multiobject configuration or no data selected"
            sys.exit()
  
        t_min, t_max = Simulation.dt1 - 2.*Simulation.spots_lifetime, Simulation.dt2 
        
        #generate an initial population of StarsimSimulation objects
        INIT_POPULATION = Glob.ncpus
        io.print_message("Inversion of Observational Data", 2, 91)
        io.print_message("Generating initial population: " + str(INIT_POPULATION) + " simulations", 2, 94)
        
        #inhibit the screen messages, for cleanness
        Glob.inhibit_msg = True
       
        Simulations = [StarsimSimulation() for i in range(INIT_POPULATION)]
        #prepare the objects
        j = 0
        for Simulation in Simulations:
            #inhibit the screen messages, for cleanness
            #Glob.inhibit_msg = True
            #load the data and specific parameters for this data, in each object of the population
            OB = Simulation.read_multifit_config('./multifit.conf')
            #Simulation.update_spectral_data(delta_t_sp=Delta_SP, delta_t_fc=Simulation.delta_t_fc)
            #check if photmetry object
            for i in range(len(OB)):
                if Simulation.Observables[i].data_type == 'ph':
                    #select a photometric zero point
                    Simulation.Observables[i].set_ph_zero_point(Simulation.Observables[i].z_ph)
                    
        pool = mp.Pool(Glob.ncpus)
        parallel_bunches = split_list(Simulations, Glob.ncpus)
     
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
        
        io.write_file("./stats.dat", SimulationsOut, ' ', '')
        
        """
        #compute the mean curves
        mean_ph_model = 0.0
        mean_rv_model = 0.0
        for simulation in SimulationsOut:
            for j in range(len(OB)):
                if simulation[0].Observables[j].data_type == 'ph':        
                    mean_ph_model += np.array(simulation[0].Observables[j].normalized_series) / len(SimulationsOut)
                    
                elif simulation[0].Observables[j].data_type == 'rv':
                    mean_rv_model += np.array(simulation[0].Observables[j].normalized_series) / len(SimulationsOut)
        
        #write data files
        for j in range(len(OB)):
            if simulation[0].Observables[j].data_type == 'rv':
                io.write_file("./output/MEAN_RV_MODEL.dat", mean_rv_model, ' ', '')            
            elif simulation[0].Observables[j].data_type == 'ph':
                io.write_file("./output/MEAN_PH_MODEL.dat", mean_ph_model, ' ', '')  
        
        #compute logL stastistic (MEAN MODEL - DATA)
        ph_residuals = []
        rv_residuals = []
        
        #MEAN CURVES
        for i in range(len(OB)):
                if Simulation.Observables[i].data_type == 'ph':
                    y_model, y_data_0 = mean_ph_model[:,1].astype(float), np.array(Simulation.Observables[i].y_data_0)
                    y_data = np.array(Simulation.Observables[i].y_data)
                    time_model = mean_ph_model[:,0].astype(float)
                    inv_sigma2 = 1.0/(Simulation.Observables[i].err**2 + Simulation.Observables[i].jitter**2)
                    inv_sigma2_0 = 1.0/(Simulation.Observables[i].err_0**2 + Simulation.Observables[i].jitter**2)
                    #interpolate highly sampled model
                    y_interpolated_model = model_interpolation(Simulation.Observables[i].obs_time_0, Simulation.Observables[i].obs_data_0, 
                                                               time_model, y_model)
                  
                    logL_0_ph = -0.5*(np.sum(((y_data - np.mean(y_data))**2)*inv_sigma2_0 + np.log(2.*np.pi) - np.log(inv_sigma2_0)) )
                    
                    logL_ph = -0.5*(np.sum(((y_data-y_interpolated_model)**2)*inv_sigma2 + np.log(2.*np.pi) - np.log(inv_sigma2)) )
     
                    for j in range(len(y_data)):
                        ph_residuals.append([Simulation.Observables[i].obs_time_0[j], y_data[j] - y_interpolated_model[j], Simulation.Observables[i].err[j] ])
                    
                    io.write_file("./output/MEAN_PH_MODEL.residuals", ph_residuals, ' ', '')
                    #write data series + offset + jitter
                    Simulation.Observables[i].write_lc_data()
                    
                    print "\n\tNull model logL:", Simulation.Observables[0].logL_null
                    print "\tlogL average photometric curve:", logL_ph
                    print "\tDlogL:", logL_ph - Simulation.Observables[0].logL_null
                    
                            
                elif Simulation.Observables[i].data_type == 'rv':
                    y_model, y_data = mean_rv_model[:,1].astype(float), Simulation.Observables[i].y_data_0
                    time_model = mean_rv_model[:,0].astype(float)
                    inv_sigma2 = 1.0/(Simulation.Observables[i].err**2 + Simulation.Observables[i].jitter**2)
                    
                    #interpolate highly sampled model
                    y_interpolated_model = model_interpolation(Simulation.Observables[i].obs_time_0, Simulation.Observables[i].obs_data_0, 
                                                               time_model, y_model)
          
                    logL_0_rv = -0.5*(np.sum(((y_data)**2)*inv_sigma2 + np.log(2.*np.pi) - np.log(inv_sigma2)) )
                    logL_rv = -0.5*(np.sum(((y_data-y_interpolated_model)**2)*inv_sigma2 + np.log(2.*np.pi) - np.log(inv_sigma2)) ) / logL_0_rv 

                    for j in range(len(y_data)):
                        rv_residuals.append([Simulation.Observables[i].obs_time_0[j], y_data[j] - y_interpolated_model[j], Simulation.Observables[i].err[j] ])
                    
                    io.write_file("./output/MEAN_RV_MODEL.residuals", rv_residuals, ' ', '')
                    #write data series + jitter
                    Simulation.Observables[i].write_rv_data()
                    
                    print "\n\tNull model logL:", Simulation.Observables[0].logL_null
                    print "\tlogL average RV curve:", logL_rv
                    print "\tDlogL:", logL_rv - Simulation.Observables[0].logL_null
        """
   
    elif num_spots != 0:
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
            dt1 = float(time_array_lc[0])
            dt2 = float(time_array_lc[-1])     
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
        print "Enter a valid mode --mode=[ph/rv/rv_hr/bootstrapping/forward...]"
        print '' 
                
        

    


























