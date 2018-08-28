import sys, os
from math import sin, cos, acos, sqrt, pi, log10, log, exp
from scipy.optimize import curve_fit, leastsq, minimize
import numpy as np
import bz2
import ConfigParser
import random
import copy
import multiprocessing as mp
from EnvGlobals import Glob
sys.path.append('./bin')
import mstarsim
import warnings


class Starsim_IO: 
    """
    """
    def __init__(self):
        #init object
        pass
    
    def print_message(self, str_, index, color):
        """
        """
        if not Glob.inhibit_msg:
            if not Glob.bw: 
                str_out = "\33["+str(index)+";"+str(color)+"m"+ str_ +"\33[1;m"
                print str_out
            else:
                print str_
                
    def read_bz2_file(self, in_file, separator):
        """read a bzip2 compressed file
        """
        with bz2.BZ2File(in_file, 'rb') as infile:
            if separator == 'single_line':
                g_list = []
                row = []
                while(True):
                    line = infile.readline()
                    if not line: break
                    if line[0] != '#'and line[0] != '\\' and line[0] != '|':
                        g_list.append(float(line))
            else:
                g_list = []
                row = []
                while(True):
                    line = infile.readline()
                    if not line: break
                    if line[0] != '#' and line[0] != '\\' and line[0] != '|':
                        if separator == ' ':
                            string = line.split()
                        else:
                            string = line.split(separator)
                        for col in range(0,len(string)):
                            row.append(float(string[col]))
                        g_list.append(row)
                        row = []
        #print "File read"
        return g_list

    def read_file(self, in_file, separator):
        with open(in_file, 'r') as infile:
            if separator == 'single_line':
                g_list = []
                row = []
                while(True):
                    line = infile.readline()
                    if not line: break
                    if line[0] != '#'and line[0] != '\\' and line[0] != '|':
                        g_list.append(float(line))
            else:
                g_list = []
                row = []
                while(True):
                    line = infile.readline()
                    if not line: break
                    if line[0] != '#' and line[0] != '\\' and line[0] != '|':
                        if separator == ' ':
                            string = line.split()
                        else:
                            string = line.split(separator)
                        for col in range(0,len(string)):
                            row.append(string[col])
                        if row != []:    
                            g_list.append(row)
                        row = []
        
        return g_list

    def write_file(self, out_file, data, separator, header):
        line = ""
        with open(out_file, 'w') as outfile:
            #outfile.write(str('#') + str(header) + "\n")
            for row in range(len(data)):
                for col in range(len(data[row])):
                    if col != len(data[row]) - 1:
                        line = line + str(data[row][col]) + str(separator)
                    elif col == len(data[row]) - 1:
                        line = line + str(data[row][col])
                outfile.write( line + "\n" )
                line = ""

    def write_file_bz2(self, out_file, data, separator, header):
        line = ""
        with bz2.BZ2File(out_file, 'wb') as outfile:
            #outfile.write(str('#') + str(header) + "\n")
            for row in range(len(data)):
                for col in range(len(data[row])):
                    if col != len(data[row]) - 1:
                        line = line + str(data[row][col]) + str(separator)
                    elif col == len(data[row]) - 1:
                        line = line + str(data[row][col])
                outfile.write( line + "\n" )
                line = ""

    def write_file_onecol(self, out_file, data, separator, header):
        line = ""
        with open(out_file, 'w') as outfile:
            #outfile.write(str('#') + str(header) + "\n")
            for row in range(len(data)):
                line = line + str(data[row]) + str(separator)
                outfile.write( line + "\n" )
                line = ""

    def get_science_data(self, file):
        """returns the time vector array of science lightcurve
        """
        lc = self.read_file(file, ' ')
        time_vector = list()
        data = list()
        sigma = list()
        
        if len(lc[0]) == 2:
            for i in range(len(lc)):
                #time_vector.append(float(lc[i][0])-297.073)
                time_vector.append(float(lc[i][0]))
                data.append([time_vector[i], float(lc[i][1]), 1.0]) #if no error column, fill with ones.
            
            return time_vector, data
        
        elif len(lc[0]) == 3:
            for i in range(len(lc)):
                #time_vector.append(float(lc[i][0])-297.073)
                time_vector.append(float(lc[i][0]))
                data.append([time_vector[i], float(lc[i][1]), float(lc[i][2])])
                
            return time_vector, data

class StarsimCommon(Starsim_IO):
    """
    """
    def __init__(self):
        #init object
        pass
    
    def flatten(self, lis):
        """return a nested list, flattened
        """
        f_lis = []
        for item in lis:
            if type(item) == type([]):
                f_lis.extend(flatten(item))
            else:
                f_lis.append(item)
        
        return f_lis
        
    def to_deg(self, radiants):
        """
        """
        return 180.0*radiants/pi
    
    def random_sin(self, z):
        """returns random number sin distributed [0:z] deg 
        """
        x = random.uniform(0,z)
        y = random.uniform(0,1)
        
        while y >= sin(pi*x/180.0):
             x = random.uniform(0,z)
             y = random.uniform(0,1)
        
        return x
            
    def impoly(self, x, in_coeffs):
        """return the image of a n-deg polynomial.
        input: x and array of n+1 coefficients, in decreasing order.
        """ 
        im = 0.0
        for order in range(0,len(in_coeffs)):
            im += in_coeffs[order]*(x**(len(in_coeffs)-order-1))
    
        return im
    
    def gauss(self, x, *p):
        """gauss function def
        """
        A, c, sigma, offset = p
        return A*np.exp(-((x-c)**2)/(2.0*sigma**2)) + offset

    def load_observational_data(self, file_path=''):
        """returns: [[time_array_lc],[lc_data]]
        """
        try:
            sys.argv[1]   #test if an argument is passed
            
        except NameError:
            #if no lc file path, use the time limits defined in starsim.conf
            self.obs_data = list()
            self.y_data = list()
            self.err = list()
            self.obs_time = np.linspace(self.dt1, self.dt2, num=int(1440.0*(self.dt2-self.dt1)/self.t_exp))
            
        else:
            try:
                #read a valid file path
                science_data_path = sys.argv[1]
                if not science_data_path[0] == '-':
                    #photometry
                    science_data = self.get_science_data(science_data_path)  
                    self.obs_time, self.obs_data = science_data[0], science_data[1]
                    self.y_data, self.err = np.array(science_data[1])[:,1], np.array(science_data[1])[:,2]
       
                else:
                    science_data = self.get_science_data(file_path)  
                    self.obs_time, self.obs_data = science_data[0], science_data[1]
                    self.y_data, self.err = np.array(science_data[1])[:,1], np.array(science_data[1])[:,2]
                    
                  
                #time limits
                self.dt1 = float(self.obs_time[0])
                self.dt2 = float(self.obs_time[-1])   
            except:
                #if no lc file path, use the time limits defined in starsim.conf
                self.obs_data = list()
                self.y_data = list()
                self.err = list()
                self.obs_time = np.linspace(self.dt1, self.dt2, num=int(1440.0*(self.dt2-self.dt1)/self.t_exp))
        #generate secundary vectors to store observational data (in case primary vectors were changed)
        self.obs_time_0, self.obs_data_0, self.y_data_0, self.err_0 = np.array(self.obs_time), np.array(self.obs_data), np.array(self.y_data), np.array(self.err)     
        
        return self.obs_time, self.obs_data, self.y_data, self.err
    
class Star(Starsim_IO):
    """star class: packs all parameters related to star parameters
    """   
    def __init__(self, conf_file):
        #init configuration file
        self.conf_file = conf_file
        #------------------
        #stellar parameters 
        #------------------
        #photosphere temperature
        self.temp_ph = float(self.conf_file.get('star','t_eff_ph'))
        #modulus 10 to round temperatures to use CCF precomputed
        self.temp_ph = self.temp_ph - (self.temp_ph % 10.0)
        #delta T spots
        self.delta_t_sp = float(self.conf_file.get('star','spot_T_contrast'))
        self.delta_t_sp = self.delta_t_sp - (self.delta_t_sp % 10.0)
        #delta T faculae
        self.delta_t_fc = float(self.conf_file.get('star','faculae_T_contrast'))
        self.delta_t_fc = self.delta_t_fc - (self.delta_t_fc % 10.0)
        #faculae / spot surface ratio
        self.Q = float(self.conf_file.get('star','q_ratio'))
        #logg stellar surface
        self.grav = float(self.conf_file.get('star','logg'))
        #rotational period
        self.P_rot = float(self.conf_file.get('star','p_rot'))
        #star axis inclination
        self.x_i = float(self.conf_file.get('star','axis_i'))
        #spot evolution rate (deg/day)
        self.evo_rate = float(self.conf_file.get('star','spot_size_evo_rate'))
        #differential rotation parameters
        self.diff_rotation = float(self.conf_file.get('star','diff_rotation')) 
        self.B_coeff = float(self.conf_file.get('star','B_rot'))
        self.C_coeff = float(self.conf_file.get('star','C_rot'))
        #DERIVED PARAMS.
        self.w = 2.0*pi/self.P_rot
        self.m_star = sqrt((self.temp_ph/5777.0)**4.0)
        #self.R_star = sqrt((self.temp_ph/5777.0)**3.0015*27440.03/10.0**self.grav)   #star radius in Rsun
        self.R_star = float(self.conf_file.get('star','R_star'))
        if self.R_star == -1:
            self.R_star = sqrt( (self.m_star**3.5) / (self.temp_ph/5777.0)**4.0  )
    
    
    def __del__(self):
        """class destructor
        """
    
class Spectra(StarsimCommon, Starsim_IO, Star):
    """
    """
    def __init__(self):
        #init 
        self.conf_file = __conf_init(self)
    
    def __conf_init(self, conf_file_path='./starsim.conf'):
        """creates an instance of class ConfigParser, and read the .conf file. Returns the object created
        ConfigParser, and read the .conf file. Returns the object created
        """
        conf_file_Object = ConfigParser.ConfigParser()
        if not conf_file_Object.read([  conf_file_path]):
            print_message( str(conf_file_path) + ": could not be read", 5, 31, inhibit_msg=False)
            sys.exit()
        else:
            return conf_file_Object
    
    def compute_cth_coeffs(self):
        """read the LD factors from the generated Kurucz interpolation
        """
        self.cth_coeffs = self.read_file("./data/cth.dat",'single_line')
   
        return self.cth_coeffs
    
    def read_precompiled_spectra(self):
        """read the spectra from file
        """   
        #read path of precompiled low-res spectra directory
        path = str(self.conf_file.get('files','precomp_spectra'))
        self.Sim.temp_ph = self.Sim.temp_ph - (self.Sim.temp_ph % 10.0)
        #temperatures of photosphere, spots and faculae
        self.temp_sp = (self.Sim.temp_ph - self.Sim.delta_t_sp) - (self.Sim.temp_ph - self.Sim.delta_t_sp) % 10.0
        self.temp_fc = (self.Sim.temp_ph + self.Sim.delta_t_fc) - (self.Sim.temp_ph + self.Sim.delta_t_fc) % 10.0
                
        if self.spectra_res == 'high':
            nspres = 63076
            #number of points depending on spectral range
            if self.wvmax >= 1350.0 and self.wvmin < 1350.0:
                n_wv = 20*(1350.0 - self.wvmin) + 10*(self.wvmax - 1350.0)
            elif self.wvmin > 2500.0:
                n_wv = 5*(wvmax - wvmin)
            elif self.wvmax < 1350.0:
                n_wv = 20*(self.wvmax - self.wvmin)
            elif self.wvmax > 1350.0 and self.wvmin >= 1350.0 and self.wvmax < 2500.0:
                n_wv = 10*(self.wvmax - self.wvmin)
            elif wvmax > 2500.0 and wvmin <= 2500.0:
                n_wv = 10*(2500.0 - self.wvmin) + 5*(self.wvmax - 2500.0)

            file_name_ph = path + "/precompiled_spectra/spec_" + str(int(self.Sim.temp_ph)) + '_' + str(self.Sim.grav) +  \
                        '_' + str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '.dat.bz2'
            file_name_sp = path + "/precompiled_spectra/spec_" + str(int(self.temp_sp)) + '_' + str(self.Sim.grav) +  \
                        '_' + str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '.dat.bz2'
            file_name_fc = path + "/precompiled_spectra/spec_" + str(int(self.temp_fc)) + '_' + str(self.Sim.grav) +  \
                        '_' + str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '.dat.bz2'
                
        elif self.spectra_res == 'low':
            
            file_name_ph = path + "/precompiled_spectra/spec_" + str(int(self.Sim.temp_ph)) + '_' + str(self.Sim.grav) + '_' + str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '_lres.dat.bz2'
            file_name_sp = path + "/precompiled_spectra/spec_" + str(int(self.temp_sp)) + '_' + str(self.Sim.grav) + '_' + str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '_lres.dat.bz2'
            file_name_fc = path + "/precompiled_spectra/spec_" + str(int(self.temp_fc)) + '_' + str(self.Sim.grav) + '_' + str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '_lres.dat.bz2'
        
        #try to read precompiled spectra, if exists. If not, generate and read again
        try:
            spectra_ph = self.read_bz2_file(file_name_ph, ' ')
        except IOError:
            self.gen_spectra(self.Sim.temp_ph, self.Sim.grav, 'ph', self.conf_file)
            spectra_ph = self.read_bz2_file(file_name_ph, ' ') 
        try:
            spectra_sp = self.read_bz2_file(file_name_sp, ' ')
        except IOError:
            self.gen_spectra(self.temp_sp, self.Sim.grav, 'sp', self.conf_file)
            spectra_sp = self.read_bz2_file(file_name_sp, ' ')
        try:
            spectra_fc = self.read_bz2_file(file_name_fc, ' ')
        except IOError:
            self.gen_spectra(self.temp_fc, self.Sim.grav, 'fc', self.conf_file)
            spectra_fc = self.read_bz2_file(file_name_fc, ' ')
        
        self.print_message("", 0, 94)
        self.print_message("Reading .bz2 low-res. precompiled spectra...", 0, 94)
        self.print_message("    Read photosphere spectra:" + str(file_name_ph), 3, 32)
        self.print_message("    Read spot spectra:" + str(file_name_sp), 3, 32)
        self.print_message("    Read faculae spectra:" + str(file_name_fc), 3, 32)
        
        #global wv, flnp_ph, flnp_sp, flnp_fc, flpk_ph, flpk_sp, flpk_fc, dlnp_ph, dlnp_sp, dlnp_fc
        
        if self.spectra_res == 'high':
            self.wv = np.zeros(int(n_wv), dtype=float)
            
        elif self.spectra_res == 'low':
            n_wv_ph, n_wv_sp, n_wv_fc  = len(spectra_ph), len(spectra_sp), len(spectra_fc)
            #find the min number of points in all spectras
            n_wv = min(n_wv_ph,n_wv_sp,n_wv_fc)
            self.print_message("    (Number of points in the spectras: " + str(n_wv) + ')', index=0, color=94)
            self.wv = np.zeros(int(n_wv), dtype=float)
        
        self.flnp_ph, self.flnp_sp, self.flnp_fc = np.zeros(int(n_wv), dtype=float), np.zeros(int(n_wv), dtype=float), np.zeros(int(n_wv), dtype=float)
        self.flpk_ph, self.flpk_sp, self.flpk_fc = np.zeros(int(n_wv), dtype=float), np.zeros(int(n_wv), dtype=float), np.zeros(int(n_wv), dtype=float)
        self.dlnp_ph, self.dlnp_sp, self.dlnp_fc = np.zeros((n_wv,18), dtype=float), np.zeros((n_wv,18), dtype=float), np.zeros((n_wv,18), dtype=float)
        
        for i in range(int(n_wv)):
            self.wv[i] = spectra_ph[i][0]
            self.flnp_ph[i], self.flnp_sp[i], self.flnp_fc[i] = spectra_ph[i][1], spectra_sp[i][1], spectra_fc[i][1]
            self.flpk_ph[i], self.flpk_sp[i], self.flpk_fc[i] = spectra_ph[i][2], spectra_sp[i][2], spectra_fc[i][2]
            for j in range(18):
                self.dlnp_ph[i][j], self.dlnp_sp[i][j], self.dlnp_fc[i][j] = spectra_ph[i][j+3], spectra_sp[i][j+3], spectra_fc[i][j+3]

    def btsettl_interpolation(self, BTSettl_params, conf_file):
        """returns an array of interpolated_BTSettl spectra given a temperature and logg
        """
        temp, grav, nm = BTSettl_params[0],BTSettl_params[1],BTSettl_params[2]
        #config params
        write_tmp_files = True
        #contant block
        ntemp = 44         #we've 44 different temperatures in files
        ngrav = 12
        maxwaves = 630760  #number of points in BT-Settl files
        out_spectra = np.zeros([maxwaves/10,2])
        #array declaration
        tmod, gmod, dgr = np.zeros(int(ntemp)), np.zeros(int(ngrav)), np.zeros(int(ngrav))
        wav, fint = np.zeros(int(maxwaves)), np.zeros(int(maxwaves))
        #build the grids
        for j in range(ntemp):
            tmod[j] = 2800.0 + 100.0*(j-1)
        for j in range(ngrav):
            gmod[j] = 0.5*(j-1)
            dgr[j] = abs(gmod[j] - grav)
        gmof = gmod[0]
        #find minimum of grav. Which file is closer to a given logg
        for j in range(1,ngrav):
            if dgr[j] < dgr[j-1]: gmof = gmod[j]
        #finds the 2 files limiting the given temperature
        for j in range(ntemp):
            if temp >= tmod[j] and temp < tmod[j+1]:
                ite1_f = tmod[j]
                ite2_f = tmod[j+1]
        path = str(conf_file.get('files','bt_settl_path'))
        #define the file names
        if self.spectra_res == 'high':
            name_ite1 =  path + '/lte0' + str(int(ite1_f)/100) + '-' + str(gmof) + '-0.0a+0.0.BT-Settl.spec.7_mod'
            name_ite2 = path + '/lte0' + str(int(ite2_f)/100) + '-' + str(gmof) + '-0.0a+0.0.BT-Settl.spec.7_mod'
        elif self.spectra_res == 'low':
            name_ite1 =  path + '/lte0' + str(int(ite1_f)/100) + '-' + str(gmof) + '-0.0a+0.0.BT-Settl.spec.7_mod_lres'
            name_ite2 = path + '/lte0' + str(int(ite2_f)/100) + '-' + str(gmof) + '-0.0a+0.0.BT-Settl.spec.7_mod_lres'
        #and read it
        ite1 = self.read_file(name_ite1, ' ')
        ite2 = self.read_file(name_ite2, ' ')
    
        if self.spectra_res == 'low':
            wav, fint = np.zeros(int(min(len(ite1), len(ite2)))), np.zeros(int(min(len(ite1), len(ite2))))
            out_spectra = np.zeros([min(len(ite1)/10,len(ite2)/10),2])
        
        elif self.spectra_res == 'high':
            wav, fint = np.zeros(int(min(len(ite1), len(ite2)))), np.zeros(int(min(len(ite1), len(ite2))))
            out_spectra = np.zeros([len(ite2)/10,2])
        
        #interpolation
        for j in range(min(len(ite1), len(ite2))):
            wav[j] = float(ite1[j][0])/10.0   #wavelength in nanometers
            fsint1, fsint2 = 10.0**float(ite1[j][1]), 10.0**float(ite2[j][1])
            fint[j] = fsint1 + ( (temp - ite1_f)/(ite2_f - ite1_f) )*(fsint2 - fsint1)
        
        for i in range(min(len(ite1), len(ite2))/10):
            wvsum = 0.0
            flsum = 0.0
            for j in range(10):
                wvsum = wav[10*(i)+j] + wvsum
                flsum = fint[10*(i)+j] + flsum
            out_spectra[i][0] = wvsum/10.0
            out_spectra[i][1] = flsum/10.0
            
        #if want to write .tmp files with spectra
        if write_tmp_files:
            if nm == 'ph':
                self.write_file('./data/specph.tmp', out_spectra, ' ','')
            elif nm == 'sp':
                self.write_file('./data/specsp.tmp', out_spectra, ' ','')
            elif nm == 'fc':
                self.write_file('./data/specfc.tmp', out_spectra, ' ','')

        return out_spectra
   
    def phoenix_interpolation(self, Phoenix_params, conf_file):
        """interpolate Phoenix-HR spectra for Temp and logg and type (ph,sp,fc)
        """
        import warnings
        warnings.filterwarnings("ignore")
        
        write_tmp_files = True
        temp, grav, nm = Phoenix_params[0], Phoenix_params[1], Phoenix_params[2]
        maxwaves = 1569128  #number of points in PhoenixHR
        ntemp = 73          #number of different temperatures
        ngrav = 3
        out_spectra = np.zeros( [maxwaves,2] )
        #array declaration
        tmod, gmod, dgr = np.zeros( int(ntemp)), np.zeros(int(ngrav)), np.zeros(int(ngrav) )
        wav, fint = np.zeros( int(maxwaves)), np.zeros(int(maxwaves) )
        #build the grids
        for j in range(0,ntemp):
            if j < 48:
                tmod[j] = 2300.0 + 100.0*(j)
            else:
                tmod[j] = 7000.0 + 200.0*(j-49)
                
        for j in range(0,ngrav):
            gmod[j] = 4.0 + 0.5*(j)
            dgr[j] = abs(gmod[j] - grav)   
        #find minimum of grav. Which file is closer to a given logg
        gmof = gmod[0]
        for j in range(1,ngrav):
            if dgr[j] < dgr[j-1]: gmof = gmod[j]
        #finds the 2 files limiting the given temperature
        for j in range(0,ntemp):
            if temp >= tmod[j] and temp < tmod[j+1]:
                ite1_f = tmod[j]
                ite2_f = tmod[j+1]
        path = str(conf_file.get('files','phoenix_path'))
        #define the file names
        name_ite1 =  path + '/lte0' + str(int(ite1_f)) + '-' + str(gmof) + str('0') + '-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes2.dat.bz2'
        name_ite2 = path + '/lte0' + str(int(ite2_f)) + '-' + str(gmof) + str('0') +'-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes2.dat.bz2'
        #and read it 
        ite1 = self.read_bz2_file(name_ite1, ' ')
        ite2 = self.read_bz2_file(name_ite2, ' ')
       
        for j in range(len(min(ite1,ite2))):
            wav[j] = float(ite1[j][0])/10.0
            #fsint1, fsint2 = 10.0**float(ite1[j][1]), 10.0**float(ite2[j][1])
            fsint1, fsint2 = float(ite1[j][1]), float(ite2[j][1])
            fint[j] = fsint1 + ( (temp - ite1_f)/(ite2_f - ite1_f) )*(fsint2 - fsint1)
        
        #vacuum-to-air correction
        for i in range(len(wav)):
            if self.spectral_mask == 'IR':
                out_spectra[i][0] = wav[i]
                out_spectra[i][1] = fint[i]
            else:
                out_spectra[i][0] = wav[i]/(1.0 + 5.792105e-2/(238.01850-((1.e3/wav[i])**2.0)) + 1.67917e-3/(57.3620-((1.e3/wav[i])**2.0)))
                out_spectra[i][1] = fint[i]
        """
        normalize and range spectra in the following way:
                1- Cut the spectra with wavelength limits
                2- Divide the spectral range in bins
                3- Search the maximum of each bin
                4- Fit a 6-deg polinomial to this points
                5- Normalize the spectra with respect to this polynomial
        """
        out_spectra_cut = list()
        for k in range(len(out_spectra)):
            if out_spectra[k][0] > self.wvmin and out_spectra[k][0] < self.wvmax:
                out_spectra_cut.append([out_spectra[k][0], out_spectra[k][1]])
        
        #write_file('./data/spec_cut.dat', out_spectra_cut, ' ', '')
        n_points_spec = len(out_spectra_cut)
        #select the number of bins
        if n_points_spec > 100000:
            nbin = 1000
        elif n_points_spec < 10000:
            nbin = 50
        else:
            nbin = 100
        
        n_points_bin = (n_points_spec / nbin) - 1 
        x_bin = []
        y_bin = []
        #width in lambda of each bin
        lambda_bin_width = (self.wvmax - self.wvmin) / nbin
        
        for j in range(nbin):
            wv_max_bin = self.wvmin + (j+1)*lambda_bin_width
            wv_min_bin = wv_max_bin - lambda_bin_width
            fmax, wlenmax = 0.0, 0.0
            for i in range(len(out_spectra_cut)):
                if out_spectra_cut[i][0] <= wv_max_bin and out_spectra_cut[i][0] >= wv_min_bin:
                    if out_spectra_cut[i][1] > fmax:
                        wlenmax, fmax = out_spectra_cut[i][:]
            x_bin.append(wlenmax)
            y_bin.append(fmax)    
        #fit a polynomial and subtract the black body continuum         
        bb_fit_coeffs = np.polyfit(np.array(x_bin), np.array(y_bin), deg=6, rcond=None, full=False, w=None, cov=False) 
        out_spectra_normalized = []
        #normalize with respect to the fit
        for j in range(len(out_spectra_cut)):
            out_spectra_normalized.append([out_spectra_cut[j][0],out_spectra_cut[j][1]/self.impoly(out_spectra_cut[j][0], bb_fit_coeffs)])
        
        #if want to write .tmp files with spectra
        if write_tmp_files:
            ref_file = str(int(random.uniform(11,99)))
            if nm == 'ph':
                self.write_file('./data/tmp/phoenixHR_ph_norm.tmp' + ref_file, out_spectra_normalized, ' ','')
            elif nm == 'sp':
                self.write_file('./data/tmp/phoenixHR_sp_norm.tmp' + ref_file, out_spectra_normalized, ' ','')
            elif nm == 'fc':
                self.write_file('./data/tmp/phoenixHR_fc_norm.tmp' + ref_file, out_spectra_normalized, ' ','')
                
        return out_spectra_cut, ref_file

    def gen_spectra(self, temperature, grav, surf_type, conf_file):
        """generates the spectra files
        """ 
        self.print_message("Spectra not found in ./precompiled_spectra. Generating...", 0, 91)
        self.print_message("Reading BT-Settl low resolution spectra...", 3, 32)
        
        if self.spectra_res == 'high':
            nspres = 63076
            if self.wvmax >= 1350.0 and self.wvmin < 1350.0:
                n_wv = 20*(1350.0 - self.wvmin) + 10*(self.wvmax-1350.0)
            elif self.wvmin >= 2500.0:
                n_wv = 5*(self.wvmax - self.wvmin)
            elif self.wvmax < 1350.0:   
                n_wv = 20*(self.wvmax - self.wvmin)
            elif self.wvmax > 1350.0 and self.wvmin >= 1350.0 and self.wvmax < 2500.0:
                n_wv = 10*(self.wvmax - self.wvmin)
            elif self.wvmax >= 2500.0 and self.wvmin < 2500.0:
                n_wv = 10*(2500.0 - self.wvmin) + 5*(self.wvmax - 2500.0)
                
        #Generation the SEDs for photosphere, spots and faculae with BTSettl
        self.print_message("Generating synthetic spectra...", 0, 91)
        #print "With", n_wv ,'points'

        self.print_message("Low Resolution BT-settl data base", 0, 91)
       
        nm = 'ph'
        interpolated_BTSettl_ph = self.btsettl_interpolation([temperature, grav, surf_type], conf_file)

        #Generation of Kurucz SED with the LD coefficients
        self.print_message("Generation of Kurucz spectra with LD coefficients", 0, 91)
        nm = 'ph'
        if temperature < 3500.0:
            temp_photk = 3500.0
            self.print_message("LD for photosphere temperature out of ATLAS9 range", 0, 91)
        else:
            temp_photk = temperature
        
        kurucz_phot = mstarsim.kurucz_interpol(temp_photk, grav, surf_type)
        #reads the LD factors from the generated Kurucz interpolation
        cth_coeffs = self.read_file("./data/cth.dat",'single_line')
        acd = np.arccos(cth_coeffs)

        #generates the BTSettl interpolated SEDs with LD factors from Kurucz_interpol

        #for photosphere
        if self.spectra_res == 'high':
            wv, dlnp_ph, flpk_ph, flnp_ph = mstarsim.limb_darkening_to_btsettl(surf_type, nspres, self.wvmin, n_wv, acd)
            
        elif self.spectra_res == 'low':
            wv_, dlnp_ph_, flpk_ph_, flnp_ph_ = mstarsim.limb_darkening_to_btsettl_fast(surf_type, len(interpolated_BTSettl_ph), self.wvmin, self.wvmax, acd)
            n_wv = 0
            wv = []
            dlnp_ph = []
            flpk_ph = []
            flnp_ph = []
            
            for i in range(len(wv_)):
                if wv_[i] != 0.0: 
                    n_wv += 1
                    wv.append(wv_[i])
                    dlnp_ph.append(dlnp_ph_[i])
                    flpk_ph.append(flpk_ph_[i])
                    flnp_ph.append(flnp_ph_[i])
    
        self.print_message("Number of points in ranged spectra: " + str(n_wv-1), 0, 91)
        
        data_block = list()
       
        for i in range(n_wv-1):
            row_dlnp = ''
            for j in range(18):
                row_dlnp = row_dlnp + str(' ') + str(dlnp_ph[i][j])
            data_block.append([wv[i], flnp_ph[i], flpk_ph[i], row_dlnp])

        path = str(conf_file.get('files','precomp_spectra'))
        
        if self.spectra_res == 'low':
            file_name = path + '/precompiled_spectra/spec_' + str(int(temperature)) + '_' + str(self.Sim.grav) + '_' + str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '_lres.dat.bz2'
    
        elif self.spectra_res == 'high':
            file_name = path + '/precompiled_spectra/spec_' + str(int(temperature)) + '_' +  \
                        str(self.Sim.grav) + '_' + str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '.dat.bz2'
        #write this new file to disk
        self.write_file_bz2(file_name, data_block, ' ','')
    
    def compute_immaculate_photosphere(self):
        """read precalculated immaculate (quiet) photosphere if in file, otherwise compute it. 
        """
        
        #name of the file
        if self.spectra_res == 'high':
            flt_immaculate_ph_file = './data/flt_immaculate_ph_' + str(int(self.Sim.temp_ph)) + '_' +  \
                                    str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '.dat'
        elif self.spectra_res == 'low':
            flt_immaculate_ph_file = './data/flt_immaculate_ph_' + str(int(self.Sim.temp_ph)) + '_' +  \
                                    str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '_lres.dat'
        try:
            #check if immaculate photosphere file exists
            self.flt_immaculate_ph = self.read_file(flt_immaculate_ph_file, ' ')
            #self.print_message("    Read immaculate photosphere file", index=3, color=32)
            
            return self.flt_immaculate_ph
        
        except IOError:
            #if not in file, compute it
            self.print_message("    Computing immaculate photosphere", 0, 91)
            
            star_params = [self.Sim.x_i, self.Sim.R_star, self.Sim.w, self.Sim.B_coeff, self.Sim.C_coeff, 
                           self.Sim.temp_ph, self.Sim.delta_t_sp, self.Sim.delta_t_fc]
            self.flt_immaculate_ph, self.flxph, self.S_index_im = mstarsim.compute_immaculate_photosphere(star_params, self.cth_coeffs, 
                                            self.flpk_ph, self.dlnp_ph, self.flnp_ph, n_wv=len(self.wv))
            #write immaculate photosphere spectra to a file (flt_immaculate_ph_file)
            self.write_file_onecol(flt_immaculate_ph_file , self.flt_immaculate_ph, ' ', '')
        
            return self.flt_immaculate_ph
    
    def immaculate_photosphere_flux(self):
        """
        """
        star_params = [self.Sim.x_i, self.Sim.R_star, self.Sim.w, self.Sim.B_coeff, self.Sim.C_coeff, 
                       self.Sim.temp_ph, self.Sim.delta_t_sp, self.Sim.delta_t_fc]
        
        self.flt_immaculate_ph, self.flxph, self.S_index_im = mstarsim.compute_immaculate_photosphere(star_params, self.cth_coeffs, 
                                            self.flpk_ph, self.dlnp_ph, self.flnp_ph, n_wv=len(self.wv))
        #write immaculate photosphere spectra to a file (flt_immaculate_ph_file)
        #self.write_file_onecol(flt_immaculate_ph_file , self.flt_immaculate_ph, ' ', '')
        
        return self.flt_immaculate_ph, self.S_index_im

    def init_ccfs(self):
        """read or compute the ccf's for the given temperatures.
        """
        self.print_message('',index=7, color=33)
        #self.print_message('Radial Velocity series calculation...', index=7, color=31)
        #interpolate Phoenix HR spectra with the given temperatures
        self.print_message('Computing / Reading CCFs...', index=6, color=35)
        
        #mask used for CCF generation
        self.spectral_mask = str(self.conf_file.get('rv','spectral_mask'))
        if self.spectral_mask == 'IR':
            base_path = './precompiled_ccf_VIS_IR/ccf_'
        else:
            base_path = './precompiled_ccf/ccf_'
        #compute interpolated spectras for different surface types
        for surf_type in ('ph','sp','fc'):
            if surf_type == 'ph':
                temperature, name = self.Sim.temp_ph, 'photosphere'
                file_ph = base_path + str(self.Sim.temp_ph) + '_' + str(self.wvmin) + '_' + str(self.wvmax) + '_' + str(self.ccf_vel_range) + '.dat'
                try:
                    ccf_data = self.read_file(file_ph, ' ')
                    ccf_ph = [float(ccf_data[i][1]) for i in range(len(ccf_data))]
                    rv = [float(ccf_data[i][0]) for i in range(len(ccf_data))]
                except IOError:
                    self.print_message('    Interpolating Phoenix HR spectra for ' + str(name), index=6, color=34)
                    phoenix_out_spectra, ref_file = self.phoenix_interpolation([temperature, self.Sim.grav, surf_type], self.conf_file)
                    self.print_message('    Computing CCF for ' + str(name), index=6, color=33) 
                    rv, ccf_ph = mstarsim.hcorr(str(surf_type), self.spectral_mask, ccfrng=self.ccf_vel_range, stp=0.25, nrv=2*int(self.ccf_vel_range/0.25), ref_file=str(ref_file))
                    #write the CCF on disk
                    self.write_file(file_ph, zip(rv,ccf_ph), ' ', '')
                #fit ccf bisectors
                bisector, bis_coeffs = self.fit_ccf_bisector(rv, ccf_ph)    
                bisector_ph, bis_coeffs_ph = bisector, bis_coeffs
            
            elif surf_type == 'sp':
                temperature, name = self.Sim.temp_ph - self.Sim.delta_t_sp, 'spots'
                file_sp = base_path + str(self.Sim.temp_ph - self.Sim.delta_t_sp) + '_' + str(self.wvmin) + '_' + str(self.wvmax) + '_' + str(self.ccf_vel_range) + '.dat'
                try:
                    ccf_data = self.read_file(file_sp, ' ')
                    ccf_sp = [float(ccf_data[i][1]) for i in range(len(ccf_data))]
                    rv = [float(ccf_data[i][0]) for i in range(len(ccf_data))]
                except IOError:
                    self.print_message('    Interpolating Phoenix HR spectra for ' + str(name), index=6, color=34)
                    phoenix_out_spectra, ref_file = self.phoenix_interpolation([temperature, self.Sim.grav, surf_type], self.conf_file)
                    self.print_message('    Computing CCF for ' + str(name), index=6, color=33) 
                    rv, ccf_sp = mstarsim.hcorr(str(surf_type), self.spectral_mask, ccfrng=self.ccf_vel_range, stp=0.25, nrv=2*int(self.ccf_vel_range/0.25), ref_file=str(ref_file))
                    #write the CCF on disk
                    self.write_file(file_sp, zip(rv,ccf_sp), ' ', '')
                #fit ccf bisectors
                bisector, bis_coeffs = self.fit_ccf_bisector(rv, ccf_sp)    
                bisector_sp, bis_coeffs_sp = bisector, bis_coeffs
                
            elif surf_type == 'fc':
                temperature, name = self.Sim.temp_ph + self.Sim.delta_t_fc, 'faculae' 
                file_fc = base_path + str(self.Sim.temp_ph + self.Sim.delta_t_fc) + '_' + str(self.wvmin) + '_' + str(self.wvmax) + '_' + str(self.ccf_vel_range) + '.dat'    
                try:
                    ccf_data = self.read_file(file_fc, ' ')
                    ccf_fc = [float(ccf_data[i][1]) for i in range(len(ccf_data))]
                    rv = [float(ccf_data[i][0]) for i in range(len(ccf_data))]
                except IOError:
                    self.print_message('    Interpolating Phoenix HR spectra for ' + str(name), index=6, color=34)
                    phoenix_out_spectra, ref_file = self.phoenix_interpolation([temperature, self.Sim.grav, surf_type], self.conf_file)
                    self.print_message('    Computing CCF for ' + str(name), index=6, color=33) 
                    rv, ccf_fc = mstarsim.hcorr(str(surf_type), self.spectral_mask, ccfrng=self.ccf_vel_range, stp=0.25, nrv=2*int(self.ccf_vel_range/0.25), ref_file=str(ref_file))
                    #write the CCF on disk
                    self.write_file(file_fc, zip(rv,ccf_fc), ' ', '')
                #fit ccf bisectors
                bisector, bis_coeffs = self.fit_ccf_bisector(rv, ccf_fc)    
                bisector_fc, bis_coeffs_fc = bisector, bis_coeffs     
            
        #compute the symmetric CCF
        self.print_message('    Computing symmetric CCF for photosphere...', index=6, color=33)
        self.rv_sym, self.ccf_sym_ph = self.gen_sym_ccf(rv, ccf_ph, bis_coeffs_ph)
        #write_file_onecol('./data/ccf_sym_ph.dat', ccf_sym_ph, ' ', '')
        self.print_message('    Computing symmetric CCF for spots...', index=6, color=33)
        self.rv_sym, self.ccf_sym_sp = self.gen_sym_ccf(rv, ccf_sp, bis_coeffs_sp)
        #write_file_onecol('./data/ccf_sym_sp.dat', ccf_sym_sp, ' ', '')
        self.print_message('    Computing symmetric CCF for faculae...', index=6, color=33)
        self.rv_sym, self.ccf_sym_fc = self.gen_sym_ccf(rv, ccf_fc, bis_coeffs_fc)
        
        #get the flux of immaculate photosphere
        self.print_message('Computing CCF of rotating immaculate photosphere', index=6, color=34)
        self.flt_immaculate_ph = self.compute_immaculate_photosphere()
        self.flxph = self.immaculate_photosphere_flux()
        
        #compute the immaculate photosphere ccf
        self.Sim.w = 2.0*pi/self.Sim.P_rot
        star_params = [self.Sim.x_i, self.Sim.R_star, self.Sim.w, self.Sim.B_coeff, self.Sim.C_coeff, 
                       self.Sim.temp_ph, self.Sim.delta_t_sp, self.Sim.delta_t_fc]
        self.rv_im, self.ccf_im = mstarsim.compute_ccf_immaculate_photosphere(star_params[0:3], self.cth_coeffs, self.flpk_ph, self.dlnp_ph, 
                                                                        self.flnp_ph, self.flxph, self.rv_sym, self.ccf_sym_ph, len(self.flpk_ph), len(self.ccf_sym_ph))
       
        ccf_immaculate_ph = list()
        for i in range(len(self.ccf_im)): ccf_immaculate_ph.append([self.rv_im[i], self.ccf_im[i]])
        
        #systemic velocity (RV offset) 
        self.sys_vel = self.systemic_vel(self.rv_sym, self.ccf_im)
    
        return self.flt_immaculate_ph, self.flxph, self.rv_sym, self.rv_im, self.ccf_im, self.ccf_sym_ph, self.ccf_sym_sp, self.ccf_sym_fc, self.sys_vel
       
class Spots(StarsimCommon, Starsim_IO, Star):        
    """
    """
   
    def __init__(self, conf_file):
        """init configuration file
        """
        self.conf_file = conf_file
        #spots 
        self.spots_lifetime = float(self.conf_file.get('spots','spots_lifetime'))
        self.sigma_spots_lifetime = float(self.conf_file.get('spots','spots_lifetime_sigma'))
        self.colat_limits = [0.0, 180.0]
        self.long_limits = [0.0, 360.0]
        self.r_limits = [0.0, 8.0]
        self.n_spots = 0
                
    def __conf_init(self, conf_file_path):
        """creates an instance of class ConfigParser, and read the .conf file. Returns the object created
        ConfigParser, and read the .conf file. Returns the object created
        """
        conf_file_Object = ConfigParser.ConfigParser()
        if not conf_file_Object.read([  conf_file_path]):
            print_message( str(conf_file_path) + ": could not be read", 5, 31, inhibit_msg=False)
            sys.exit()
        else:
            return conf_file_Object
    
    def read_spot_params(self, spot_list_file):
        """Load the number of spots and read parameters
        """        
        try:
            spot_map = self.read_file(spot_list_file, ' ')
            spot_map = np.array(spot_map).astype(float)
            #in case of end file
            if str(spot_map[0][0]) == 'END': return False
            self.n_spots = len(spot_map)
            n_col = len(spot_map[0])
            self.spot_params =  [ spot_map, self.evo_rate, self.diff_rotation, self.n_spots, self.Q ]
            #self.print_message( "Process id: " + str(self.id), index=10, color=31)
            self.print_message( "Reading SpotMap file... ", index=10, color=36)
            self.print_message( "    Number of spots on star: " +  str(self.n_spots), index=10, color=37)
      
        except IOError:
            self.print_message(str(8*' ') + "IO Error in reading spot file. Apparently there is not any file ./data/spot_list.cat", index=5, color=91)
            sys.exit()
        
        return self.spot_params
    
    def add_spot(self):
        """add a random spot in the current map
        """
        self.spot_map[0].extend( self.gen_rnd_map(num_spots=1) )
        self.n_spots += 1
    
    def map_perturbation(self):
        """
        """
        row = int(random.uniform(0,len(self.spot_map[0])))
        col = int(random.uniform(0,5))
        """
        new_spot = [random.uniform(self.dt1-self.spots_lifetime, self.dt2), 
                    abs(random.gauss(self.spots_lifetime, self.sigma_spots_lifetime)),
                    random.uniform(0., 360.), random.uniform(0.0, 180.0),
                    random.uniform(0.0, 8.0) ]
        """
        if col == 0:
            new_t_ini = random.uniform(self.dt1-self.spots_lifetime, self.dt2)    
            self.spot_map[0][row][0] = new_t_ini
            #spots lifetime     
        elif col == 1: 
            new_t_life = abs(random.gauss(self.spots_lifetime, self.sigma_spots_lifetime))
            self.spot_map[0][row][1] = new_t_life
            #if colatitude is chosen
        elif col == 2:
            colat = random.uniform(0., 180.)
            self.spot_map[0][row][2] = colat
            #if longitude is chosen
        elif col == 3:
            long = random.uniform(0., 360.) 
            self.spot_map[0][row][3] = long
            #if radius is chosen
        elif col == 4:
            r = random.uniform(0.0, 8.0)
            self.spot_map[0][row][4] = r
        
        while not self.check_params_compatibility(self.spot_map): 
                if col == 0:
                    new_t_ini = random.uniform(self.dt1-self.spots_lifetime, self.dt2)    
                    self.spot_map[0][row][0] = new_t_ini
                #spots lifetime     
                elif col == 1: 
                    new_t_life = abs(random.gauss(self.spots_lifetime, self.sigma_spots_lifetime))
                    self.spot_map[0][row][1] = new_t_life
                #if colatitude is chosen
                elif col == 2:
                    colat = random.uniform(0., 180.)
                    self.spot_map[0][row][2] = colat
                #if longitude is chosen
                elif col == 3:
                    long = random.uniform(0., 360.) 
                    self.spot_map[0][row][3] = long
                #if radius is chosen
                elif col == 4:
                    r = random.uniform(0.0, 8.0)
                    self.spot_map[0][row][4] = r
                    
                """
                new_spot = [random.uniform(self.dt1-self.spots_lifetime, self.dt2), 
                            abs(random.gauss(self.spots_lifetime, self.sigma_spots_lifetime)),
                            random.uniform(0., 360.), 
                            random.uniform(0.0, 180.0),
                            random.uniform(0.0, 8.0) ]
                self.spot_map[0][spot] = new_spot
                """
        return self.spot_map[0]
                                
    def gen_rnd_map(self, num_spots):
        """returns a random spot map
        """
        try:
            #find init time and final time of simulation (dt1 and dt2)
            #take the widest range possible
            self.dt1 = min([min(object.obs_time) for object in self.Observables])  
            self.dt2 = max([max(object.obs_time) for object in self.Observables])  
        
        except:
            pass
         
        data_block = []
        for i in range(int(num_spots)):
            t_init = random.uniform(self.dt1-2.*self.spots_lifetime, self.dt2)#+self.spots_lifetime)
            t_life = abs(random.gauss(self.spots_lifetime, self.sigma_spots_lifetime))
            #generate spot longitude ([0:360] deg)
            x_long = random.uniform(0., 360.)
            #generate spot colatitude ([20:160] deg)
            x_colat = random.uniform(0.0, 180.0)
            #x_colat = random.gauss(90.0, 60.0)
            spot_radius = random.uniform(0.0, 4.0)  
            #append data to output matrix
            data_block.append([t_init, t_life, x_colat, x_long, spot_radius])

        return data_block
    
    def check_params_compatibility(self, spot_params_in):
        """retruns False if a spot parameter array has any pair of spots overlapping (is rotation period dependent!)
        """  
        _spot_map_ = spot_params_in[0]
        _spot_params_ = spot_params_in[1:]  #(evo. rate, diff_rot, n_spots, Q )

        if mstarsim.check_spot_config(_spot_map_, _spot_params_, len(_spot_map_)) == 1:
            return True
        else:
            return False

    def avoid_overlaps(self, spot_params_in):
        """retruns a slightly modified map with no spot overlaps
        """  
        _spot_map_ = spot_params_in[0]
        _spot_params_ = spot_params_in[1:]  #(evo. rate, diff_rot, n_spots, Q )

        is_ok, sp1, sp2 = mstarsim.check_spot_config_w(_spot_map_, _spot_params_, len(_spot_map_))
        
        if not is_ok:
            counter = 0
            while(1):
                counter += 1
                _spot_map_[sp1-1][2] = float(_spot_map_[sp1-1][2])  \
                                                   + random.gauss(0.0,5.0)
                _spot_map_[sp1-1][3] = float(_spot_map_[sp1-1][3])  \
                                                   + random.gauss(0.0,5.0) 
                _spot_map_[sp1-1][4] = float(_spot_map_[sp1-1][4])  \
                                                   - abs(random.gauss(0.0,1.0))
                                               
                is_ok, sp1, sp2 = mstarsim.check_spot_config_w(_spot_map_, _spot_params_, len(_spot_map_))
                
                if is_ok:
                    return _spot_map_ + _spot_params_
                
                elif counter > 10000:
                    print "Warning: Error in avoiding overlaps"
                    return _spot_map_ + _spot_params_
        else:
            return _spot_map_ + _spot_params_
    
    def density_map(self, colat, lon, sign):
        """
        """
        
        
    
    def gen_modified_state(self, old_state, random_walk):
        """return a perturbed spot config
           canonical spot configuration (t_ini, t_life, colatitude, longitude, radii)
        """
        m_state = copy.deepcopy(old_state)
        
        #shuffle spot & parameter to perturbe
        if random_walk:
            #only possible to chanch coordinates & radii
            row = int(random.uniform(0,len(m_state[0])))
            col = int(random.uniform(2,5))
        else:
            row = int(random.uniform(0,len(m_state[0])))
            col = int(random.uniform(0,5))

        #spot init time
        if col == 0:
            new_t_ini = random.uniform(self.dt1 - 2.*self.spots_lifetime, self.dt2)# , self.dt1+self.spots_lifetime) #update the new value    
            m_state[0][row][0] = new_t_ini
        #spots lifetime     
        elif col == 1: 
            new_spot_lifetime = abs(random.gauss(self.spots_lifetime, self.sigma_spots_lifetime))
            m_state[0][row][1] = new_spot_lifetime
        #if colatitude is chosen
        elif col == 2:
            if random_walk:
                new_colat = float(old_state[0][row][2]) + random.gauss(0.0, 0.2)   #throw_rnd_gauss(0.0, 1.0, 0.05)
            else:
                #new_colat = self.throw_rnd_unif(self.colat_limits[0], self.colat_limits[1], 0.1)
                new_colat = self.random_sin(180.0)
            #update the new value    
            m_state[0][row][2] = new_colat
        #if longitude is chosen
        elif col == 3:
            if random_walk:
                new_long = float(old_state[0][row][3]) + random.gauss(0.0, 0.2)    #throw_rnd_gauss(0.0, 1.0, 0.05)
            else:
                new_long = self.throw_rnd_unif(self.long_limits[0], self.long_limits[1], 0.1)
            #update the new value
            m_state[0][row][3] = new_long
        #if radius is chosen
        elif col == 4:
            if random_walk:
                new_r = float(old_state[0][row][4]) + random.gauss(0.0, 0.02)      #throw_rnd_gauss(0.0, 0.05, 0.001)
            else:
                #new_r = self.throw_rnd_unif(r_limits[0], r_limits[1], 0.01)
                new_r = random.uniform(self.r_limits[0], self.r_limits[1])

            while not self.in_range(new_r, self.r_limits):
                new_r = (self.r_limits[1] - self.r_limits[0])/2.0 + random.gauss(0.0, 0.05)  #throw_rnd_gauss(0.0, 0.5, 0.01)
            #update the new value
            m_state[0][row][4] = new_r
            
        return m_state

    def spherical_distance(self, th_0, ph_0, th, ph):
        """
        """
        try:
            to_radians = pi/180.0
            to_deg = 1.0/to_radians
            cosD = sin(pi/2.0 - th_0*to_radians)*sin(pi/2.0 - th*to_radians) + \
                cos(pi/2.0 - th_0*to_radians)*cos(pi/2.0 - th*to_radians)*cos(abs(ph_0 - ph)*to_radians)
            distance = to_deg*acos(cosD)
            return distance
        
        except ValueError:
            
            return 0.0
    
    def throw_rnd_gauss(self, mu, std, step):
        """return gaussian distributed number of (mu,sigma) params in steps of 'step'
        e.g: for i in range(5): print throw_rnd_gauss(0.0,1.0,0.1) --> {-2.6, 1.2, -1.0, 0.3, 0.4}
        """
        nbr = random.gauss(mu, std)
        c_nbr = int((1./step)*nbr)
        r = c_nbr*step
        if r == 0.0:
            return throw_rnd_gauss(mu,std,step)
        else:
            return r

    def throw_rnd_unif(self, inf, sup, step):
        """return uniform distributed number in range (inf,sup) in steps of 'step'
        """
        nbr = random.uniform(inf, sup)
        c_nbr = int((1./step)*nbr)
        r = c_nbr*step

        return r

    def in_range(self, nbr, limits):
        """
        """
        if nbr >= limits[0] and nbr <= limits[1]:
            return True
        else:
            return False

class Planet:
    """packs all information related to the transiting planet
    """
    def __init__(self, conf_file): 
        #init configuration file
        self.conf_file = conf_file
        #------------------
        #planet parameters 
        #------------------ 
        #planet rotation period
        self.P_planet = float(self.conf_file.get('planet','t_planet'))
        #planet T0-ephemerid
        self.T0_planet = float(self.conf_file.get('planet','planet_eph_t0'))
        #parameter of impact
        self.bim = float(self.conf_file.get('planet','planet_impact_param'))
        #spin-orbit angle (deg)
        self.spin_orbit = float(self.conf_file.get('planet','spin_orbit_angle'))
        #planet radius (R* units)
        self.R_planet = float(self.conf_file.get('planet','planet_radius'))
        #planet semiaxis in R_star
        self.A_planet = (((self.m_star*0.05229*1.4248*1.e5)*self.P_planet**2)/(4.0*pi**2))**0.33333   
    
        
    def __del__(self):
        """class destructor
        """   

class CPhotometry(Starsim_IO, Spectra, Star):
    """packs all information related to a single photometry band
    """    
    def __init__(self, SimulationEnvironment, conf_file_path='./starsim.conf', inhibit_load_spectra=False):
        
        #init configuration file
        #self = SimulationEnvironment
        self.conf_file_path = conf_file_path
        self.inhibit_load_spectra = inhibit_load_spectra
        self.conf_file = self.__conf_init(conf_file_path)
        #class composition (to access the attributes of StarsimSimulation objects)
        self.Sim = SimulationEnvironment
        #object id
        self.id = 'Ph_' + str(int(random.uniform(1,2048)))
        self.mean_normalized_series = -1.0
        self.normalized_series = []
        #------------------
        #sim. parameters 
        #------------------
        #grid size
        self.ang_res = float(self.conf_file.get('general','ang_res'))
        #initial time of simulation
        self.dt1 = float(self.conf_file.get('general','init_time_sim'))
        #final time of simulation
        self.dt2 = float(self.conf_file.get('general','final_time_sim'))
        #exposure time cadence
        self.t_exp = float(self.conf_file.get('general','time_cadence'))
        #min spectral wavelength
        self.wvmin = float(self.conf_file.get('general','min_spec_range_gen'))
        #max spectral wavelength
        self.wvmax = float(self.conf_file.get('general','max_spec_range_gen'))
        #cth coeffs
        self.cth_coeffs = self.compute_cth_coeffs()
        #type of data (ph/rv/...)
        self.data_type = ''
        #spectra resolution
        self.spectra_res = str(self.conf_file.get('general','spectra_resolution'))
        #-------------------------------------------------------------------------
        #INIT ROUTINES
        #-------------------------------------------------------------------------
        if not self.inhibit_load_spectra:
            #read generated spectra or create if not found
            self.read_precompiled_spectra()
            #immaculate photosphere or read from disk
            self.compute_immaculate_photosphere()
            self.immaculate_photosphere_flux()
            #filter
            self.gen_filter()
            
            #load time arrays (from .conf file)
            self.load_observational_data(file_path='')
            
            #flux of immaculate photosphere filtered
            self.immaculate_flux_filtered = mstarsim.filter_convolution(self.flt_immaculate_ph, self.interpolated_filter_coeffs, len(self.wv))     
            #integrate immaculate flux
            self.immaculate_flux = mstarsim.trapezoid_integration(self.wv, self.immaculate_flux_filtered, len(self.wv))
        
        #photometry "zero point"
        self.z_ph = 1.0
        #photometric jitter
        self.jitter = 0.0
               
    def __del__(self):
        """class destructor
        """ 
        
    def __conf_init(self, conf_file_path):
        """creates an instance of class ConfigParser, and read the .conf file. Returns the object created
        ConfigParser, and read the .conf file. Returns the object created
        """
        conf_file_Object = ConfigParser.ConfigParser()
        if not conf_file_Object.read([  conf_file_path]):
            print_message( str(conf_file_path) + ": could not be read", 5, 31, inhibit_msg=False)
            sys.exit()
        else:
            return conf_file_Object
    
    def filter_interpol(self, filter_coeffs):
        """returns an array with lambda, interpolated filter coeff
        """
        interpolated_filter_coeffs = np.zeros(len(self.wv), dtype=float)
        for j in range(len(self.wv)):
            for filter_index in range(len(filter_coeffs)-1):
                if float(filter_coeffs[filter_index][0])*1000.0 < self.wv[j] and \
                float(filter_coeffs[filter_index+1][0])*1000.0 > self.wv[j] :
                    interpolated_filter_coeffs[j] = float(filter_coeffs[filter_index][1]) +  \
                            ((float(filter_coeffs[filter_index+1][1]) - \
                            float(filter_coeffs[filter_index][1]))/(float(filter_coeffs[filter_index+1][0])*1000.0 - \
                            float(filter_coeffs[filter_index][0])*1000.0))*(self.wv[j]- float(filter_coeffs[filter_index][0])*1000.0 )
                    break

        return interpolated_filter_coeffs

    def gen_filter(self):
        """
        """
        try:
            try:
                test = self.filter_name
            except AttributeError:     
                self.filter_name = self.conf_file.get('files','filter_path').split('/')[-1]
            
            if self.spectra_res == 'high':
                file = './data/interpolated_filter_coeffs_' + str(self.filter_name) + '_' + str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '.dat'
            elif self.spectra_res == 'low':
                file = './data/interpolated_filter_coeffs_' + str(self.filter_name) + '_' + str(int(self.wvmin)) + '_' + str(int(self.wvmax)) + '_lres.dat'
           
            self.interpolated_filter_coeffs = self.read_file(file, ' ')
            self.print_message('    Reading interpolated filter coefficients...', index=3, color=32)
        
        except IOError:
            filter_coeffs = self.read_file(self.conf_file.get('files','filter_path'), ' ')
            print '    Interpolating filter coefficients...'
            self.interpolated_filter_coeffs = self.filter_interpol(filter_coeffs)
            self.print_message( '    Writing file...' + file, index=3, color=31)
            self.write_file_onecol(file, self.interpolated_filter_coeffs, ' ', '')
            
        return self.filter_name, self.interpolated_filter_coeffs
               
    def logL_null_model(self):
        """return logL of data, centered to its mean
        """
        def compute_logL_null(jitter_null):
            """
            """
            inv_sigma2_null = 1.0/(self.err_0**2 + jitter_null**2)
            logL_null = -0.5*(np.sum(((self.y_data_0 - np.mean(self.y_data_0))**2)*inv_sigma2_null + np.log(2.*np.pi) - np.log(inv_sigma2_null))) 
         
            return -logL_null 
        
        
        init_jitter_0 = 0.01
        res = minimize(compute_logL_null, init_jitter_0, method='Nelder-Mead', tol=1e-6)
       
        self.jitter_null = res.x
        self.logL_null = -compute_logL_null(self.jitter_null)
           
        return self.logL_null, self.jitter_null
    
    def model_logL(self):
        """compute likelihood of model-data
        """
        inv_sigma2 = 1.0/(self.err**2 + self.jitter**2)
       
        self.logL = -0.5*(np.sum(((self.compute_lc_model() - self.y_data)**2)*inv_sigma2 + np.log(2.*np.pi) - np.log(inv_sigma2))) 
                    
        return self.logL
    
    def compute_lc_model(self):
        """
        Compute the lightcurve given a Simulation.spotmap 
        Simulation: an object containing the info of the simulation and where spotmap is defined)
        write_lc: write lc on disk
        """ 
        def normalize_time_series(tflux_series, normalization):
            """filter convolution and normalization over immaculate photosphere flux
            """
            flux_series = [ [tflux_series[i][0], tflux_series[i][1]/normalization] for i in range(len(tflux_series))]
                
            return flux_series
              
        #compute rotating spotted photosphere
        raw_flux_series, S_index, ff_sp, ff_fc = self.compute_rotating_spotted_photosphere()
        #normalize the series with respect to immaculate flux
        normalized_series = normalize_time_series(raw_flux_series, self.immaculate_flux)
        
        #return only y-values
        return np.array(normalized_series)[:,1] 
    
    def compute_lc_series(self, write_lc=False, mark=None):
        """
        Compute the lightcurve given a Simulation.spotmap 
        Simulation: an object containing the info of the simulation and where spotmap is defined)
        write_lc: write lc on disk
        """ 
        def normalize_time_series(tflux_series, normalization):
            """filter convolution and normalization over immaculate photosphere flux
            """
            flux_series = [ [tflux_series[i][0], tflux_series[i][1]/normalization] for i in range(len(tflux_series))]
                
            return flux_series
              
        #compute rotating spotted photosphere
        if not Glob.planet:
            raw_flux_series, S_index, ff_sp, ff_fc = self.compute_rotating_spotted_photosphere()
        else:
            raw_flux_series = self.compute_rotating_spotted_photosphere()
        
        #normalize the series with respect to immaculate flux
        self.normalized_series = normalize_time_series(raw_flux_series, self.immaculate_flux)
        #mean of ph series
        normalized_series_ = np.array(self.normalized_series)
        self.mean_normalized_series = normalized_series_[:,1].astype(float).mean()
        
        if self.obs_data != []:
            logL_lc = self.model_logL()
              
        else:
            logL_lc = 'not defined'
       
        #write lc file?
        if write_lc:
            if mark == None: 
                mark = ''
            else:
                mark = mark + '_'
                
            if self.obs_data != []:
                lc_model_residuals = []
                for i in range(len(self.obs_time)):
                    lc_model_residuals.append( [ self.obs_time[i], self.normalized_series[i][1], self.obs_data[i][2], self.normalized_series[i][1] - self.obs_data[i][1] ] ) 
                
                
                self.write_file('./output/lc_series_' + str(mark) + str(self.Sim.id) + '_' + str(self.wvmin) + '_' + str(self.wvmax) + '_' + str(self.filter_name) 
                                + '.dat', lc_model_residuals, ' ', '')
            
            else:
                
                self.write_file('./output/lc_series_' + str(mark) + str(self.Sim.id) + '_' + str(self.wvmin) +'_' + str(self.wvmax) + '_' + str(self.filter_name) 
                                + '.dat', self.normalized_series, ' ', '')
            
            
            self.print_message( str(4*' ') + './output/lc_series_' + str(mark) + str(self.Sim.id) + '_' + str(self.wvmin) + '_' + str(self.wvmax) + 
                                '_' + str(self.filter_name) + '.dat' + " created", 2, 31)        
       
            #S-index
            """
            normalized_series_S_index = normalize_time_series(S_index, self.S_index_im)
            
            self.write_file('./output/S-index_series_' + str(self.wvmin) +'_' + str(self.wvmax) + '.dat', 
                        normalized_series_S_index, ' ', '')
            """
            #Filling factors
            if not Glob.planet:
                try:
                    self.write_file('./output/spot_ff_' + str(self.wvmin) +'_' + str(self.wvmax) + '.dat', 
                                ff_sp, ' ', '')
                    self.write_file('./output/faculae_ff_' + str(self.wvmin) +'_' + str(self.wvmax) + '.dat', 
                                ff_fc, ' ', '')
                except:
                    print "Error writing filling factor curves"
                    pass
            
            
        return logL_lc

    def compute_rotating_spotted_photosphere(self):
        """input: Band 
            a StarsimSimulation object, containing all properties and methods
        """  
        #call external library
        if Glob.planet:
            
            self.print_message("Computing light curve + transits...", 2, 31)
            
            self.Sim.w = 2.0*pi/self.Sim.P_rot
            star_params = [self.Sim.x_i, self.Sim.R_star, self.Sim.w, self.Sim.B_coeff, self.Sim.C_coeff, 
                           self.Sim.temp_ph, self.Sim.delta_t_sp, self.Sim.delta_t_fc]
            planet_params = [self.Sim.T0_planet, self.Sim.P_planet, self.Sim.A_planet, self.Sim.spin_orbit, self.Sim.bim, self.Sim.R_planet]
            star_params = star_params + planet_params
            
            flt_bolometric_arr = mstarsim.compute_rotating_spotted_photosphere_planet(self.obs_time, len(self.obs_time), \
                                                                                    self.ang_res, star_params, self.Sim.spot_map[0], \
                                                                                    self.Sim.spot_map[1:], self.cth_coeffs, self.interpolated_filter_coeffs, \
                                                                                    self.wv, self.dlnp_ph, self.flpk_ph, self.flnp_ph, \
                                                                                    self.dlnp_sp, self.dlnp_fc, self.flpk_sp, self.flnp_sp, \
                                                                                    self.flnp_fc, self.flt_immaculate_ph, len(self.wv), len(self.Sim.spot_map[0]))
            return flt_bolometric_arr
        
        else:
     
            self.Sim.w = 2.0*pi/self.Sim.P_rot
            star_params = [self.Sim.x_i, self.Sim.R_star, self.Sim.w, self.Sim.B_coeff, self.Sim.C_coeff, 
                           self.Sim.temp_ph, self.Sim.delta_t_sp, self.Sim.delta_t_fc]
        
            flt_bolometric_arr, S_index, ff_sp, ff_fc = mstarsim.compute_rotating_spotted_photosphere(self.obs_time, len(self.obs_time), 
                                                                            star_params, self.Sim.spot_map[0], 
                                                                            self.Sim.spot_map[1:], self.cth_coeffs,  
                                                                            self.interpolated_filter_coeffs, 
                                                                            self.wv, self.dlnp_ph, self.flpk_ph, self.flnp_ph, 
                                                                            self.dlnp_sp, self.flpk_sp, self.flnp_sp, 
                                                                            self.flnp_fc, self.flt_immaculate_ph, self.S_index_im, len(self.wv), len(self.Sim.spot_map[0]))
        
            return flt_bolometric_arr, S_index, ff_sp, ff_fc

    def set_wv_range(self, wvmin, wvmax):
        """setter for changing the wvmin,wvmax and read the new spectras, immac. photosphere and filter 
        """
        self.wvmin, self.wvmax = wvmin, wvmax
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        self.immaculate_photosphere_flux()
        #filter
        self.gen_filter()
        #flux of immaculate photosphere filtered
        self.immaculate_flux_filtered = mstarsim.filter_convolution(self.flt_immaculate_ph, 
                                                                    self.interpolated_filter_coeffs, len(self.wv))     
        #integrate immaculate flux
        self.immaculate_flux = mstarsim.trapezoid_integration(self.wv, self.immaculate_flux_filtered, len(self.wv))
      
    def set_delta_t(self):
        """
        """
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #filter
        self.gen_filter()
        #flux of immaculate photosphere filtered
        self.immaculate_flux_filtered = mstarsim.filter_convolution(self.flt_immaculate_ph, 
                                                                    self.interpolated_filter_coeffs, len(self.wv))     
        #integrate immaculate flux
        self.immaculate_flux = mstarsim.trapezoid_integration(self.wv, self.immaculate_flux_filtered, len(self.wv))
        
    def set_ph_zero_point(self, ph_zero_point):
        """offset the lc by summation of a constant in all of data points
        """
        path = self.info
        self.load_observational_data(path)
        self.z_ph = ph_zero_point 
        self.err = self.err_0 / self.z_ph
        for j in range(len(self.obs_data)):
            self.obs_data[j][1] /= self.z_ph
            self.y_data[j] /= self.z_ph
            try:
                self.obs_data[j][2] /= self.z_ph   
            except:
                raise
            
    def write_lc_data(self):
        """write the lc data curve with z_ph offset and jitter
        """
        self.obs_data, self.err = np.array(self.obs_data), np.array(self.err)
        self.obs_data[:,2] = np.sqrt(self.err**2 + self.jitter**2)
    
        self.write_file('./output/lc_data_' + str(self.wvmin) +'_' + str(self.wvmax) + '_' + str(self.filter_name) + 
                        '_' + str(self.z_ph) + '_' + str(self.jitter), self.obs_data, ' ', '')
        self.print_message( str(4*' ') + "./output/lc_data_" + str(self.wvmin) +'_' + str(self.wvmax) + '_' + 
                            str(self.filter_name) + '_' + str(self.z_ph) + '_' + str(self.jitter) + " created", 2, 31) 
    
class CRadialVelocity(StarsimCommon, Starsim_IO, Spectra, Star):
    """
    """
    def __init__(self, SimulationEnvironment, conf_file_path='./starsim.conf', 
                 inhibit_load_spectra=False):
        
        warnings.filterwarnings('ignore')
        #init configuration file
        self.conf_file_path = conf_file_path
        self.inhibit_load_spectra = inhibit_load_spectra
        self.conf_file = self.__conf_init(conf_file_path)
        #class composition (to access the attributes of StarsimSimulation objects)
        self.Sim = SimulationEnvironment
        #object id
        self.id = 'RV_' + str(int(random.uniform(1,2048)))
        self.id_n = 0
        #-------------------------------------------------------------------------
        #sim. parameters 
        #-------------------------------------------------------------------------
        #grid size
        self.ang_res = float(self.conf_file.get('general','ang_res'))
        #initial time of simulation
        self.dt1 = float(self.conf_file.get('general','init_time_sim'))
        #final time of simulation
        self.dt2 = float(self.conf_file.get('general','final_time_sim'))
        #exposure time cadence
        self.t_exp = float(self.conf_file.get('general','time_cadence'))
        #min spectral wavelength
        self.wvmin = float(self.conf_file.get('general','min_spec_range_gen'))
        #max spectral wavelength
        self.wvmax = float(self.conf_file.get('general','max_spec_range_gen'))
        #cth coeffs
        self.cth_coeffs = self.compute_cth_coeffs()
        #type of data (ph/rv/...)
        self.data_type = ''
        #radial velocity CCF range 
        self.ccf_vel_range = float(self.conf_file.get('rv','ccf_vel_range'))
        #spectra resolution
        self.spectra_res = str(self.conf_file.get('general','spectra_resolution'))
        #type of data (ph/rv/...)
        self.data_type = ''
        #mask
        self.mask = str(self.conf_file.get('rv','spectral_mask'))
        self.normalized_series = []
        #-------------------------------------------------------------------------
        #INIT ROUTINES
        #-------------------------------------------------------------------------
        #load time arrays (from .conf file)
        if not self.inhibit_load_spectra:
            #read generated spectra or create if not found
            self.read_precompiled_spectra()
            #immaculate photosphere or read from disk
            self.compute_immaculate_photosphere()
            self.load_observational_data(file_path='')
            #init CCFs
            self.init_ccfs()
                 
        #init point of CCF fitting (A, center, sigma, offset)
        self.init_p = (30.0,-0.5,1.0,-5.0)
        #RV jitter
        self.jitter = 0.0
        
    def __del__(self):
        """class destructor
        """ 
        
    def __conf_init(self, conf_file_path):
        """creates an instance of class ConfigParser, and read the .conf file. Returns the object created
        ConfigParser, and read the .conf file. Returns the object created
        """
        conf_file_Object = ConfigParser.ConfigParser()
        if not conf_file_Object.read([  conf_file_path]):
            print_message( str(conf_file_path) + ": could not be read", 5, 31, inhibit_msg=False)
            sys.exit()
        else:
            return conf_file_Object

    def ccf_indices(self, time_array, rv, ccf_rotating):
        """calculate RV, BIS, FWHM, CONTRAST
        """
        indices_array = []
        
        for i in range(len(time_array)):
            y = ccf_rotating[i][:]
            #fit gaussian
            coeffs, var_matrix = curve_fit(self.gauss, rv, y, p0=self.init_p)
            self.init_p = coeffs
            amplitude, center, sigma, offset = coeffs
            FWHM = 2.0*sigma*sqrt(2.0*log(2))
            CONTRAST = 0.0
            #bisector 
            bisector_low = np.array( self.compute_bisector(rv,y,lims=[0.10,0.40]) )
            bisector_high = np.array( self.compute_bisector(rv,y,lims=[0.60,0.90]) )
            BIS = -(bisector_high[:,0].astype(float).mean() - bisector_low[:,0].astype(float).mean()) 
            #write_file('./bisector_' + str(i), bisector, ' ', '' )
            indices_array.append([center, BIS, FWHM, CONTRAST])     
        
        return time_array, indices_array 
        
    def compute_bisectors(self, coeffs, rv, ccfs, lims):
        """
        """
        N_BIS_POINTS = 80
        bisector_series = []
        
        for i in range(len(ccfs)):
            ccf = ccfs[i]
            a, mean, sigma, offset = coeffs[i,0], coeffs[i,1], coeffs[i,2], coeffs[i,3]  
            rv_bis = list()
            y_ccf_l = list()
            bisector = []
            
            max_gaussian_fit = self.gauss(mean, a, mean, sigma, offset)
            MAX_BIS, MIN_BIS = lims[1]*(max_gaussian_fit - offset) , lims[0]*(max_gaussian_fit - offset) 
          
            for i in range(N_BIS_POINTS):
                y_ccf = MIN_BIS + i*(MAX_BIS - MIN_BIS) / N_BIS_POINTS
                for j in range(len(ccf)-1):
                    if (ccf[j] < y_ccf) and (ccf[j+1] > y_ccf) and (rv[j] < mean):
                        interpolated_rv_1 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))
                  
                    elif (ccf[j] > y_ccf) and (ccf[j+1] < y_ccf) and (rv[j] > mean):
                        interpolated_rv_2 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))  
                
                rv_bis.append((interpolated_rv_2 + interpolated_rv_1)/2.0) 
                y_ccf_l.append(y_ccf)
                bisector.append([(interpolated_rv_2 + interpolated_rv_1)/2.0, y_ccf])
            
            bisector_series.append(bisector)   
        
        return np.array(bisector_series)
    
    def set_wv_range(self, wvmin, wvmax):
        """setter for changing the wvmin,wvmax and read the new spectras, immac. photosphere and ccfs 
        """
        self.wvmin, self.wvmax = wvmin, wvmax
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #init CCFs
        self.init_ccfs()

    def set_delta_t(self):
        """
        """
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #init CCFs
        self.init_ccfs()

    def systemic_vel(self, rv_sym, ccf_im):
        """fit a gaussian to CCF immaculate photosphere
        """     
        coeff, var_matrix = curve_fit(self.gauss, rv_sym, ccf_im, p0=(50.0,-0.5,-1.0,-5.0), ftol=1e-09, xtol=1e-09)
        self.sys_vel = coeff[1]     
            
        return self.sys_vel

    def logL_null_model(self):
        """return logL of data, centered to its mean
        """
        def compute_logL_null(jitter_null):
            """
            """
            inv_sigma2_null = 1.0/(self.err_0**2 + jitter_null**2)
            logL_null = -0.5*(np.sum(((self.y_data_0 - np.mean(self.y_data_0))**2)*inv_sigma2_null + np.log(2.*np.pi) - np.log(inv_sigma2_null))) 
         
            return -logL_null 
        
        init_jitter_0 = 0.01
        res = minimize(compute_logL_null, init_jitter_0, method='Nelder-Mead', tol=1e-6)
       
        self.jitter_null = res.x
        self.logL_null = -compute_logL_null(self.jitter_null)
           
        return self.logL_null, self.jitter_null
           
    def model_logL(self):
        """compute likelihood of model-data
        """
        #RV model (only y values)
        y_rv_model = self.compute_rv_model()[:,1].astype(float)
        #renormalize the error column
        errs = self.err / max(self.err)
        inv_sigma2 = 1.0/(errs*self.jitter**2)
        #inv_sigma2 = 1.0/(self.err**2 + self.jitter**2)
        #inv_sigma2 = 1.0/(self.jitter**2)
        self.logL = -0.5*(np.sum(((y_rv_model - self.y_data)**2)*inv_sigma2 + np.log(2.*np.pi) - np.log(inv_sigma2))) 
                    
        return self.logL
    
    def fit_ccf_bisector(self, rv, ccf):
        """returns bis and a set of coefficients that fits 10-90% interval bisector of CCF
        """
        N_BIS_POINTS = 80
        rv_bis = list()
        y_ccf_l = list()
        bisector = []
        #get the coefficients of a gaussian fit gaussian_fit=(a,mean,sigma,offset)
        """
        gaussian_fit = fit_gaussian(rv, ccf)
        """
        #find the gaussian max.
        coeffs, mvar = curve_fit(self.gauss, rv, ccf, p0=(50.0,-0.5,-1.0,-5.0), ftol=1e-09, xtol=1e-09)
        a, mean, sigma, offset = coeffs[:]
        
        max_gaussian_fit = self.gauss(mean, a, mean, sigma, offset)
        #limits of CCF bisector
        #0.1438
        MAX_BIS, MIN_BIS = 0.8*(max_gaussian_fit - offset) , 0.2*(max_gaussian_fit - offset) 
        
        for i in range(N_BIS_POINTS):
            y_ccf = MIN_BIS + i*(MAX_BIS - MIN_BIS) / N_BIS_POINTS
            for j in range(len(ccf)-1):
                if ccf[j] < y_ccf and ccf[j+1] > y_ccf and rv[j] < mean:
                    interpolated_rv_1 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))
                elif ccf[j] > y_ccf and ccf[j+1] < y_ccf and rv[j] > mean:
                    interpolated_rv_2 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))  
            rv_bis.append((interpolated_rv_2 + interpolated_rv_1)/2.0) 
            y_ccf_l.append(y_ccf)
            bisector.append([(interpolated_rv_2 + interpolated_rv_1)/2.0, y_ccf])
        #fit the bisector to a 6-deg polinomial        
        bis_coeffs = np.polyfit(np.array(y_ccf_l), np.array(rv_bis), deg=6, rcond=None, full=False, w=None, cov=False)   
        """
        for i in range(len(bisector)):
            print bisector[i][0], bisector[i][1]
        """
        return bisector, bis_coeffs

    def gen_sym_ccf(self, rv, ccf, bis_coeffs):
        """
        """
        rv_sym = np.zeros(len(rv[:]))

        for j in range(len(ccf)):
            rv_sym[j] = rv[j] - self.impoly(ccf[j], bis_coeffs) 

        return rv_sym, ccf

    def fit_ccf(self, time_array, rv, ccf_rotating):
        """
        """ 
        coeffs_array = []
        #fit a gaussian in each ccf
        for i in range(len(time_array)):
            y = ccf_rotating[i][:]
            try:
                coeff, var_matrix = curve_fit(self.gauss, rv, y, p0=self.init_p, ftol=1e-4, xtol=1e-4)
                #coeff = [amplitude, center, sigma, offset]
                coeffs_array.append(coeff)
                self.init_p = coeff  #reset the initial guess with the last fit output
            
            except RuntimeError:
                print "Error in fitting CCF"
                raise
     
        return time_array, np.array(coeffs_array)
    
    def compute_rv_series(self, write_rv=True, mark=None):
        """
        """        
        normalized_rv_series = self.compute_rv_model()
       
        #compute logL, if observational data is provided
        if self.obs_data != []:
            logL_rv = self.model_logL()
        else:
            logL_rv = 'not defined'
        
        if write_rv:
            #write RV series&residuals on disk
            #compute the residuals
            if mark == None: 
                mark = ''
                if Glob.indices: mark = mark + 'indices_'
            elif Glob.indices:
                mark = mark + '_indices_'
            else:
                mark = mark + '_'
                
            if self.obs_data != []:
                #if observational data exists, compute residuals as well 
                rv_series = [ ['# time, rv(model), residuals(data-model), rv_err(data)'] ]
                for i in range(len(self.obs_time)):
                    rv_series.append( [ self.obs_time[i], self.normalized_series[i][1], 
                                                 self.normalized_series[i][1] - self.obs_data[i][1], self.obs_data[i][2] ] ) 
               
                self.write_file('./output/rv_series_' + str(mark) + str(self.Sim.id) + '_' + str(self.wvmin) +'_' + str(self.wvmax) + '_' + 
                                str(self.mask) + '.dat', rv_series, ' ', '')
            
            else:
                
                self.write_file('./output/rv_series_' + str(mark) + str(self.Sim.id) + '_' + str(self.wvmin) +'_' + str(self.wvmax) + '_' + 
                                str(self.mask) + '.dat', self.normalized_series, ' ', '')
            
            
            self.print_message( str(4*' ') + './output/rv_series_' + str(mark) + str(self.Sim.id) + '_' + 
                                str(self.wvmin) + '_' + str(self.wvmax) + '_' + str(self.mask) + '.dat' + " created", 2, 31)        
            
                   
        return logL_rv

    def compute_rv_model(self):
        """
        Compute the radial velocity curve given a Simulation.spotmap 
        Simulation: an object containing the info of the simulation and where spotmap is defined)
        write_rv: write lc on disk
        """ 
        self.Sim.w = 2.0*pi/self.Sim.P_rot
        star_params = [self.Sim.x_i, self.Sim.R_star, self.Sim.w, self.Sim.B_coeff, self.Sim.C_coeff, 
                       self.Sim.temp_ph, self.Sim.delta_t_sp, self.Sim.delta_t_fc]
        planet_params = [self.Sim.T0_planet, self.Sim.P_planet, self.Sim.A_planet, self.Sim.spin_orbit, self.Sim.bim, self.Sim.R_planet]
        
        star_params = star_params + planet_params
        
        spot_map = self.Sim.spot_map[:]
        spot_params_, spot_map, n_spots = spot_map[1:], spot_map[0], self.Sim.n_spots #spot_map[3]
        
        #self.print_message('Computing CCF spotted rotating photosphere', index=6, color=34)
        #rotating spotted photosphere ccf 
        ccf_rotating = mstarsim.ccf_rotating_spotted_photosphere(self.obs_time, star_params, 
                                                                  spot_params_, self.cth_coeffs, spot_map, self.dlnp_ph, 
                                                                  self.flpk_ph, self.flnp_ph, self.dlnp_sp, self.flpk_sp, self.flnp_sp, 
                                                                  self.flnp_fc, self.rv_sym, self.ccf_sym_ph, self.ccf_sym_sp, self.ccf_sym_fc, 
                                                                  self.ccf_im, self.flxph, len(self.obs_time), len(self.rv_sym), len(self.flnp_ph), n_spots) 
        
        if Glob.save_ccfs: 
            #save ccf series
            for i in range(len(ccf_rotating)):
                ccf = []
                for j in range(len(ccf_rotating[i])):
                    ccf.append( [ self.rv_sym[j], ccf_rotating[i][j] ] )
                self.write_file('./ccfs/' + str(self.obs_time[i]) + '.dat', ccf, ' ', '')    
            
 
        time_, gaussian_fit_coeffs = self.fit_ccf(self.obs_time, self.rv_sym, ccf_rotating)
        
        #fit a gaussian and get its center (rv)
        rv_series = gaussian_fit_coeffs[:,1].astype(float)
        
        #if "--indices" option set, compute CCF-derived indices
        if Glob.indices:
            #amplitude, center, sigma, offset = gaussian_fit_coeffs
            amplitude, sigma, offset = gaussian_fit_coeffs[:,0], gaussian_fit_coeffs[:,2], gaussian_fit_coeffs[:,3]
            #index contrast
            CONTRAST = 100.0 - (amplitude - offset)
            #index FWHM
            FWHM = 2.0*sigma*sqrt(2.0*log(2))
            #BIS index
            bisectors_low = self.compute_bisectors(gaussian_fit_coeffs, self.rv_sym, ccf_rotating, lims=[0.10,0.40]) 
            bisectors_high = self.compute_bisectors(gaussian_fit_coeffs, self.rv_sym, ccf_rotating, lims=[0.60,0.90])
            BIS = []
            for i in range(len(bisectors_low)):
                BIS.append(-(bisectors_high[i][:,0].astype(float).mean() - bisectors_low[i][:,0].astype(float).mean()) )
            #normalize and save the RV + indices curves       
            self.normalized_series = []
            for i in range(len(time_)): 
                self.normalized_series.append( [ time_[i], 1000.0*(rv_series[i] - self.sys_vel), BIS[i], FWHM[i], CONTRAST[i] ]) 
            
        else:
            #normalize and save the RV curve       
            self.normalized_series = []
            for i in range(len(time_)): self.normalized_series.append([time_[i], 1000.0*(rv_series[i] - self.sys_vel)])
        
        #only y-values
        #self.y_rv_model = np.array(self.normalized_series)[:,1]
        
        return np.array(self.normalized_series)
   
    def write_rv_data(self):
        """write the lc data curve with z_ph offset and jitter
        """
        self.obs_data, self.err = np.array(self.obs_data), np.array(self.err)
        self.obs_data[:,2] = np.sqrt(self.err**2 + self.jitter**2)
    
        self.write_file('./output/RV_data_' + str(self.wvmin) +'_' + str(self.wvmax) + '_' + 
                        str(self.jitter) + '.dat', self.obs_data , ' ', '')
        self.print_message( str(4*' ') + "./output/RV_data_" + str(self.wvmin) + '_' + 
                            str(self.wvmax) + '_' + str(self.jitter) + '.dat' + " created", 2, 31) 
       
class CFWHM(StarsimCommon, Starsim_IO, Spectra, Star):
    """
    """
    def __init__(self, SimulationEnvironment, conf_file_path='./starsim.conf'):
        
        #init configuration file
        self.conf_file_path = conf_file_path
        self.conf_file = self.__conf_init(conf_file_path)
        #class composition (to access the attributes of StarsimSimulation objects)
        self.Sim = SimulationEnvironment
        #object id
        self.id = 'FWHM_' + str(int(random.uniform(1,4096)))
        #-------------------------------------------------------------------------
        #sim. parameters 
        #-------------------------------------------------------------------------
        #grid size
        self.ang_res = float(self.conf_file.get('general','ang_res'))
        #initial time of simulation
        self.dt1 = float(self.conf_file.get('general','init_time_sim'))
        #final time of simulation
        self.dt2 = float(self.conf_file.get('general','final_time_sim'))
        #exposure time cadence
        self.t_exp = float(self.conf_file.get('general','time_cadence'))
        #min spectral wavelength
        self.wvmin = float(self.conf_file.get('general','min_spec_range_gen'))
        #max spectral wavelength
        self.wvmax = float(self.conf_file.get('general','max_spec_range_gen'))
        #cth coeffs
        self.cth_coeffs = self.compute_cth_coeffs()
        #type of data (ph/rv/...)
        self.data_type = ''
        #radial velocity CCF range 
        self.ccf_vel_range = float(self.conf_file.get('rv','ccf_vel_range'))
        #spectra resolution
        self.spectra_res = str(self.conf_file.get('general','spectra_resolution'))
        #type of data (ph/rv/...)
        self.data_type = ''
        #-------------------------------------------------------------------------
        #INIT ROUTINES
        #-------------------------------------------------------------------------
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #load time arrays (from .conf file)
        self.load_observational_data(file_path='')
        #init CCFs
        self.init_ccfs()
        #get the normalization factor for the rms
        self.khi2_normalization_factor()
        
        
    def __del__(self):
        """class destructor
        """ 
        
    def __conf_init(self, conf_file_path):
        """creates an instance of class ConfigParser, and read the .conf file. Returns the object created
        ConfigParser, and read the .conf file. Returns the object created
        """
        conf_file_Object = ConfigParser.ConfigParser()
        if not conf_file_Object.read([  conf_file_path]):
            print_message( str(conf_file_path) + ": could not be read", 5, 31, inhibit_msg=False)
            sys.exit()
        else:
            return conf_file_Object
          
    def set_wv_range(self, wvmin, wvmax):
        """setter for changing the wvmin,wvmax and read the new spectras, immac. photosphere and ccfs 
        """
        self.wvmin, self.wvmax = wvmin, wvmax
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #init CCFs
        self.init_ccfs()

    def set_delta_t(self):
        """
        """
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #init CCFs
        self.init_ccfs()

    def systemic_vel(self, rv_sym, ccf_im):
        """fit a gaussian to CCF immaculate photosphere
        """     
        coeff, var_matrix = curve_fit(self.gauss, rv_sym, ccf_im, p0=(0.5,0.0,1.0,0.0))
        self.sys_vel = coeff[1]     
            
        return self.sys_vel
    
    def khi2_normalization_factor(self):
        """computes the khi2 of data in case of 
        """
        #time_, ccf_indices = self.ccf_indices([0.0], self.rv_sym, [self.ccf_sym_ph])
        immaculate_ph_FWHM = 0.0#ccf_indices[0][1]   
        if self.obs_data != []:
            immaculate_fwhm_im = []
            for i in range(len(self.obs_data)):
                time = self.obs_data[i][0]
                immaculate_fwhm_im.append([time, immaculate_ph_FWHM])
            self.khi2_immaculate_fwhm = mstarsim.khi2(immaculate_fwhm_im, self.obs_data, len(immaculate_fwhm_im))
        else:
            self.khi2_immaculate_fwhm = None

    def fit_ccf_bisector(self, rv, ccf):
        """returns bis and a set of coefficients that fits 10-90% interval bisector of CCF
        """
        N_BIS_POINTS = 100
        rv_bis = list()
        y_ccf_l = list()
        bisector = []
        #get the coefficients of a gaussian fit gaussian_fit=(a,mean,sigma,offset)
        """
        gaussian_fit = fit_gaussian(rv, ccf)
        """
        #find the gaussian max.
        coeffs, mvar = curve_fit(self.gauss, rv, ccf, p0=(0.5,0.0,1.0,0.0))
        a, mean, sigma, offset = coeffs[:]
        
        max_gaussian_fit = self.gauss(mean, a, mean, sigma, offset)
        #limits of CCF bisector
        #0.1438
        MAX_BIS, MIN_BIS = 0.8*(max_gaussian_fit - offset) , 0.2*(max_gaussian_fit - offset) 
        
        for i in range(N_BIS_POINTS):
            y_ccf = MIN_BIS + i*(MAX_BIS - MIN_BIS) / N_BIS_POINTS
            for j in range(len(ccf)-1):
                if ccf[j] < y_ccf and ccf[j+1] > y_ccf and rv[j] < mean:
                    interpolated_rv_1 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))
                elif ccf[j] > y_ccf and ccf[j+1] < y_ccf and rv[j] > mean:
                    interpolated_rv_2 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))  
            rv_bis.append((interpolated_rv_2 + interpolated_rv_1)/2.0) 
            y_ccf_l.append(y_ccf)
            bisector.append([(interpolated_rv_2 + interpolated_rv_1)/2.0, y_ccf])
        #fit the bisector to a 6-deg polinomial        
        bis_coeffs = np.polyfit(np.array(y_ccf_l), np.array(rv_bis), deg=6, rcond=None, full=False, w=None, cov=False)   
        
        return bisector, bis_coeffs

    def gen_sym_ccf(self, rv, ccf, bis_coeffs):
        """
        """
        rv_sym = rv[:]
        for j in range(len(ccf)):
            rv_sym[j] = rv[j] - self.impoly(ccf[j], bis_coeffs) 
            
        return rv_sym, ccf

    def fit_ccf(self, time_array, rv, ccf_rotating):
        """
        """ 
        #initial guess (A,c,sigma,offset)
        p0 = (1.0,0.0,1.0,1.0)
        coeffs_array = np.zeros(len(time_array), dtype=np.float)
        #fit a gaussian in each ccf
        for i in range(len(time_array)):
            y = ccf_rotating[i][:]
            coeff, var_matrix = curve_fit(self.gauss, rv, y, p0, ftol=1e-05)
            coeffs_array[i] = coeff[1]
            p0 = coeff  #reset the initial guess with the last fit output
       
        return time_array, coeffs_array
   
    def compute_bisector(self, rv,ccf, lims):
        """lims = [min_%, max%]
        """
        N_BIS_POINTS = 50
        rv_bis = list()
        y_ccf_l = list()
        bisector = []
        
        coeffs, mvar = curve_fit(self.gauss, rv, ccf, p0=(0.5,0.0,1.0,0.0))
        a, mean, sigma, offset = coeffs[:]
        
        max_gaussian_fit = self.gauss(mean, a, mean, sigma, offset)
        #limits of CCF bisector
        #0.1438
        MAX_BIS, MIN_BIS = lims[1]*(max_gaussian_fit - offset) , lims[0]*(max_gaussian_fit - offset) 
        
        for i in range(N_BIS_POINTS):
            y_ccf = MIN_BIS + i*(MAX_BIS - MIN_BIS) / N_BIS_POINTS
            for j in range(len(ccf)-1):
                if ccf[j] < y_ccf and ccf[j+1] > y_ccf and rv[j] < mean:
                    interpolated_rv_1 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))
                elif ccf[j] > y_ccf and ccf[j+1] < y_ccf and rv[j] > mean:
                    interpolated_rv_2 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))  
            rv_bis.append((interpolated_rv_2 + interpolated_rv_1)/2.0) 
            y_ccf_l.append(y_ccf)
            bisector.append([(interpolated_rv_2 + interpolated_rv_1)/2.0, y_ccf])
            
        return bisector
    
    def ccf_indices(self, time_array, rv, ccf_rotating):
        """calculate RV, FWHM, BIS, CONTRAST
        """
        indices_array = []
        p0=(0.5,0.0,1.0,0.0)
        
        for i in range(len(time_array)):
            y = ccf_rotating[i][:]
            #fit gaussian
            coeffs, var_matrix = curve_fit(self.gauss, rv, y, p0)
            amplitude, center, sigma, offset = coeffs
            FWHM = 2.0*sigma*sqrt(2.0*log(2))
            CONTRAST = 1.0/(amplitude - offset)
            #bisector 
            #bisector_low = np.array( self.compute_bisector(rv,y,lims=[0.10,0.40]) )
            #bisector_high = np.array( self.compute_bisector(rv,y,lims=[0.60,0.90]) )
            #BIS = -(bisector_high[:,0].astype(float).mean() - bisector_low[:,0].astype(float).mean()) 
            BIS = 0.0
            #write_file('./bisector_' + str(i), bisector, ' ', '' )
            indices_array.append([center, FWHM, CONTRAST, BIS])     
            p0 = coeffs
        
        return time_array, indices_array
    
    def compute_FWHM_series(self, write_rv=True):
        """
        """    
        self.Sim.w = 2.0*pi/self.Sim.P_rot
        star_params = [self.Sim.x_i, self.Sim.R_star, self.Sim.w, self.Sim.B_coeff, self.Sim.C_coeff, 
                       self.Sim.temp_ph, self.Sim.delta_t_sp, self.Sim.delta_t_fc]
        planet_params = [self.Sim.T0_planet, self.Sim.P_planet, self.Sim.A_planet, self.Sim.spin_orbit, self.Sim.bim, self.Sim.R_planet]
        
        star_params = star_params + planet_params
        
        spot_map = self.Sim.spot_map[:]
        spot_params_, spot_map, n_spots = spot_map[1:], spot_map[0], spot_map[3]
        
        #self.print_message('Computing CCF spotted rotating photosphere', index=6, color=34)
        #rotating spotted photosphere ccf 
        ccf_rotating = mstarsim.ccf_rotating_spotted_photosphere(self.obs_time, star_params, 
                                                                  spot_params_, self.cth_coeffs, spot_map, self.dlnp_ph, 
                                                                  self.flpk_ph, self.flnp_ph, self.dlnp_sp, self.flpk_sp, self.flnp_sp, 
                                                                  self.flnp_fc, self.rv_sym, self.ccf_sym_ph, self.ccf_sym_sp, self.ccf_sym_fc, 
                                                                  self.ccf_im, self.flxph, len(self.obs_time), len(self.rv_sym), len(self.flnp_ph), n_spots) 
        
    
        #compute CCF indices
        time_, ccf_indices = self.ccf_indices(self.obs_time, self.rv_sym, ccf_rotating)
        
        #normalize and save the FWHM curve and center the series with respect to its mean 
        ccf_indices = np.array(ccf_indices)
        mean_FWHM = ccf_indices[:,1].astype(float).mean()
        normalized_FWHM_series = []
        for i in range(len(time_)): normalized_FWHM_series.append([time_[i], ccf_indices[i][1] - mean_FWHM])
 
        """
        #center the data with respect to its mean
        self.obs_data = np.array(self.obs_data)
        for i in range(len(self.obs_data)): self.obs_data[i][1] -= self.obs_data[:,1].astype(float).mean()
        """
        
        #compute khi2, if observational data is provided
        if self.obs_data != []:
            normalized_khi2_fwhm = mstarsim.khi2(normalized_FWHM_series, self.obs_data, len(normalized_FWHM_series)) / self.khi2_immaculate_fwhm     
        else:
            normalized_khi2_fwhm = 'not defined'
        
        if write_rv:
            #write RV series on disk
            self.write_file('./output/FWHM_P=' + str(self.Sim.P_rot) + '_' + str(self.Sim.id) + '_Tph=' + str(int(self.Sim.temp_ph)) + '_dTsp=' + str(int(self.Sim.delta_t_sp)) + 
                            '_Q=' + str(self.Sim.Q) + '_diff=' + str(self.Sim.diff_rotation) + '.dat', normalized_FWHM_series, ' ', '')
            
            self.print_message( str(4*' ') + 'File ./output/FWHM_P=' + str(self.Sim.P_rot) + '_Tph=' + str(int(self.Sim.temp_ph)) + '_dTsp=' + str(int(self.Sim.delta_t_sp)) + 
                            '_Q=' + str(self.Sim.Q) + '_diff=' + str(self.Sim.diff_rotation) + '.dat' + " created", 2, 31)            
      
            #compute the residuals
            if self.obs_data != []:
                rv_residuals = []
                for i in range(len(self.obs_time)):
                    rv_residuals.append([self.obs_time[i], normalized_BIS_series[i][1] - self.obs_data[i][1], self.obs_data[i][2] ])
                self.write_file('./output/FWHM_P=' + str(self.Sim.P_rot) + '_Tph=' + str(int(self.Sim.temp_ph)) + '_dTsp=' + str(int(self.Sim.delta_t_sp)) + 
                                '_Q=' + str(self.Sim.Q) + '_diff=' + str(self.Sim.diff_rotation) + '.residuals', rv_residuals, ' ', '')
                
        return normalized_khi2_fwhm

class CCONTRAST(StarsimCommon, Starsim_IO, Spectra, Star):
    """
    """
    def __init__(self, SimulationEnvironment, conf_file_path='./starsim.conf'):
        
        #init configuration file
        self.conf_file_path = conf_file_path
        self.conf_file = self.__conf_init(conf_file_path)
        #class composition (to access the attributes of StarsimSimulation objects)
        self.Sim = SimulationEnvironment
        #object id
        self.id = 'CONTRAST_' + str(int(random.uniform(1,4096)))
        #-------------------------------------------------------------------------
        #sim. parameters 
        #-------------------------------------------------------------------------
        #grid size
        self.ang_res = float(self.conf_file.get('general','ang_res'))
        #initial time of simulation
        self.dt1 = float(self.conf_file.get('general','init_time_sim'))
        #final time of simulation
        self.dt2 = float(self.conf_file.get('general','final_time_sim'))
        #exposure time cadence
        self.t_exp = float(self.conf_file.get('general','time_cadence'))
        #min spectral wavelength
        self.wvmin = float(self.conf_file.get('general','min_spec_range_gen'))
        #max spectral wavelength
        self.wvmax = float(self.conf_file.get('general','max_spec_range_gen'))
        #cth coeffs
        self.cth_coeffs = self.compute_cth_coeffs()
        #type of data (ph/rv/...)
        self.data_type = ''
        #radial velocity CCF range 
        self.ccf_vel_range = float(self.conf_file.get('rv','ccf_vel_range'))
        #spectra resolution
        self.spectra_res = str(self.conf_file.get('general','spectra_resolution'))
        #type of data (ph/rv/...)
        self.data_type = ''
        #-------------------------------------------------------------------------
        #INIT ROUTINES
        #-------------------------------------------------------------------------
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #load time arrays (from .conf file)
        self.load_observational_data(file_path='')
        #init CCFs
        self.init_ccfs()
        #get the normalization factor for the rms
        self.khi2_normalization_factor()
        
        
    def __del__(self):
        """class destructor
        """ 
        
    def __conf_init(self, conf_file_path):
        """creates an instance of class ConfigParser, and read the .conf file. Returns the object created
        ConfigParser, and read the .conf file. Returns the object created
        """
        conf_file_Object = ConfigParser.ConfigParser()
        if not conf_file_Object.read([  conf_file_path]):
            print_message( str(conf_file_path) + ": could not be read", 5, 31, inhibit_msg=False)
            sys.exit()
        else:
            return conf_file_Object
          
    def set_wv_range(self, wvmin, wvmax):
        """setter for changing the wvmin,wvmax and read the new spectras, immac. photosphere and ccfs 
        """
        self.wvmin, self.wvmax = wvmin, wvmax
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #init CCFs
        self.init_ccfs()

    def set_delta_t(self):
        """
        """
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #init CCFs
        self.init_ccfs()

    def systemic_vel(self, rv_sym, ccf_im):
        """fit a gaussian to CCF immaculate photosphere
        """     
        coeff, var_matrix = curve_fit(self.gauss, rv_sym, ccf_im, p0=(0.5,0.0,1.0,0.0))
        self.sys_vel = coeff[1]     
            
        return self.sys_vel
    
    def khi2_normalization_factor(self):
        """computes the khi2 of data in case of 
        """
        #time_, ccf_indices = self.ccf_indices([0.0], self.rv_sym, [self.ccf_sym_ph])
        immaculate_ph_CON = 0.0#ccf_indices[0][1]   
        if self.obs_data != []:
            immaculate_contrast_im = []
            for i in range(len(self.obs_data)):
                time = self.obs_data[i][0]
                immaculate_contrast_im.append([time, immaculate_ph_CON])
            self.khi2_immaculate_contrast = mstarsim.khi2(immaculate_contrast_im, self.obs_data, len(immaculate_contrast_im))
        else:
            self.khi2_immaculate_contrast = None

    def fit_ccf_bisector(self, rv, ccf):
        """returns bis and a set of coefficients that fits 10-90% interval bisector of CCF
        """
        N_BIS_POINTS = 100
        rv_bis = list()
        y_ccf_l = list()
        bisector = []
        #get the coefficients of a gaussian fit gaussian_fit=(a,mean,sigma,offset)
        """
        gaussian_fit = fit_gaussian(rv, ccf)
        """
        #find the gaussian max.
        coeffs, mvar = curve_fit(self.gauss, rv, ccf, p0=(0.5,0.0,1.0,0.0))
        a, mean, sigma, offset = coeffs[:]
        
        max_gaussian_fit = self.gauss(mean, a, mean, sigma, offset)
        #limits of CCF bisector
        #0.1438
        MAX_BIS, MIN_BIS = 0.8*(max_gaussian_fit - offset) , 0.2*(max_gaussian_fit - offset) 
        
        for i in range(N_BIS_POINTS):
            y_ccf = MIN_BIS + i*(MAX_BIS - MIN_BIS) / N_BIS_POINTS
            for j in range(len(ccf)-1):
                if ccf[j] < y_ccf and ccf[j+1] > y_ccf and rv[j] < mean:
                    interpolated_rv_1 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))
                elif ccf[j] > y_ccf and ccf[j+1] < y_ccf and rv[j] > mean:
                    interpolated_rv_2 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))  
            rv_bis.append((interpolated_rv_2 + interpolated_rv_1)/2.0) 
            y_ccf_l.append(y_ccf)
            bisector.append([(interpolated_rv_2 + interpolated_rv_1)/2.0, y_ccf])
        #fit the bisector to a 6-deg polinomial        
        bis_coeffs = np.polyfit(np.array(y_ccf_l), np.array(rv_bis), deg=6, rcond=None, full=False, w=None, cov=False)   
        
        return bisector, bis_coeffs

    def gen_sym_ccf(self, rv, ccf, bis_coeffs):
        """
        """
        rv_sym = rv[:]
        for j in range(len(ccf)):
            rv_sym[j] = rv[j] - self.impoly(ccf[j], bis_coeffs) 
            
        return rv_sym, ccf

    def fit_ccf(self, time_array, rv, ccf_rotating):
        """
        """ 
        #initial guess (A,c,sigma,offset)
        p0 = (1.0,0.0,1.0,1.0)
        coeffs_array = np.zeros(len(time_array), dtype=np.float)
        #fit a gaussian in each ccf
        for i in range(len(time_array)):
            y = ccf_rotating[i][:]
            coeff, var_matrix = curve_fit(self.gauss, rv, y, p0, ftol=1e-05)
            coeffs_array[i] = coeff[1]
            p0 = coeff  #reset the initial guess with the last fit output
       
        return time_array, coeffs_array
     
    def compute_bisector(self, rv,ccf, lims):
        """lims = [min_%, max%]
        """
        N_BIS_POINTS = 50
        rv_bis = list()
        y_ccf_l = list()
        bisector = []
        
        coeffs, mvar = curve_fit(self.gauss, rv, ccf, p0=(0.5,0.0,1.0,0.0))
        a, mean, sigma, offset = coeffs[:]
        
        max_gaussian_fit = self.gauss(mean, a, mean, sigma, offset)
        #limits of CCF bisector
        #0.1438
        MAX_BIS, MIN_BIS = lims[1]*(max_gaussian_fit - offset) , lims[0]*(max_gaussian_fit - offset) 
        
        for i in range(N_BIS_POINTS):
            y_ccf = MIN_BIS + i*(MAX_BIS - MIN_BIS) / N_BIS_POINTS
            for j in range(len(ccf)-1):
                if ccf[j] < y_ccf and ccf[j+1] > y_ccf and rv[j] < mean:
                    interpolated_rv_1 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))
                elif ccf[j] > y_ccf and ccf[j+1] < y_ccf and rv[j] > mean:
                    interpolated_rv_2 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))  
            rv_bis.append((interpolated_rv_2 + interpolated_rv_1)/2.0) 
            y_ccf_l.append(y_ccf)
            bisector.append([(interpolated_rv_2 + interpolated_rv_1)/2.0, y_ccf])
            
        return bisector
    
    def ccf_indices(self, time_array, rv, ccf_rotating):
        """calculate RV, FWHM, BIS, CONTRAST
        """
        indices_array = []
        p0=(0.5,0.0,1.0,0.0)
        
        for i in range(len(time_array)):
            y = ccf_rotating[i][:]
            #fit gaussian
            coeffs, var_matrix = curve_fit(self.gauss, rv, y, p0)
            amplitude, center, sigma, offset = coeffs
            #FWHM = 2.0*sigma*sqrt(2.0*log(2))
            CONTRAST = 100.0 - (amplitude - offset)
            #bisector 
            #bisector_low = np.array( self.compute_bisector(rv,y,lims=[0.10,0.40]) )
            #bisector_high = np.array( self.compute_bisector(rv,y,lims=[0.60,0.90]) )
            #BIS = -(bisector_high[:,0].astype(float).mean() - bisector_low[:,0].astype(float).mean()) 
            BIS = 0.0
            FWHM = 0.0
            #write_file('./bisector_' + str(i), bisector, ' ', '' )
            indices_array.append([center, FWHM, CONTRAST, BIS])     
            p0 = coeffs
        
        return time_array, indices_array
    
    def compute_CONTRAST_series(self, write_rv=True):
        """
        """    
        self.Sim.w = 2.0*pi/self.Sim.P_rot
        star_params = [self.Sim.x_i, self.Sim.R_star, self.Sim.w, self.Sim.B_coeff, self.Sim.C_coeff, 
                       self.Sim.temp_ph, self.Sim.delta_t_sp, self.Sim.delta_t_fc]
        planet_params = [self.Sim.T0_planet, self.Sim.P_planet, self.Sim.A_planet, self.Sim.spin_orbit, self.Sim.bim, self.Sim.R_planet]
        
        star_params = star_params + planet_params
        
        spot_map = self.Sim.spot_map[:]
        spot_params_, spot_map, n_spots = spot_map[1:], spot_map[0], spot_map[3]
        
        #self.print_message('Computing CCF spotted rotating photosphere', index=6, color=34)
        #rotating spotted photosphere ccf 
        ccf_rotating = mstarsim.ccf_rotating_spotted_photosphere(self.obs_time, star_params, 
                                                                  spot_params_, self.cth_coeffs, spot_map, self.dlnp_ph, 
                                                                  self.flpk_ph, self.flnp_ph, self.dlnp_sp, self.flpk_sp, self.flnp_sp, 
                                                                  self.flnp_fc, self.rv_sym, self.ccf_sym_ph, self.ccf_sym_sp, self.ccf_sym_fc, 
                                                                  self.ccf_im, self.flxph, len(self.obs_time), len(self.rv_sym), len(self.flnp_ph), n_spots) 
        
    
        #compute CCF indices
        time_, ccf_indices = self.ccf_indices(self.obs_time, self.rv_sym, ccf_rotating)
        
        #normalize and save the FWHM curve and center the series with respect to its mean 
        ccf_indices = np.array(ccf_indices)
        mean_CON = ccf_indices[:,2].astype(float).mean()
        print "mean CONT:", mean_CON
        normalized_CON_series = []
        for i in range(len(time_)): normalized_CON_series.append([time_[i], ccf_indices[i][2] - mean_CON])
 
        """
        #center the data with respect to its mean
        self.obs_data = np.array(self.obs_data)
        for i in range(len(self.obs_data)): self.obs_data[i][1] -= self.obs_data[:,1].astype(float).mean()
        """
        
        #compute khi2, if observational data is provided
        if self.obs_data != []:
            normalized_khi2_contrast = mstarsim.khi2(normalized_CON_series, self.obs_data, len(normalized_CON_series)) / self.khi2_immaculate_contrast     
        else:
            normalized_khi2_contrast = 'not defined'
        
        if write_rv:
            #write RV series on disk
            self.write_file('./output/CONTRAST_P=' + str(self.Sim.P_rot) + '_' + str(self.Sim.id) + '_Tph=' + str(int(self.Sim.temp_ph)) + '_dTsp=' + str(int(self.Sim.delta_t_sp)) + 
                            '_Q=' + str(self.Sim.Q) + '_diff=' + str(self.Sim.diff_rotation) + '.dat', normalized_CON_series, ' ', '')
            
            self.print_message( str(4*' ') + 'File ./output/CONTRAST_P=' + str(self.Sim.P_rot) + '_Tph=' + str(int(self.Sim.temp_ph)) + '_dTsp=' + str(int(self.Sim.delta_t_sp)) + 
                            '_Q=' + str(self.Sim.Q) + '_diff=' + str(self.Sim.diff_rotation) + '.dat' + " created", 2, 31)            
      
            #compute the residuals
            if self.obs_data != []:
                rv_residuals = []
                for i in range(len(self.obs_time)):
                    rv_residuals.append([self.obs_time[i], normalized_CON_series[i][1] - self.obs_data[i][1], self.obs_data[i][2] ])
                self.write_file('./output/CONTRAST_P=' + str(self.Sim.P_rot) + '_Tph=' + str(int(self.Sim.temp_ph)) + '_dTsp=' + str(int(self.Sim.delta_t_sp)) + 
                                '_Q=' + str(self.Sim.Q) + '_diff=' + str(self.Sim.diff_rotation) + '.residuals', rv_residuals, ' ', '')
                
        return normalized_khi2_contrast

class CBIS(StarsimCommon, Starsim_IO, Spectra, Star):
    """
    """
    def __init__(self, SimulationEnvironment, conf_file_path='./starsim.conf'):
        
        #init configuration file
        self.conf_file_path = conf_file_path
        self.conf_file = self.__conf_init(conf_file_path)
        #class composition (to access the attributes of StarsimSimulation objects)
        self.Sim = SimulationEnvironment
        #object id
        self.id = 'BIS_' + str(int(random.uniform(1,2048)))
        #-------------------------------------------------------------------------
        #sim. parameters 
        #-------------------------------------------------------------------------
        #grid size
        self.ang_res = float(self.conf_file.get('general','ang_res'))
        #initial time of simulation
        self.dt1 = float(self.conf_file.get('general','init_time_sim'))
        #final time of simulation
        self.dt2 = float(self.conf_file.get('general','final_time_sim'))
        #exposure time cadence
        self.t_exp = float(self.conf_file.get('general','time_cadence'))
        #min spectral wavelength
        self.wvmin = float(self.conf_file.get('general','min_spec_range_gen'))
        #max spectral wavelength
        self.wvmax = float(self.conf_file.get('general','max_spec_range_gen'))
        #cth coeffs
        self.cth_coeffs = self.compute_cth_coeffs()
        #type of data (ph/rv/...)
        self.data_type = ''
        #radial velocity CCF range 
        self.ccf_vel_range = float(self.conf_file.get('rv','ccf_vel_range'))
        #spectra resolution
        self.spectra_res = str(self.conf_file.get('general','spectra_resolution'))
        #type of data (ph/rv/...)
        self.data_type = ''
        #-------------------------------------------------------------------------
        #INIT ROUTINES
        #-------------------------------------------------------------------------
        #load time arrays (from .conf file)
        if not self.inhibit_load_spectra:
            #read generated spectra or create if not found
            self.read_precompiled_spectra()
            #immaculate photosphere or read from disk
            self.compute_immaculate_photosphere()
            self.load_observational_data(file_path='')
            #init CCFs
            self.init_ccfs()
            
        #init point of CCF fitting
        self.init_p = (30.0,-0.5,1.0,-5.0)
        #RV jitter
        self.jitter = 0.0
        
        
    def __del__(self):
        """class destructor
        """ 
        
    def __conf_init(self, conf_file_path):
        """creates an instance of class ConfigParser, and read the .conf file. Returns the object created
        ConfigParser, and read the .conf file. Returns the object created
        """
        conf_file_Object = ConfigParser.ConfigParser()
        if not conf_file_Object.read([  conf_file_path]):
            print_message( str(conf_file_path) + ": could not be read", 5, 31, inhibit_msg=False)
            sys.exit()
        else:
            return conf_file_Object
          
    def set_wv_range(self, wvmin, wvmax):
        """setter for changing the wvmin,wvmax and read the new spectras, immac. photosphere and ccfs 
        """
        self.wvmin, self.wvmax = wvmin, wvmax
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #init CCFs
        self.init_ccfs()

    def set_delta_t(self):
        """
        """
        #read generated spectra or create if not found
        self.read_precompiled_spectra()
        #immaculate photosphere or read from disk
        self.compute_immaculate_photosphere()
        #init CCFs
        self.init_ccfs()

    def systemic_vel(self, rv_sym, ccf_im):
        """fit a gaussian to CCF immaculate photosphere
        """     
        coeff, var_matrix = curve_fit(self.gauss, rv_sym, ccf_im, p0=(0.5,0.0,1.0,0.0), ftol=1e-05)
        self.sys_vel = coeff[1]     
            
        return self.sys_vel

    def logL_null_model(self):
        """return logL of data, centered to its mean
        """
        def compute_logL_null(jitter_null):
            """
            """
            inv_sigma2_null = 1.0/(self.err_0**2 + jitter_null**2)
            logL_null = -0.5*(np.sum(((self.y_data_0 - np.mean(self.y_data_0))**2)*inv_sigma2_null + np.log(2.*np.pi) - np.log(inv_sigma2_null))) 
         
            return -logL_null 
        
        init_jitter_0 = 0.01
        res = minimize(compute_logL_null, init_jitter_0, method='Nelder-Mead', tol=1e-6)
       
        self.jitter_null = res.x
        self.logL_null = -compute_logL_null(self.jitter_null)
           
        return self.logL_null, self.jitter_null
           
    def model_logL(self):
        """compute likelihood of model-data
        """
        inv_sigma2 = 1.0/(self.err**2 + self.jitter**2)
        #inv_sigma2 = 1.0/(self.jitter**2)
        self.logL = -0.5*(np.sum(((self.compute_BIS_model() - self.y_data)**2)*inv_sigma2 + np.log(2.*np.pi) - np.log(inv_sigma2))) 
                    
        return self.logL

    def fit_ccf_bisector(self, rv, ccf):
        """returns bis and a set of coefficients that fits 10-90% interval bisector of CCF
        """
        N_BIS_POINTS = 80
        rv_bis = list()
        y_ccf_l = list()
        bisector = []
        #get the coefficients of a gaussian fit gaussian_fit=(a,mean,sigma,offset)
        """
        gaussian_fit = fit_gaussian(rv, ccf)
        """
        #find the gaussian max.
        coeffs, mvar = curve_fit(self.gauss, rv, ccf, p0=(0.5,0.0,1.0,0.0), ftol=1e-05)
        a, mean, sigma, offset = coeffs[:]
        
        max_gaussian_fit = self.gauss(mean, a, mean, sigma, offset)
        #limits of CCF bisector
        #0.1438
        MAX_BIS, MIN_BIS = 0.8*(max_gaussian_fit - offset) , 0.2*(max_gaussian_fit - offset) 
        
        for i in range(N_BIS_POINTS):
            y_ccf = MIN_BIS + i*(MAX_BIS - MIN_BIS) / N_BIS_POINTS
            for j in range(len(ccf)-1):
                if ccf[j] < y_ccf and ccf[j+1] > y_ccf and rv[j] < mean:
                    interpolated_rv_1 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))
                elif ccf[j] > y_ccf and ccf[j+1] < y_ccf and rv[j] > mean:
                    interpolated_rv_2 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))  
            rv_bis.append((interpolated_rv_2 + interpolated_rv_1)/2.0) 
            y_ccf_l.append(y_ccf)
            bisector.append([(interpolated_rv_2 + interpolated_rv_1)/2.0, y_ccf])
        #fit the bisector to a 6-deg polinomial        
        bis_coeffs = np.polyfit(np.array(y_ccf_l), np.array(rv_bis), deg=6, rcond=None, full=False, w=None, cov=False)   
        
        return bisector, bis_coeffs

    def gen_sym_ccf(self, rv, ccf, bis_coeffs):
        """
        """
        rv_sym = rv[:]
        for j in range(len(ccf)):
            rv_sym[j] = rv[j] - self.impoly(ccf[j], bis_coeffs) 
            
        return rv_sym, ccf

    def fit_ccf(self, time_array, rv, ccf_rotating):
        """
        """ 
        #initial guess (A,c,sigma,offset)
        #p0 = (1.0,0.0,1.0,1.0)
        coeffs_array = np.zeros(len(time_array), dtype=np.float)
        #fit a gaussian in each ccf
        for i in range(len(time_array)):
            y = ccf_rotating[i][:]
            coeff, var_matrix = curve_fit(self.gauss, rv, y, p0=self.init_p, ftol=1e-05)
            coeffs_array[i] = coeff[1]
            self.init_p = coeff  #reset the initial guess with the last fit output
       
        return time_array, coeffs_array
      
    def compute_bisector(self, rv,ccf, lims):
        """lims = [min_%, max%]
        """
        N_BIS_POINTS = 80
        rv_bis = list()
        y_ccf_l = list()
        bisector = []
        
        coeffs, mvar = curve_fit(self.gauss, rv, ccf, p0=self.init_p, ftol=1e-05)
        self.init_p = coeffs
        
        a, mean, sigma, offset = coeffs[:]
        
        max_gaussian_fit = self.gauss(mean, a, mean, sigma, offset)
        #limits of CCF bisector
        #0.1438
        MAX_BIS, MIN_BIS = lims[1]*(max_gaussian_fit - offset) , lims[0]*(max_gaussian_fit - offset) 
     
        for i in range(N_BIS_POINTS):
            y_ccf = MIN_BIS + i*(MAX_BIS - MIN_BIS) / N_BIS_POINTS
            for j in range(len(ccf)-1):
                if ccf[j] < y_ccf and ccf[j+1] > y_ccf and rv[j] < mean:
                    interpolated_rv_1 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))
                elif ccf[j] > y_ccf and ccf[j+1] < y_ccf and rv[j] > mean:
                    interpolated_rv_2 = rv[j] + (((y_ccf - ccf[j])*(rv[j+1]-rv[j])) / (ccf[j+1]-ccf[j]))  
            rv_bis.append((interpolated_rv_2 + interpolated_rv_1)/2.0) 
            y_ccf_l.append(y_ccf)
            bisector.append([(interpolated_rv_2 + interpolated_rv_1)/2.0, y_ccf])
        
        return bisector
    
    def ccf_indices(self, time_array, rv, ccf_rotating):
        """calculate RV, FWHM, BIS, CONTRAST
        """
        
        indices_array = []
        
        for i in range(len(time_array)):
            y = ccf_rotating[i][:]
            #fit gaussian
            coeffs, var_matrix = curve_fit(self.gauss, rv, y, p0=self.init_p)
            self.init_p = coeffs
            amplitude, center, sigma, offset = coeffs
            FWHM = 2.0*sigma*sqrt(2.0*log(2))
            CONTRAST = 0.0
            #bisector 
            bisector_low = np.array( self.compute_bisector(rv,y,lims=[0.10,0.40]) )
            bisector_high = np.array( self.compute_bisector(rv,y,lims=[0.60,0.90]) )
            BIS = -(bisector_high[:,0].astype(float).mean() - bisector_low[:,0].astype(float).mean()) 
            #write_file('./bisector_' + str(i), bisector, ' ', '' )
            indices_array.append([center, FWHM, CONTRAST, BIS])     
        
        return time_array, indices_array
    
    def compute_BIS_series(self, write_rv=True):
        """
        """        
        self.Sim.w = 2.0*pi/self.Sim.P_rot
        star_params = [self.Sim.x_i, self.Sim.R_star, self.Sim.w, self.Sim.B_coeff, self.Sim.C_coeff, 
                       self.Sim.temp_ph, self.Sim.delta_t_sp, self.Sim.delta_t_fc]
        planet_params = [self.Sim.T0_planet, self.Sim.P_planet, self.Sim.A_planet, self.Sim.spin_orbit, self.Sim.bim, self.Sim.R_planet]
        
        star_params = star_params + planet_params
        
        spot_map = self.Sim.spot_map[:]
        spot_params_, spot_map, n_spots = spot_map[1:], spot_map[0], spot_map[3]
        
        #self.print_message('Computing CCF spotted rotating photosphere', index=6, color=34)
        #rotating spotted photosphere ccf 
        ccf_rotating = mstarsim.ccf_rotating_spotted_photosphere(self.obs_time, star_params, 
                                                                  spot_params_, self.cth_coeffs, spot_map, self.dlnp_ph, 
                                                                  self.flpk_ph, self.flnp_ph, self.dlnp_sp, self.flpk_sp, self.flnp_sp, 
                                                                  self.flnp_fc, self.rv_sym, self.ccf_sym_ph, self.ccf_sym_sp, self.ccf_sym_fc, 
                                                                  self.ccf_im, self.flxph, len(self.obs_time), len(self.rv_sym), len(self.flnp_ph), n_spots) 
        
        """
        #fit a gaussian and get its center (rv)
        time_, rv_series = self.fit_ccf(self.obs_time, self.rv_sym, ccf_rotating)
        """
        
        #compute CCF indices
        time_, ccf_indices = self.ccf_indices(self.obs_time, self.rv_sym, ccf_rotating)
        
        
        #normalize and save the BIS curve and center the series with respect to its mean 
        ccf_indices = np.array(ccf_indices)
        mean_BIS = ccf_indices[:,3].astype(float).mean()
        normalized_BIS_series = []
        for i in range(len(time_)): normalized_BIS_series.append([time_[i], ccf_indices[i][3] - mean_BIS])
        
        """
        #center the data with respect to its mean
        self.obs_data = np.array(self.obs_data)
        for i in range(len(self.obs_data)): self.obs_data[i][1] -= self.obs_data[:,1].astype(float).mean()
        """
        
        #compute khi2, if observational data is provided
        if self.obs_data != []:
            normalized_khi2_bis = mstarsim.khi2(normalized_BIS_series, self.obs_data, len(normalized_BIS_series)) / self.khi2_immaculate_bis     
        else:
            normalized_khi2_bis = 'not defined'
        
        if write_rv:
            #write RV series on disk
            self.write_file('./output/BIS_P=' + str(self.Sim.P_rot) + '_' + str(self.Sim.id) + '_Tph=' + str(int(self.Sim.temp_ph)) + '_dTsp=' + str(int(self.Sim.delta_t_sp)) + 
                            '_Q=' + str(self.Sim.Q) + '_diff=' + str(self.Sim.diff_rotation) + '.dat', normalized_BIS_series, ' ', '')
            
            self.print_message( str(4*' ') + 'File ./output/BIS_P=' + str(self.Sim.P_rot) + '_Tph=' + str(int(self.Sim.temp_ph)) + '_dTsp=' + str(int(self.Sim.delta_t_sp)) + 
                            '_Q=' + str(self.Sim.Q) + '_diff=' + str(self.Sim.diff_rotation) + '.dat' + " created", 2, 31)            
      
            #compute the residuals
            if self.obs_data != []:
                rv_residuals = []
                for i in range(len(self.obs_time)):
                    rv_residuals.append([self.obs_time[i], normalized_BIS_series[i][1] - self.obs_data[i][1], self.obs_data[i][2] ])
                self.write_file('./output/BIS_P=' + str(self.Sim.P_rot) + '_Tph=' + str(int(self.Sim.temp_ph)) + '_dTsp=' + str(int(self.Sim.delta_t_sp)) + 
                                '_Q=' + str(self.Sim.Q) + '_diff=' + str(self.Sim.diff_rotation) + '.residuals', rv_residuals, ' ', '')
                
        return normalized_khi2_bis
                            
class SimulatedAnnealing:
    """environment for optimization 
    """
    def __init__(self, conf_file_path='./starsim.conf'):
        #init configuration file
        self.conf_file_path = conf_file_path
        self.conf_file = self.__conf_init(conf_file_path)
        #----------------------------------------------------------------------
        #sim. parameters 
        #----------------------------------------------------------------------
        #MCSA number of iterations per T step
        self.MC_iterations = int(self.conf_file.get('MCSA','MC_iterations'))
        self.alpha = float(self.conf_file.get('MCSA','MC_alpha'))
        self.T0 = float(self.conf_file.get('MCSA','MC_T0'))
        #---------------------------------------------------------------------- 
        self.RV_WEIGHT = 4.0
        self.InvertedSpotMap = []
        self.s_jstat = np.inf
        #self.jitter = 0.0000
        
              
    def __del__(self):
        """class destructor
        """  

    def __conf_init(self, conf_file_path):
        """creates an instance of class ConfigParser, and read the .conf file. Returns the object created
        ConfigParser, and read the .conf file. Returns the object created
        """
        conf_file_Object = ConfigParser.ConfigParser()
        if not conf_file_Object.read([  conf_file_path]):
            print_message( str(conf_file_path) + ": could not be read", 5, 31, inhibit_msg=False)
            sys.exit()
        else:
            return conf_file_Object
        
    def kirkpatrick_cooling(self, T, alpha):
        while True:
            yield T
            T = alpha*T
    
    def compute_joint_statistic(self, DataObject):
        """
        """
        joint_logL = 0.0
       
        try:
            for object in DataObject:
                if object.data_type == 'ph':
                    joint_logL += object.model_logL()
                      
                elif object.data_type == 'rv':
                    joint_logL += self.RV_WEIGHT*object.model_logL() 
                    
        except ValueError:    
            raise
             
        return -joint_logL 
        
    def compute_metropolis(self, DataObjects):
        """
        """
        #self.old_jitter = self.jitter
        if self.rnd_walk:
            self.MAX_ITERATIONS = 10*self.MAX_ITERATIONS
         
        iteration = 0
        while iteration < self.MAX_ITERATIONS:
            s_state = self.spot_map
            p_state = s_state
            #perturb the state
            #if random.uniform(0,1) < 0.65:
            p_state = self.gen_modified_state(s_state, self.rnd_walk)
            #else: 
            #    self.old_jitter = self.jitter
            #    self.jitter = random.uniform(0.0,0.01)    
            alarm = 0
            #check if this perturbation overlaps a pair of spots. If so,repeat generation
            while 1:
                if self.check_params_compatibility(p_state) == False:
                    p_state = self.gen_modified_state(s_state, self.rnd_walk)    
                    alarm += 1
                    if alarm > 10:
                        #try to make a tiny perturbation
                        p_state = self.gen_modified_state(s_state, True)
                        #if not, recover the last config
                        if self.check_params_compatibility(p_state) == False:
                            p_state = s_state
                            alarm = 0
                            break
                #if config compatible, break            
                else:
                    alarm = 0
                    break
           
            #test the perturbed configuration    
            self.spot_map = p_state       
            self.p_jstat = self.compute_joint_statistic(DataObjects)
            
            #imporoves?
            delta_jstat = self.p_jstat - self.s_jstat
            
            if delta_jstat < 0.0:
                #adopt the new configuration
                self.s_jstat = self.p_jstat
                
                #print stat.
                if not Glob.inhibit_msg:
                    sys.stdout.write(str() + "\t Process id: " + str(self.id) + "  Temperature: " + str(1.0/self.beta) +  \
                                    "  |=====> Accepted joint logL: %f \r" % (-self.s_jstat))
                    sys.stdout.flush()
                
            #otherwise do a random decision
            elif random.uniform(0,1) < exp(-self.beta*delta_jstat):
                #print "Entered with T=", 1.0/self.beta
                self.s_jstat = self.p_jstat
     
            else:
                #recover old spot configuration --> no change in spotmap
                self.spot_map = s_state[:]
                #self.jitter = self.old_jitter
                
            iteration += 1
        
        return self.spot_map    
       
    def simulated_annealing(self, DataObjects):
        """
        """
        #os.system("clear")
        self.init_temp = 0
        self.temp_steps = 25
        #set a temperature cooling schedule. (initial T, alpha)
        temperatures = self.kirkpatrick_cooling(self.T0, self.alpha)
        #find init time and final time of simulation (dt1 and dt2)
        #take the widest range possible
        self.dt1 = min([min(object.obs_time) for object in DataObjects])  
        self.dt2 = max([max(object.obs_time) for object in DataObjects])  
        
        for i in range(self.init_temp): self.beta = 1.0/temperatures.next()
        #init temperature sweep
        for self.p_temp in range(self.init_temp, self.temp_steps):
            
            self.MAX_ITERATIONS = int(self.MC_iterations) #*sqrt(self.p_temp+1))
            self.beta = 1.0/temperatures.next()
            #print "Temp. step:", p_temp, "Beta:", self.beta, "Iterations:", self.MAX_ITERATIONS              
            if self.p_temp == self.temp_steps - 1:
                self.is_last_beta = True
                self.rnd_walk = True
            else:
                self.rnd_walk = False
                self.is_last_beta = False
            
            state_out = self.compute_metropolis(DataObjects)
            
        self.Ejoint_statistic = self.compute_joint_statistic(DataObjects)
                
        if not Glob.inhibit_msg:
            self.write_file(str(DataObjects[0].conf_file.get('files','output_path')) + 
                            '/spot_map_' + str(self.id) + '.dat', state_out[0], ' ', '')
            
        self.InvertedSpotMap = state_out[0]
        
        return -self.Ejoint_statistic, self.id, self.info, self.RV_WEIGHT 
               
class StarsimSimulation(Star, Spots, Planet, SimulatedAnnealing):
    """packs all information related to simulation config
    """
    def __init__(self, conf_file_path='./starsim.conf'):
        #init configuration file
        self.conf_file_path = conf_file_path
        self.conf_file = self.__conf_init(conf_file_path)
        #instantiate the classes inherited to share their class variables
        Star.__init__(self, self.conf_file)
        Planet.__init__(self, self.conf_file)
        Spots.__init__(self, self.conf_file)
        SimulatedAnnealing.__init__(self)
        #object id
        self.id = int(random.uniform(1,1e7))
        #simulation info
        self.info = ''
        #spectra resolution
        self.spectra_res = str(self.conf_file.get('general','spectra_resolution'))
        #init spotmap
        self.init_spotmap()
        self.dt1 = float(self.conf_file.get('general','init_time_sim'))
        self.dt2 = float(self.conf_file.get('general','final_time_sim'))
        self.t_exp = float(self.conf_file.get('general','time_cadence'))
        self.n_spots = int(self.spot_map[3])
        self.IterCounter = 0
    
    def __del__(self):
        """class destructor
        """  
       
    def init_spotmap(self):
        """read the spotmap in file
        """  
        self.spot_map = self.read_spot_params(spot_list_file='./data/spot_list.dat')
        
    def print_start_info(self):
        """
        """
        self.print_message('', 7,94)
        self.print_message('                     STARSIM 2.0                       ', 7, 90)
        #self.print_message('                                                       ', 7, 94)
        self.print_message('     SIMULATION OF STELLAR PHOTOMETRIC/RV TIME SERIES  ', 7, 94)                                     
        #self.print_message('                                                       ', 7, 94)
        self.print_message('                                                       ', 7, 90)
        self.print_message('', 7,90)
        self.print_message("Read configuration file " + str(self.conf_file_path), index=2, color=31)
        #self.print_message( 'Filter TF data file: ' + str(conf_file.get('files','filter_path')), 3, 37)
        #self.print_message( '', 7,37)
        #self.print_message( 'Number of pixels on the star: ' + str(n1), 3, 37)
        #self.print_message( 'Size of the pixel: ' + str(da1) + 'deg', 3, 37)
        #self.print_message( 'Radius of the star: ' + str(star.R_star) + 'Rsun', 3, 37)
        if Glob.ncpus != 0:
            self.print_message( 'Detected CPUs / using CPUs: ' + str(mp.cpu_count()) + "/" + str(Glob.ncpus), 5, 95)
            self.print_message('', 5,95)
            self.print_message('INIT. STARSIM...', 6,90)
    
    def __conf_init(self, conf_file_path):
        """creates an instance of class ConfigParser, and read the .conf file. Returns the object created
        ConfigParser, and read the .conf file. Returns the object created
        """
        conf_file_Object = ConfigParser.ConfigParser()
        if not conf_file_Object.read([conf_file_path]):
            print_message( str(conf_file_path) + ": could not be read", 5, 31, inhibit_msg=False)
            sys.exit()
        else:
            return conf_file_Object

    def update_spectral_data(self, delta_t_sp, delta_t_fc, Observables=None):
        """set the new spot/faculae temperatures to the observables objects 
        """

        self.delta_t_sp, self.delta_t_fc = delta_t_sp, delta_t_fc
        if Observables == None:
            Observables = self.Observables
        if type(Observables) is list:
            for ObservableObject in Observables:
                ObservableObject.set_delta_t() 
        else:
            Observables.set_delta_t()

    def set_evo_rate(self, evo_rate):
        """setter for star inclination
        """
        self.evo_rate = evo_rate
        self.spot_map[1] = evo_rate

    def set_Q(self, Q):
        """setter for Q factor
        """
        self.Q = Q
        self.spot_map[4] = Q
        
    def set_diff_rotation(self, df_rot):
        """setter for differential rotation
        """
        self.df_rot = df_rot
        self.spot_map[2] = df_rot
    
    def get_diff_rotation(self):
        """getter for differential rotation
        """
        return self.spot_map[2]
        
    def read_multifit_config(self, multifit_file_path):
        """build an array of Photometry, RV objects 
        """
        self.Observables = []  #the list of objects (can be Ph, RV or indices)
        #read conf. file
        try:
            mfit_data = self.read_file(multifit_file_path, ' ')
        except IOError:
            print "Multifit .conf file could not be read!"
        #read the multifit configuration file
        for j in range(len(mfit_data)):
            #labeled as '1' and is photometry
            if str(mfit_data[j][2]) == '1' and str(mfit_data[j][1]) == 'ph':
                #generate a photometry object
                ph_object = CPhotometry(self, inhibit_load_spectra=True)
                #set wavelength range
                wvmin, wvmax = float(mfit_data[j][3]), float(mfit_data[j][4])
                ph_object.set_wv_range(wvmin, wvmax)
                #time / data arrays
                ph_object.science_data = file_path=mfit_data[j][0]
                ph_object.load_observational_data(file_path=mfit_data[j][0])
                #type label
                ph_object.data_type = mfit_data[j][1]
                #info field (data path)
                ph_object.info = str(mfit_data[j][0]) 
                #filter name
                ph_object.filter_name = str( mfit_data[j][5]).split('/')[-1].split('.dat')[0]
                #photometric offset
                ph_object.z_ph = float(mfit_data[j][6])
                #photometric jitter
                try:
                    ph_object.jitter = float(mfit_data[j][7])
                except:
                    ph_object.jitter = 0.0 
                #null model logL
                ph_object.logL_null_model()
                #finally, append this object to the list
                self.Observables.append( ph_object )
            
            #labeled as '1' and is rv   
            elif str(mfit_data[j][2]) == '1' and str(mfit_data[j][1]) == 'rv':
                #generate a RV object
                rv_object = CRadialVelocity(self, inhibit_load_spectra=True)
                #set wavelength range
                wvmin, wvmax = float(mfit_data[j][3]), float(mfit_data[j][4])
                rv_object.set_wv_range(wvmin, wvmax)
                #time / data arrays
                rv_object.science_data = file_path=mfit_data[j][0]
                rv_object.load_observational_data(file_path=mfit_data[j][0])
                #type label
                rv_object.data_type = mfit_data[j][1]
                #RV jitter
                try:
                    rv_object.jitter = float(mfit_data[j][7])
                except:
                    rv_object.jitter = 0.0
                #null model logL
                rv_object.logL_null_model()
                #finally, append it to the list
                self.Observables.append( rv_object )
                
            #labeled as '1' and is BIS   
            elif str(mfit_data[j][2]) == '1' and str(mfit_data[j][1]) == 'bis':
                #generate a RV object
                BIS_object = CBIS(self)
                #set wavelength range
                wvmin, wvmax = float(mfit_data[j][3]), float(mfit_data[j][4])
                BIS_object.set_wv_range(wvmin, wvmax)
                #time / data arrays
                BIS_object.science_data = file_path=mfit_data[j][0]
                BIS_object.load_observational_data(file_path=mfit_data[j][0])
                BIS_object.khi2_normalization_factor()
                #type label
                BIS_object.data_type = mfit_data[j][1]
                #finally, append it to the list
                self.Observables.append( BIS_object )
                
            #labeled as '1' and is FWHM   
            elif str(mfit_data[j][2]) == '1' and str(mfit_data[j][1]) == 'fwhm':
                #generate a RV object
                FWHM_object = CFWHM(self)
                #set wavelength range
                wvmin, wvmax = float(mfit_data[j][3]), float(mfit_data[j][4])
                FWHM_object.set_wv_range(wvmin, wvmax)
                #time / data arrays
                FWHM_object.science_data = file_path=mfit_data[j][0]
                FWHM_object.load_observational_data(file_path=mfit_data[j][0])
                FWHM_object.khi2_normalization_factor()
                #type label
                FWHM_object.data_type = mfit_data[j][1]
                #finally, append it to the list
                self.Observables.append( FWHM_object )
                
            #labeled as '1' and is contrast   
            elif str(mfit_data[j][2]) == '1' and str(mfit_data[j][1]) == 'contrast':
                #generate a RV object
                CONTRAST_object = CCONTRAST(self)
                #set wavelength range
                wvmin, wvmax = float(mfit_data[j][3]), float(mfit_data[j][4])
                CONTRAST_object.set_wv_range(wvmin, wvmax)
                #time / data arrays
                CONTRAST_object.science_data = file_path=mfit_data[j][0]
                CONTRAST_object.load_observational_data(file_path=mfit_data[j][0])
                CONTRAST_object.khi2_normalization_factor()
                #type label
                CONTRAST_object.data_type = mfit_data[j][1]
                #finally, append it to the list
                self.Observables.append( CONTRAST_object )
                
        #find min time and max time for current datasets
        self.dt1 = min([min(object.obs_time) for object in self.Observables])  
        self.dt2 = max([max(object.obs_time) for object in self.Observables]) 
         
        return self.Observables
     

          
           
