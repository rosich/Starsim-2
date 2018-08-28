import sys
from math import sqrt, log
import numpy as np
sys.path.append('./src')
import starsim_generic
from scipy.optimize import curve_fit

def gauss(x, *p):
    """gauss function def
    """
    A, c, sigma, offset = p
    
    return A*np.exp(-((x-c)**2)/(2.0*sigma**2)) + offset
    
def fit_gaussian(x,y):
    """
    """
    fit_coeffs, covar = curve_fit(gauss, x, y)
    
    return fit_coeffs

def fit_ccf_bisector(rv,ccf):
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
    coeffs, mvar = curve_fit(gauss, rv, ccf, p0=(0.5,0.0,1.0,0.0))
    a, mean, sigma, offset = coeffs[:]
    
    max_gaussian_fit = gauss(mean, a, mean, sigma, offset)
    
    
    #limits of CCF bisector
    #0.1438
    MAX_BIS, MIN_BIS = 0.8*(max_gaussian_fit - offset) , 0.15*(max_gaussian_fit - offset) 
    
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

def gen_sym_ccf(rv, ccf, bis_coeffs):
    """
    """
    rv_sym = rv[:]
    for j in range(len(ccf)):
        rv_sym[j] = rv[j] - starsim_generic.impoly(ccf[j], bis_coeffs) 
        
    return rv_sym, ccf

def fit_ccf(time_array, rv, ccf_rotating):
    """
    """    
    coeffs_array = []
    for i in range(len(time_array)):
        y = ccf_rotating[i][:]
        coeff, var_matrix = curve_fit(gauss, rv, y, p0=(0.5,0.0,1.0,0.0))
        coeffs_array.append(coeff[1])     
        
    return time_array, coeffs_array

def ccf_indicators(time_array, rv, ccf_rotating):
    """calculate RV, FWHM, BIS, CONTRAST
    """
    indicators_array = []
    for i in range(len(time_array)):
        y = ccf_rotating[i][:]
        #fit gaussian
        coeff, var_matrix = curve_fit( gauss, rv, y, p0=(0.5,0.0,1.0,0.0) )
        amplitude, center, sigma, offset = coeff
        FWHM = 2.0*sigma*sqrt(2.0*log(2))
        CONTRAST = amplitude - offset
        BIS = 'Not_yet_implemented' 
        """
        
        
        
        
        
        """
        indicators_array.append([center, FWHM, CONTRAST, BIS])     
        
    return time_array, indicators_array 
    
def systemic_vel(rv_sym, ccf_im):
    """fit a gaussian to CCF immaculate photosphere
    """     
    coeff, var_matrix = curve_fit(gauss, rv_sym, ccf_im, p0=(0.5,0.0,1.0,0.0))
    sys_vel = coeff[1]     
        
    return sys_vel
    
    
    
    
    
    
    
    
    
    
    
    
    
    