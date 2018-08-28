import numpy as np
from scipy.optimize import curve_fit
from math import pi


def gaussian(x_val, a, mean, sigma, offset):
   
    return a*np.exp(-(x_val-mean)**2/(2.0*sigma**2)) + offset


def fit_gaussian(x,y):
    
    fit_coeffs, covar = curve_fit(gaussian, x, y)
    
    return fit_coeffs


def impoly(x, in_coeffs):
    """return the image of a n-deg polynomial.
       input: x and array of n+1 coefficients, in decreasing order.
    """ 
    im = 0.0
    for order in range(0,len(in_coeffs)):
        im += in_coeffs[order]*(x**(len(in_coeffs)-order-1))
   
    return im
        