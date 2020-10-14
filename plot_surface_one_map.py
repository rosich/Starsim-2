#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import numpy as np 
import random
from math import sin, cos, acos, log10
import sys

def to_deg(rad):
    """radians --> degrees
    """
    return rad*180.0/np.pi

def to_radians(deg):
    """degrees --> radians
    """
    return deg*np.pi/180.0

def read_file(in_file, separator):
    """
    """
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
                    g_list.append(row)
                    row = []
    return g_list

def write_file(out_file, data, separator, header):
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
            
def spherical_distance(th_0, ph_0, th, ph):
    """
    """
    try:
        to_radians = np.pi/180.0
        to_deg = 1.0/to_radians
        cosD = sin(np.pi/2.0 - th_0*to_radians)*sin(np.pi/2.0 - th*to_radians) + \
               cos(np.pi/2.0 - th_0*to_radians)*cos(np.pi/2.0 - th*to_radians)*cos(abs(ph_0 - ph)*to_radians)
        distance = to_deg*acos(cosD)
        return distance
    
    except ValueError:
        
        return 0.0

def plot(Map, time):
    """
    """
    alpha = []
    delta = []
    probability = []
    for i in range(len(Map)):
        alpha.append(float(Map[i][0]))
        delta.append(float(Map[i][1]))
        probability.append(float(Map[i][2]))
    
    #colors = [0.75 for i in range(n)]
    cm = plt.cm.get_cmap('gist_yarg')
    #make figure
    fig = plt.figure(1,figsize=(8,6))
    plt.subplots_adjust(left=0.05, right=0.85, top=0.98, bottom=0.05, hspace=0.001)
    ax = fig.add_subplot(111, projection='lambert' ) # or mollweide, aitoff, or lambert. [example](http://matplotlib.org/examples/pylab_examples/geo_demo.html) 
    ax.grid(True)
    spotmap = ax.scatter(alpha,delta, marker='s', c=probability, s=0.35, edgecolors='none', cmap=cm, norm=colors.PowerNorm(gamma=1./10.5))
    ax.set_xticklabels(['','','','','','','','','','','']) #use if you want to change to a better version of RA. 
    #ax.set_xticklabels(['','240$^\circ$','270$^\circ$','300$^\circ$','330$^\circ$','0$^\circ$','30$^\circ$','60$^\circ$','90$^\circ$','120$^\circ$',''], fontsize=10)
    #ax.set_xticklabels(['','0$^\circ$','30$^\circ$','60$^\circ$','90$^\circ$','120$^\circ$','150$^\circ$','180$^\circ$','210$^\circ$','240$^\circ$',''], fontsize=10)
    #plt.xlabel(r'$\lambda$')
    #plt.ylabel(r'$\delta$')
    #cb = plt.colorbar(spotmap,ax=ax, fraction=0.035, pad=0.05, aspect=20)
    plt.show()
    #save .png
    #fig.savefig('s_' + str(time) + '.png', dpi=600)

def reconstruct_surface(time):
    """
    """
    P_rot = 15.0
    phi_t = (360.0/P_rot)*time
    #scan the surface
    N_div = 1500
    Surface = []
    for i in range(2*N_div):
        alpha = -np.pi + 2*np.pi*float(i)/(2.*N_div)
        print alpha
        for j in range(N_div):
            delta = -np.pi/2 + np.pi*(float(j)/N_div)
            #sum over all spot maps
            probability = 0.0
            for SpotMap in SpotMaps:
                #SpotMap = SpotMap[0]
                for spot in range(len(SpotMap)):
                    t_ini, t_life = float(SpotMap[spot][0]), float(SpotMap[spot][1])
                    radii = float(SpotMap[spot][4])
                    th_0, ph_0 = float(SpotMap[spot][2]), float(SpotMap[spot][3]) + phi_t#120.0
                    #check time condition
                    if (t_ini < time) and (t_ini + t_life > time): 
                        deg_alpha = to_deg(alpha)
                        if deg_alpha >= 0.0 and deg_alpha <= 180.0:
                            ph = deg_alpha
                        else:
                            ph = 360.0 + deg_alpha
                        th = 90.0 - to_deg(delta)
                        #check distance condition
                        if spherical_distance(th_0, ph_0, th, ph) < radii:
                            probability += 1.0/len(SpotMaps)
            
            #append pixel to the global spot map            
            Surface.append([alpha, delta, probability])

    return Surface

def plot_active_longitudes(t, N_div):
    """
    """
    #scan the surface
    Surface = []
    for i in range(2*N_div):
        #longitude
        alpha = -np.pi + 2*np.pi*float(i)/(2.*N_div)
        #sum over all spot maps
        probability = 0.0
        for SpotMap in SpotMaps:
            for spot in range(len(SpotMap)):
                t_ini, t_life = float(SpotMap[spot][0]), float(SpotMap[spot][1])
                radii = float(SpotMap[spot][4])
                th_0, ph_0 = float(SpotMap[spot][2]), float(SpotMap[spot][3])
                #check time condition
                if (t_ini < t) and (t_ini + t_life > t): 
                    deg_alpha = to_deg(alpha)
                    if deg_alpha >= 0.0 and deg_alpha <= 180.0:
                        ph = deg_alpha
                    else:
                        ph = 360.0 + deg_alpha
                    th = 90.0 #90.0 - to_deg(delta)
                    th_0 = 90.0
                    #check distance condition
                    if spherical_distance(th_0, ph_0, th, ph) < radii:
                        probability += 1.0/len(SpotMaps)
        
        #append pixel to the global spot map            
        Surface.append([alpha, t, probability])

    return Surface
    
       
# get random values to plot
#needs to be in radians from (-pi,pi) & (-pi/2, pi/2)
if __name__ == '__main__':
    """
    """
    only_plot = False

    #read spotmap files
    import glob, shutil

    # a single map
    map_files = glob.glob(sys.argv[1])    
    qfile = read_file(sys.argv[1], ' ')

    SpotMaps = []
    for map_file in map_files:
        SpotMap = read_file(map_file, ' ')
        SpotMaps.append(SpotMap)
    
    #use the first map to compute the mean time
    SpotMap = SpotMaps[0]
    t_init, t_life, colat, lon, r = zip(*SpotMap)
    t_init = map(float,t_init)
    
    if only_plot:
        t_init, t_max = min(t_init), max(t_init)
        t_mean = (t_max + t_init)/2.0
        SpotMap = read_file('ReconstructedSurface_PhRV.dat', ' ')
        plot(SpotMap, t_max)
        sys.exit()
    
    
    t_min, t_max = min(t_init), max(t_init)
    t_mean = (t_max + t_min)/2.0
    
    time = 25.36
    
    print "Time:", time
    #t_init, t_max = 0.0, 50.0
    Surface = reconstruct_surface(time)
    write_file('ReconstructedSurface_PhRV.dat', Surface, ' ', '')
    #plot it
    plot(Surface, time)
        
      
                
                
                
                
                





