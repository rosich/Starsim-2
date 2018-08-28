#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys, os
from math import sqrt, sin, cos, pi, tan, atan
import random


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
            
            
if __name__ == '__main__':
    """
    """
    #read data from data file
    Data = read_file(sys.argv[1], ' ')
    
    data_noisy = []
    for i in range(len(Data)):
        #sigma = abs(random.uniform(0.00005,0.0002)) #ph
        sigma = abs(random.uniform(1.2, 2.5))      #rv
        data_noisy.append([float(Data[i][0]), \
                         random.gauss(float(Data[i][1]),sigma), sigma ])
        
    write_file(str(sys.argv[1]) + '.noisy', data_noisy, ' ', '')
    
    