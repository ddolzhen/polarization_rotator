# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 09:14:27 2021

@author: liamm

Reads in the out_data_jun14.csv and reformats the file and makes a histogram of 
DOP and Power

INPUT FORMAT:
    
    <ROTATOR>,<PADDLE1>,<PADDLE2>,<PADDLE3>,<S1>,<S2>,<S3>,<DOP>,<POWER>
    -each 4 consecutive lines a repeats for same rotator and paddle settings.
    
OUTPUT FORMAT:
    <ROTATOR>,<PADDLE1>,<PADDLE2>,<PADDLE3>,<S1>,<S2>,<S3>,<DOP>,<POWER>
    -each trial is now averaged and combined into 1 line
    -lines come paired by paddle settings but different rotator settings

"""

import sys
import numpy as np
import matplotlib.pyplot as pyplot




def main():
    
    if len(sys.argv) != 3:
        print("the program requires 3 arguments")
        print("format: <program_file> <input_file> <output_file>")
        quit()
    else: 
        infile_name = sys.argv[1]
        outfile_name = sys.argv[2]
        
        
    infile = open(infile_name, "r")
    
    content = infile.readlines()
    
    del content[0]

    
    for i in range(len(content)):
        content[i] = content[i].split(",")
        

    out_data = np.zeros(shape = (len(content[1]), int(len(content)/4 )) )
    
    
    for j in range(len(content[0])):
        for k in range(int(len(content)/4)):
            out_data[j][k] = (float(content[4*k][j]) + float(content[4*k+1][j]) + float(content[4*k+2][j]) + float(content[4*k+3][j]))/4    
    
    """
    #out_data has format:
        - each row is one of the variables (rotator, s1, s2, etc...)
        -each column is corresponding value for that run.
    
    """
    #define the resulting data arrays
    
    #first order the data based on paddle_1 settings
    sorting = np.lexsort((out_data[3], out_data[2], out_data[1]))
    
    reformatted = np.ndarray(shape = (len(out_data), len(out_data[0])))
    for l in range(len(out_data)):
        for m in range(len(sorting)):
            reformatted[l][m] = out_data[l][sorting[m]]
    #print(reformatted)
    #print(reformatted[0])
    
    
    outfile = open(outfile_name, "w")
    outfile.write("#Rotator,#Paddle1,#Paddle2,#Paddle3,#S1,#S2,#S3,#DOP,#Power\n")
    for a in range(len(reformatted[0])):
        outfile.write("{0:f},{1:f},{2:f},{3:f},{4:f},{5:f},{6:f},{7:f},{8:f}\n".format(reformatted[0][a],reformatted[1][a], reformatted[2][a], reformatted[3][a],
                                                                                     reformatted[4][a],reformatted[5][a],reformatted[6][a],reformatted[7][a],
                                                                                     reformatted[8][a]))
    outfile.close()

    #Now create histograms of DOP and Power.
    #1st DOP histogram
    pyplot.hist(reformatted[7], bins = 100, range= (0.7, 1.2))
    pyplot.title("Histogram of Degree of Polarization values")
    
    pyplot.savefig("DOP_histo.png", dpi = 1000)
    pyplot.show()
    
    #2nd power histo separated by H and D
    powerH = np.zeros(int(len(reformatted[8])/2))
    powerD = np.zeros(int(len(reformatted[8])/2))
    
    for c in range(int(len(reformatted[8])/2)):
        powerH[c] = reformatted[8][2*c]
        powerD[c] = reformatted[8][2*c+1]
        
    pyplot.hist(powerH, bins = 100, range = (0.0, 2e-5), label = "H polarized")
    pyplot.hist(powerD, bins = 100, range = (0.0, 2e-5), label = "D polarized")
    pyplot.title("Histogram of power of H and D polarized photons")
    pyplot.legend(loc='upper right')
    pyplot.savefig("power_histo.png", dpi = 1000)
    pyplot.show()
    
    inner_prod = np.zeros(int(len(reformatted[0])/2))
    for x in range(int(len(reformatted[0])/2)):
        first = np.array([reformatted[4][2*x], reformatted[5][2*x], reformatted[6][2*x]])
        second = np.array([reformatted[4][2*x + 1], reformatted[5][2*x + 1], reformatted[6][2*x + 1]])
        inner_prod[x] = np.dot(first, second)
    
    pyplot.hist(inner_prod, bins = 100)
    pyplot.title("Histogram of inner products for H and D vectrors")
    pyplot.savefig("innerProd.png")
    pyplot.show()

main()

