# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 09:45:35 2021

@author: liamm
"""

import sys
import numpy as np
import matplotlib.pyplot as pyplot

from math import radians,acos,pi,sqrt
from py_pol.stokes import Stokes
from scipy.spatial.transform import Rotation as R

from polarization_rotator import findRotation

import warnings
warnings.filterwarnings("ignore")


    
# Below is just data_reform.py

if len(sys.argv) != 4:
    print("the program requires 3 arguments")
    print("format: <program_file> <input_file> <reform_file> <paddle-transform_file>")
    quit()
else: 
    infile_name = sys.argv[1]
    reformfile_name = sys.argv[2]
    paddlefile_name = sys.argv[3]
    
infile = open(infile_name, "r")

content = infile.readlines()

del content[0]


for i in range(len(content)):
    content[i] = content[i].split(",")
    

out_data = np.zeros(shape = (len(content[1]), int(len(content)/4 )) )
std_dev_list = np.array([])

for j in range(len(content[0])):
    for k in range(int(len(content)/4)):
        data = np.array([float(content[4*k][j]), float(content[4*k+1][j]), float(content[4*k+2][j]), float(content[4*k+3][j])])
        out_data[j][k] = np.mean(data) 
        if j == 4 or j == 5 or j == 6:
            std_dev_list = np.append(std_dev_list, np.std(data))
        
        
"""
#out_data has format:
    - each row is one of the variables (rotator, s1, s2, etc...)
    -each column is corresponding value for that run.

"""
#define the resulting data arrays

#first order the data based on paddle_3, then paddle_2, then paddle_1 settings
sorting = np.lexsort((out_data[3], out_data[2], out_data[1]))

reformatted = np.ndarray(shape = (len(out_data), len(out_data[0])))
for l in range(len(out_data)):
    for m in range(len(sorting)):
        reformatted[l][m] = out_data[l][sorting[m]]
        


outfile = open(reformfile_name, "w")
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

pyplot.savefig("DOP_histo-09-07.png", dpi = 1000)
pyplot.show()

#histo of std dev between stokes vector components between repeat trials.
pyplot.hist(std_dev_list, bins = 100, range = (0, 0.02))

pyplot.title("Histogram of Standard Deviation of Stokes vectors between trials")

pyplot.savefig("STD_histo-09-07.png", dpi = 1000)
pyplot.show()


#2nd power histo separated by H and D
powerH = np.zeros(int(len(reformatted[8])/2))
powerD = np.zeros(int(len(reformatted[8])/2))

for c in range(int(len(reformatted[8])/2)):
    powerH[c] = reformatted[8][2*c]
    powerD[c] = reformatted[8][2*c+1]
    
pyplot.hist(powerH, bins = 100, range = (0.0, 3e-5), label = "H polarized")
pyplot.hist(powerD, bins = 100, range = (0.0, 3e-5), label = "D polarized")
pyplot.title("Histogram of power of H and D polarized photons")
pyplot.legend(loc='upper right')
pyplot.xlabel("Power / W")
pyplot.savefig("power_histo-09-07.png", dpi = 1000)
pyplot.show()

inner_prod = np.zeros(int(len(reformatted[0])/2))
for x in range(int(len(reformatted[0])/2)):
    first = np.array([reformatted[4][2*x], reformatted[5][2*x], reformatted[6][2*x]])
    second = np.array([reformatted[4][2*x + 1], reformatted[5][2*x + 1], reformatted[6][2*x + 1]])
    inner_prod[x] = np.dot(first, second)

pyplot.hist(inner_prod, bins = 100)
pyplot.title("Histogram of inner products for H and D vectrors")
pyplot.savefig("innerProd-09-07.png")
pyplot.show()







#Below is just paddle-transform_mapper.py


#Define desired polarization points



#read input data then apply the transformation line by line
infile = open(reformfile_name, "r")
    
content1 = infile.readlines()

if content1[0][0] == "#":
    del content1[0]

    
for i in range(len(content1)):
    content1[i] = content1[i].split(",")
     
#prepare output files for writing
out_file = open(paddlefile_name, "w")
out_file.write("#Paddle1,#Paddle2,#Paddle3,#H_S1,#H_S2,#H_S3,#D_S1,#D_S2,#D_S3,#AXIS_S1,#AXIS_S2,#AXIS_S3,#ROT_ANGLE\n")



#This block of code is repeated for all polarization pairs in the data file.
for lines in range(int(len(content1)/2)):
    
    #H
    stokesH = Stokes("H")
    stokesH.from_components((1,1,0,0))
    
    #D
    stokesD = Stokes("D")
    stokesD.from_components((1,0,1,0))
    
    componentsH = stokesH.parameters.components()[1:]
    componentsD = stokesD.parameters.components()[1:]
    
    
    
    #define out_data for this line
    out_data = np.zeros(13)
    #set out_data paddle settings to those of infile
    out_data[0],out_data[1],out_data[2] = content1[2*lines][1], content1[2*lines][2], content1[2*lines][3] 
    
    """Stokes Parameters from observed values in the Data files"""

    #Find stokes components of H
    H_s1= float(content1[2*lines][4])
    H_s2= float(content1[2*lines][5])
    H_s3= float(content1[2*lines][6])
    
    
    #Find stokes components of D
    D_s1 = float(content1[2*lines+1][4])
    D_s2 = float(content1[2*lines+1][5])
    D_s3 = float(content1[2*lines+1][6])
    
    
    
    #Define A' and B' from files, these are then corrected to be orthogonal.
    stokesA=Stokes("A")
    stokesA.from_components((sqrt(H_s1**2+H_s2**2+H_s3**2),H_s1,H_s2,H_s3)).normalize()
    
    
    #Define B' as orthagonal to A'
    stokesB=Stokes("B")
    stokesB.from_components((sqrt(D_s1**2+D_s2**2+D_s3**2),D_s1,D_s2,D_s3)).normalize()
    
    
    #Get stokes components A and B
    componentsA=stokesA.parameters.components()[1:]
    componentsB=stokesB.parameters.components()[1:]
    
    #How far are A and B from being perpendicular, divide by 2 as this is applied to both vectors
    rot_angle = (acos(np.dot(componentsA, componentsB))/pi*180 - 90)/2
    
    #axis perpendicular to plane of A and B. Later used for determining direction of rotation 2.
    rot_axis = np.cross(componentsA, componentsB)
    rot_axis /= np.linalg.norm(rot_axis)
    
    #defining rotations and applying them to A and B, updating stokes vectors
    rot_A = R.from_rotvec(rot_axis*radians(rot_angle))
    rot_B = R.from_rotvec(rot_axis*radians(-1*rot_angle))
     
    componentsA = rot_A.apply(componentsA)
    componentsB = rot_B.apply(componentsB)
    
    
    
    stokesA.from_components(np.insert(componentsA,0,1))
    stokesB.from_components(np.insert(componentsB,0,1))
    
    #write in these starting values to an output file.
    out_data[3],out_data[4],out_data[5] = componentsA[0], componentsA[1], componentsA[2]
    out_data[6],out_data[7],out_data[8] = componentsB[0], componentsB[1], componentsB[2]
    
    
    
    
    total_rotation = findRotation(stokesA, stokesB)
    
    #Represent it as rotation vector (its axis is the axis of rotation and its length is angle of rotation in radians)
    rotvec = total_rotation.as_rotvec()
    angle = np.linalg.norm(rotvec)/pi*180
    rot_axis = rotvec/np.linalg.norm(rotvec)
    
    
    out_data[9], out_data[10], out_data[11], out_data[12] = str(rot_axis[0]), str(rot_axis[1]), str(rot_axis[2]), str(angle)
    
    #Now use total rotation to put H and D into A' and B' respectively
    #We are now using 1 combined rotation instead of 2 separate rotations, but the effect should be the same 
    
    componentsH_final = total_rotation.apply(componentsH)
    componentsD_final = total_rotation.apply(componentsD)
    
    stokesH_final=Stokes("H final2")
    stokesD_final=Stokes("D final2")
    
    
    stokesH_final.from_components(np.insert(componentsH_final,0,1))
    stokesD_final.from_components(np.insert(componentsD_final,0,1))

    
    
    out_file.write("{0:f},{1:f},{2:f},{3:f},{4:f},{5:f},{6:f},{7:f},{8:f},{9:f},{10:f},{11:f},{12:f}\n".format(out_data[0],out_data[1],out_data[2],out_data[3],
                                                                                                               out_data[4],out_data[5],out_data[6],out_data[7],
                                                                                                               out_data[8],out_data[9],out_data[10],out_data[11],
                                                                                                            out_data[12]))
    
out_file.close()