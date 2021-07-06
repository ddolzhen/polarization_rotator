# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 09:28:38 2021

assigns a transformation (rotation angle and axis of rotation) for a 
particular paddle setting by reading in the input data.

INPUT FORMAT: 
    <ROTATOR>,<PADDLE1>,<PADDLE2>,<PADDLE3>,<S1>,<S2>,<S3>,<DOP>,<POWER>
    -lines come paired up with same paddle settings but different rotator settings
    -first line is designated H, second is D.
    
OUTPUT FORMAT: 
    <PADDLE1>,<PADDLE2>,<PADDLE3>,<H_S1>,<H_S2>,<H_S3>,<D_S1>,<D_S2>,<D_S3>,<ROT_AXIS1>,<ROT_AXIS2>,<ROT_AXIS3>,<ROT_ANGLE>


@author: liamm
"""

import sys
import numpy as np
from math import radians,acos,pi,sqrt
from py_pol.stokes import Stokes
from scipy.spatial.transform import Rotation as R


from polarization_rotator import findRotation

import warnings
warnings.filterwarnings("ignore")

#Define desired polarization points


#takes in 3 file names as input
# 1) reformatted data file as input
# 2) paddle_transform_data.csv
# 3) error_file.csv
if len(sys.argv) != 3:
    print("the program requires 3 arguments")
    print("format: <program_file> <input_file> <output_file>")
    quit()
else: 
    #infile is reformatted data file
    infile_name = sys.argv[1]
    outfile_name = sys.argv[2]

#read input data then apply the transformation line by line
infile = open(infile_name, "r")
    
content = infile.readlines()

if content[0][0] == "#":
    del content[0]

    
for i in range(len(content)):
    content[i] = content[i].split(",")
     
#prepare output files for writing
out_file = open(outfile_name, "w")
out_file.write("#Paddle1,#Paddle2,#Paddle3,#H_S1,#H_S2,#H_S3,#D_S1,#D_S2,#D_S3,#AXIS_S1,#AXIS_S2,#AXIS_S3,#ROT_ANGLE\n")



#This block of code is repeated for all polarization pairs in the data file.
for lines in range(int(len(content)/2)):
    
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
    out_data[0],out_data[1],out_data[2] = content[2*lines][1], content[2*lines][2], content[2*lines][3] 
    
    """Stokes Parameters from observed values in the Data files"""

    #Find stokes components of H
    H_s1= float(content[2*lines][4])
    H_s2= float(content[2*lines][5])
    H_s3= float(content[2*lines][6])
    
    
    #Find stokes components of D
    D_s1 = float(content[2*lines+1][4])
    D_s2 = float(content[2*lines+1][5])
    D_s3 = float(content[2*lines+1][6])
    
    
    
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
    
    
    out_data[3],out_data[4],out_data[5] = componentsA[0], componentsA[1], componentsA[2]
    out_data[6],out_data[7],out_data[8] = componentsB[0], componentsB[1], componentsB[2]
    
    
    
    
    total_rotation = findRotation(stokesA, stokesB)
    
    #Represent it as rotation vector (its axis is the axis of rotation and its length is angle of rotation in radians)
    rotvec = total_rotation.as_rotvec()
    angle = np.linalg.norm(rotvec)/pi*180
    rot_axis = rotvec/np.linalg.norm(rotvec)
    #print("Final rotation: rotate ",angle ," degrees around vector ",rot_axis)
    
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