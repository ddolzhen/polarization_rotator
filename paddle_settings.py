# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 09:07:10 2021

a function that finds the paddle setting that best describes a given 
transformation. Requires the transformation to be given as a vector with the
direction describing the axis of rotation and the magnitude describes the angle
of rotation in radians.

INPUT FORMAT:
    <Paddle1>,<Paddle2>,<Paddle3>,<H_S1>,<H_S2>,<H_S3>,<D_S1>,<D_S2>,<D_S3>,<AXIS_S1>,<AXIS_S2>,<AXIS_S3>,<ROT_ANGLE>

RETURNS: numpy array with the paddle settings of the best transformation
       - [paddle1, paddle2, paddle3]
@author: liamm
"""

import sys
import numpy as np
from math import radians

def paddle_settings(transformation):
    #takes in the paddle_transform_data.csv file as input or the coarse version.
    if len(sys.argv) != 2:
        print("the program requires 2 arguments")
        print("format: <program_file> <input_file>")
        quit()
    else: 
        #infile is reformatted data file
        infile_name = sys.argv[1]
    
    #read input data then apply the transformation line by line
    infile = open(infile_name, "r")
        
    content = infile.readlines()
    
    if content[0][0] == "#":
        del content[0]
        
    for i in range(len(content)):
        content[i] = content[i].split(",")
          
        
    #define the closest vector to the given transform
    closest_transform = np.zeros(3)
    transform_diff = np.linalg.norm(transformation - closest_transform)
    #Keep track of the best paddle settings
    paddles = np.zeros(3)
    
    for lines in range(len(content)):
        #compare proximity of axes by subtraction
        
        #makes a transform vector in the same form as the input based on input data
        new_transform = radians(float(content[lines][12]))*np.array([float(content[lines][9]), float(content[lines][10]), float(content[lines][11])])
        if np.linalg.norm(transformation - new_transform) < transform_diff:
            closest_transform = new_transform
            transform_diff = np.linalg.norm(transformation - closest_transform)
            paddles = np.array([float(content[lines][0]), float(content[lines][1]), float(content[lines][2])])
         
    return paddles