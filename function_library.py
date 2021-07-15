"""
Created on Tue Jul  6 09:07:10 2021

@author: liamm
"""

import numpy as np
from scipy.spatial.transform import Rotation as R
from math import cos,sin,radians,acos,asin,degrees,pi,sqrt
from py_pol.stokes import Stokes
from py_pol.utils import degrees
from py_pol.drawings import draw_poincare,draw_stokes_points,draw_empty_sphere
import random


"""a function that finds the paddle setting that best describes a given 
transformation. Requires the transformation to be given as a vector with the
direction describing the axis of rotation and the magnitude describes the angle
of rotation in radians.
INPUT FORMAT:
    <Paddle1>,<Paddle2>,<Paddle3>,<H_S1>,<H_S2>,<H_S3>,<D_S1>,<D_S2>,<D_S3>,<AXIS_S1>,<AXIS_S2>,<AXIS_S3>,<ROT_ANGLE>
RETURNS: numpy array with the paddle settings of the best transformation
       - [paddle1, paddle2, paddle3]
"""
def paddle_settings(transformation, infile_name):
    #takes in the paddle_transform_data.csv file as input or the coarse version.
    
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
	
""" Returns the Inverse of the input transformation. Uses the scipy.spatial.transform.Rotation library"""
def inverse(transform):
		return transform.inv()
	

"""Takes 2 input Stoke's vectors and returns the rotation that maps them onto the H and D states as a rotation vector
stokesA_prime maps to H, stokesB_prime maps to D"""

def findRotation(stokesA_prime,stokesB_prime):

    #H
    stokesA = Stokes("A")
    stokesA.from_components((1,1,0,0))

    #D
    stokesB = Stokes("B")
    stokesB.from_components((1,0,1,0))




    componentsA=stokesA.parameters.components()[1:]
    componentsB=stokesB.parameters.components()[1:]

    componentsA_prime=stokesA_prime.parameters.components()[1:]
    componentsB_prime=stokesB_prime.parameters.components()[1:]

    #Rotation 1
    rotation_angle=acos(np.dot(componentsA_prime,componentsA))/pi*180
    axis=np.cross(componentsA_prime,componentsA)
    axis=axis/np.linalg.norm(axis)

    rotation1=R.from_rotvec(axis*radians(rotation_angle))

    componentsA_primeprime=rotation1.apply(componentsA_prime)
    componentsB_primeprime=rotation1.apply(componentsB_prime)

    stokesA_primeprime=Stokes("A''")
    stokesB_primeprime=Stokes("B''")

    stokesA_primeprime.from_components(np.insert(componentsA_primeprime,0,1))
    stokesB_primeprime.from_components(np.insert(componentsB_primeprime,0,1))

    #Rotation 2

    rotation_angle=acos(np.dot(componentsB_primeprime,componentsB))/pi*180
    axis=np.cross(componentsB_primeprime,componentsB)
    axis=axis/np.linalg.norm(axis)

    if axis[0]>0:
        rotation2=R.from_rotvec(componentsA_primeprime*radians(rotation_angle))
    else:
        rotation2=R.from_rotvec(componentsA_primeprime*radians(-1*rotation_angle))


    componentsB_final=rotation2.apply(componentsB_primeprime)

    # A_final should be the same as A'', so this rotation should have no effect
    componentsA_final=rotation2.apply(componentsA_primeprime)


    stokesA_final=Stokes("A final")
    stokesB_final=Stokes("B final")


    stokesA_final.from_components(np.insert(componentsA_final,0,1))
    stokesB_final.from_components(np.insert(componentsB_final,0,1))



    total_rotation=rotation2*rotation1

    #Represent it as rotation vector (its axis is the axis of rotation and its length is angle of rotation in radians)
    rotvec=total_rotation.as_rotvec()
    print("Final rotation: rotate ",np.linalg.norm(rotvec)/pi*180," degrees around vector ",rotvec/np.linalg.norm(rotvec))

    return total_rotation
