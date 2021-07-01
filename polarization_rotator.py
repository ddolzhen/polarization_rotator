import numpy as np
from math import cos,sin,radians,acos,asin,degrees,pi,sqrt
from py_pol.stokes import Stokes
from py_pol.utils import degrees
from py_pol.drawings import draw_poincare,draw_stokes_points,draw_empty_sphere
from scipy.spatial.transform import Rotation as R
import random



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