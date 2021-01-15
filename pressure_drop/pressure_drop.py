# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 21:06:28 2021

@author: Darius Blaszyk
"""

import math
import random

def reynolds(u, dh, v):
    """
    Calculate Reynolds.
    
    Parameters
    ----------
    u : float
        Velocity based on the actial cross section of the duct or pipe [m/s]
    dh : float
        Hydraulic diameter [m]
    v : float
        Kinematic viscosity of the media [m2/s]
    """
    return u * dh / v

def ColebrookWhite(D, Re, k, eps = 0.0001, MaxIter = 1000):
    """
    Calculates the Darcy-Weisbach friction factor for pressure loss calculations
    via solving the implicit Colebrook&White expression, details in https://doi.org/10.1098/rspa.1937.0150
    by Tol,Hakan Ibrahim from the PhD study at Technical University of Denmark
    PhD Topic: District Heating in Areas with Low-Energy Houses
    
    Parameters
    ----------
    D : float
        Inner diameter of the pipe [mm]
    Re : float
        Reynolds Number [-]
    k : float
        Absolute roughness of pipe [mm]
    eps : float
        Termination Tolerance(Iteration) [-]
    MaxIter : float
        Max. limit (Iteration) [-]
    """

    #Initializing the Iteration
    Err = 10    #Iteration error
    IterNum = 0 #Iteration steps number
    
    #Initial estimate (making use of SwameeJain algorithm)
    X0 = random.random() #Random starting point
    X0 = 0.5
    #X0 = f_Clamond(Re, aRou / D)
    
    #Fasten your seat belts, iteration starts
    
    while ((Err > eps) and (IterNum < MaxIter)):
        IterNum = IterNum + 1
        X1 = (2 * math.log10((k / D) / 3.7 + 2.51 / (Re * X0 ** 0.5))) ** -2
        Err = abs(X1 - X0)
        X0 = X1
    
    return X1

# calcuate pressure loss in pipe via the Darcy-Weisbach equation
def pressure_loss(length, diam, flow, dens, k, v):
    """
    Calculate pressure loss for a rough pipe.
    
    Parameters
    ----------
    length : float
        Length of the pipe [m]
    diam : float
        Hydraulic diameter of the pipe [mm]
    flow : float
        Volumetric flow thorugh the pipe [m3/hr]
    dens : float
        Media density [kg/m3]
    k : float
        Roughness of pipe [mm]
    v : float
        Kinematic viscosity of the media [m2/s] 
    """
    
    O = ((diam / 1000) / 2) ** 2 * math.pi
    u = flow / 3600 / O
    
    Re = reynolds(u, diam / 1000, v)
    
    # Re<2000, we use λ=64/Re
    # Re>2500, we use λ=Darcy-Weisbach
    # linear interpolation in the area 2000-2500
    if Re < 2000:
        labda = 64 / Re
    else:
        if Re > 2500:
            labda = ColebrookWhite(diam / 1000, Re, k)
        else:
            labda_2500 = ColebrookWhite(diam / 1000, 2500, k)
            labda_2000 = 64 / Re
            
            labda = (labda_2500 - labda_2000)/(2500 - 2000) * (Re - 2000) + labda_2000
    
    return labda * (length / (diam / 1000)) * (dens * u ** 2 / 2)


# print reynolds
Re = reynolds(
    0.282942, 
    0.05,
    0.00000113)
print('Re =', Re)

print('λ =', ColebrookWhite(0.05, Re, 0.0001))

# print pressure loss
dp = pressure_loss(
    100,
    50,
    2,
    998.3,
    0.0001,
    0.00000113)

print('Δp =', dp)