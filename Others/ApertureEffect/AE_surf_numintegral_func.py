import os
import obspy
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
import math
import pandas as pd
from tqdm import tqdm
import warnings
from scipy import integrate
import scipy

import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

def stf_cosine(t, TR, fz):
    '''
        source time function of cosine wavelet
        https://tktmyd.github.io/OpenSWPC/English/2._Parameter_Settings/0207_source/
        Argument: 
            t:: time vector
    '''
    stf = np.zeros(len(t))
    for i, tt in enumerate(t):
        if 0<tt and tt<TR:
            stf[i] = (fz/TR) * (1 - np.cos((2*np.pi*tt/TR)))
        else:
            stf[i] = 0
            
    return stf
    
def y_cosine_analytic(t, v, theta, TR, R):
    y = np.zeros(len(t))
    for i, tt in enumerate(t):
        if theta==0:
            y[i] = np.cos(2*np.pi*tt/TR)
        else:
            va = v/np.sin(theta)
            y[i] = ((va*TR)/(np.pi*R)) * np.cos(theta) * mp.besselj(1, (2*np.pi*R)/(va*TR)) *  np.cos(2*np.pi*tt/TR)

    return y


def y_cosine_numeric(t, v, theta, TR, R):
    y = np.zeros(len(t))
    for i, tt in enumerate(t):
        if theta==0:
            y[i] = np.cos(2*np.pi*tt/TR)
        else:
            va = v/np.sin(theta)
            # compute integral
            f_int = lambda phi, r: np.cos(2*np.pi*((tt + r*np.cos(phi)/va)/TR)) * r
            y[i] = (np.cos(theta)/(np.pi*R**2))*integrate.dblquad(f_int, 0, R, lambda x: 0, lambda x: 2*np.pi)[0]
            
    return y    


# Compute amplitude response with incident angle
def incidentangle_scalingfactor_analytic(v, theta, TR, R):
    if theta==0:
        return 1.0
    else:
        va = v/np.sin(theta)
        J1 = mp.besselj(1, (2*np.pi*R)/(va*TR))
        return  ((va * TR)/(np.pi*R)) * J1


def stf_herzian_mclaskey2009(t, rho, R, v, E1, nu1, E2, nu2):
    '''
        source time function of Herzian solution
    '''
    def get_delta(E1, mu1, E2, mu2):
        del1 = (1-nu1**2)/(np.pi*E1)
        del2 = (1-nu2**2)/(np.pi*E2)
        return (del1, del2)
    
    
    def get_contact_time(rho, del1, del2, R, v):
        return 4.53*(( 4*rho  * np.pi * (del1 + del2) /3 )**(2/5)) * R * (v **(-1/5))

    def get_maximum_force(rho, del1, del2, R, v):
        return 1.917 * (rho**(3/5)) * ((del1 + del2)**(-2/5)) * (R**2) * (v**(6/5))

        
    del1, del2 = get_delta(E1, nu1, E2, nu2)
    
#     print(del1, del2, rho, R, v)
    
    tc = get_contact_time(rho, del1, del2, R, v)
    fmax = get_maximum_force(rho, del1, del2, R, v)
#     print(tc, fmax)
    
    stf = np.zeros(len(t))
    for i, tt in enumerate(t):
        if 0<tt and tt < tc:
            stf[i] = fmax*(np.sin(np.pi*tt/tc))**(3/2)
        else:
            stf[i] = 0
            
    return stf, fmax, tc

# Ball constants
rho1 = 7781.1
cp  = 5900
cs  = 3230

#shear modulus
G1 = rho1 * (cs **2) #[Pa]
# poisson's ratio of rock
vpvs = cp/cs
nu1  = 0.5 * (((vpvs**2) - 2.0) / ((vpvs**2) - 1.0))
# Young's modulus of rock
E1 = 2.0 * G1 * (1+nu1) # [Pa]

print(E1/1e9, nu1)

# Rock block constants
rho2 = 2980
cp  = 6200 #6919
cs  = 3590 #3631

#shear modulus
G2 = rho2 * (cs **2) #[Pa]
# poisson's ratio of rock
vpvs = cp/cs
nu2  = 0.5 * (((vpvs**2) - 2.0) / ((vpvs**2) - 1.0))
# Young's modulus of rock
E2 = 2.0 * G2 * (1+nu2) # [Pa]

print(E2/1e9, nu2, G2/1e9)

R_ball = 3.0e-3/2# 0.5e-3#(3.0e-3) /2#(3.0e-3)/2 # ball with diameter of 1.5mm:  for tha case of 4m ball drop.
h = 0.500 # 0.127 #0.500  # [m] dropped 500mm for tha case of 4m ball drop.
g = 9.80665
v_impact = np.sqrt(2*g*h)


def stf_balldrop(t):
    return stf_herzian_mclaskey2009(t, rho1, R_ball, v_impact, E1, nu1, E2, nu2)

def y_balldrop_numeric(t, v, theta, R):
    y = np.zeros(len(t))
    for i, t1 in enumerate(tqdm(t)):
        if theta==0:
            y[i] = stf_balldrop([t1])[0]
        else:
            va = v/np.sin(theta)
            # compute integral
            f_int = lambda phi, r: r*stf_balldrop([t1 + (r*np.cos(phi)/va)])[0]
            y[i] = (np.cos(theta)/(np.pi*R**2))*integrate.dblquad(f_int, 0, R, lambda x: 0, lambda x: 2*np.pi)[0] 
            
    return y    



def stf_herzian_mclaskey2009_scaled(t, TR):
    stf = np.zeros(len(t))
    for i, tt in enumerate(t):
        if 0<tt and tt < TR:
            stf[i] = (np.sin(np.pi*tt/TR))**(3/2)
        else:
            stf[i] = 0
            
    return stf

def y_balldrop_scaled_numeric(t, v, theta, TR, R):
    y = np.zeros(len(t))
    for i, t1 in enumerate(t):
        if theta==0:
            y[i] = stf_herzian_mclaskey2009_scaled([t1], TR)
        else:
            va = v/np.sin(theta)
            # compute integral
            f_int = lambda phi, r: r*stf_herzian_mclaskey2009_scaled([t1 + (r*np.cos(phi)/va)], TR)
            y[i] = (np.cos(theta)/(np.pi*R**2))*integrate.dblquad(f_int, 0, R, lambda x: 0, lambda x: 2*np.pi)[0]
            
    return y    


# Compute amplitude response with incident angle
def incidentangle_scalingfactor_balldrop(v, theta, TR, R):
    if theta==0:
        return 1.0
    else:
        tmax_balldrop = TR/2
        return y_balldrop_scaled_numeric([tmax_balldrop], v, theta, TR, R)/np.cos(theta)


#---------------------------------------------------#
# Functions for fitting scaling factor
#---------------------------------------------------#
# def f_obj(param, A, theta_scaling, v, R):
#     """
#     Objective function to find best power coefficient of fitting.
#     param: Array with [TR, b] to be optimized
#     """
    
#     # compute analytical response
#     if len(A) != len(theta_scaling):
#         print("error: length of theta scaling must be equal to the size of A.")
    
#     TR, b = param
#     k_scale_base = np.zeros(len(theta_scaling))
#     for i, theta1 in enumerate(theta_scaling):
#         k_scale_base[i] = incidentangle_scalingfactor_analytic(v, theta1, TR, R)

#     return np.linalg.norm(A-k_scale_base**b, ord=2)

def f_obj(param, A, theta_scaling, v, R):
    """
    Objective function to find best power coefficient of fitting.
    hatTR to be optimized
    """
    hatTR = param[0]
    
    # compute analytical response
    if len(A) != len(theta_scaling):
        print("error: length of theta scaling must be equal to the size of A.")
    
    k_scale_base = np.zeros(len(theta_scaling))
    for i, theta1 in enumerate(theta_scaling):
        k_scale_base[i] = incidentangle_scalingfactor_analytic(v, theta1, hatTR, R)

    return np.linalg.norm(A-k_scale_base, ord=2)


def f_obj_debug(param, A, theta_scaling, v, R):
    """
    Objective function to find best power coefficient of fitting.
    param: Array with [TR, b] to be optimized
    """
    
    # compute analytical response
    if len(A) != len(theta_scaling):
        print("error: length of theta scaling must be equal to the size of A.")
    
    TR, b = param
    k_scale_base = np.zeros(len(theta_scaling))
    for i, theta1 in enumerate(theta_scaling):
        k_scale_base[i] = incidentangle_scalingfactor_analytic(v, theta1, TR, R)

    return A, k_scale_base #**b


def incidentangle_scalingfactor_fitwithanalytic(v, theta, TR, R):
    if theta==0:
        return 1.0
    else:
        return incidentangle_scalingfactor_analytic(v, theta, TR, R)



