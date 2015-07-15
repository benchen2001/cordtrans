# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 11:43:46 2014
A rigid body rotation solver for angular velocity along body. 
Input: 
Output: 
@author: user
"""

import numpy as np

def coinEOMgeneral(wi,torquei_along_body):
    InerTen = np.array([[Ix,0,0],[0,Iy,0],[0,0,Iz]])
    angL=np.dot(InerTen,wi)
    crossprod = np.cross(wi,angL)
    fderiv = torquei_along_body - crossprod
    fderiv = np.dot(np.array([[1/Ix,0,0],
                              [0,1/Iy,0],
                              [0,0,1/Iz]]), fderiv)
    return fderiv

def topRK(wi,torquei,h):
    K1 = h*coinEOMgeneral(wi,torquei)
    K2 = h*coinEOMgeneral(wi+K1/2,torquei)
    K3 = h*coinEOMgeneral(wi+K2/2,torquei)
    K4 = h*coinEOMgeneral(wi+K3,torquei)
    nextw = wi + (K1+2*K2+2*K3+K4)/6
    return nextw
    
tn=1.55                    # end of simulation time, take longer time if tn > 10 seconds
t0=0                        # start of time zero
samplerate = 3000           # rate of iteration in Hz
N = int(round((tn-t0)*samplerate))    #number of steps
h = (tn-t0)/N 

w = np.zeros([N+1,3]) 

Iy= 0.002   #moment of inertial, substitute disk formula if needed Iy=0.25*M*R**2 + 1/12*M*L**2 +M*arm**2
Iz= 0.0008  #Iz = 0.5*M*R**2
Ix= 0.002


for j in range(N):
    Tau_s_j = np.cross(R*C[j,:,1],np.array([0,0,-m*g]))
    Tau_b_j = np.dot(C[j,:,:],Tau_s_j)