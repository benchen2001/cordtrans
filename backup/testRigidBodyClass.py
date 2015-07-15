# -*- coding: utf-8 -*-
"""
Created on Fri Nov 07 15:15:20 2014

@author: user
"""

import numpy as np

class RigidBodyObject:
    def __init__(self):
        self.data=[]

    ##_________ Parameters of top
    description = "Description not defined yet"
    M = 1       #mass in kg
    #R = 0.025  #radius of disk in m
    #L = 0.01   #width of disk in m
    arm = 0.04  #level arm of Center of mass in m
    Iy= 0.002   #moment of inertial, substitute disk formula if needed Iy=0.25*M*R**2 + 1/12*M*L**2 +M*arm**2
    Iz= 0.0008  #Iz = 0.5*M*R**2
    Ix= 0.002
    g = 9.8                     #gravity m/s^2
    F = np.array([0,0,-M*g])    #gravity f
    freq = 15                 #top revolution speed in hertz, along symmetric axis
    tn=1                      # end of simulation time, take longer time if tn > 10 seconds
    t0=0.0                      # start of time zero
    samplerate = 2000           # rate of iteration in Hz
    N = int(round((tn-t0)*samplerate))    #number of steps
    h = (tn-t0)/N

    #__________________________________________ Initial condition
    orien = np.array([-np.pi/3,0,0])            #starting orientation vector of top from lab xyz
    #orien = np.array([-np.radians(5),0,0])     #starting orientation turn of top from space xyz
    w = np.zeros([N+1,3])                       # w_s
    w_b = np.zeros([N+1,3])                     # w_lab
    w[0,2]= 2*np.pi*freq                        #wz is the symetry axis



a = RigidBodyObject
a.M=2
print(a.description)