# -*- coding: utf-8 -*-
"""
Created on Thu Mar 05 11:29:52 2015

@author: user
"""
import numpy as np

N = 5
omega = np.array([[ 1.0, 2, 2],
                 [ -1, 1,0],
                 [ 2,  2,  3],
                 [ -1,  1,  0],
                 [ 1,  2 , 3]])
h = 0.1
#print omega
#alpha = np.zeros((N,3))

#def alphafunc(omegalist):    
#    for i in range(N):
#        alpha[i,:]=np.sum(omega[:i+1,:],axis=0)
#alphafunc(omega)
alpha=np.array([np.sum(omega[:i+1,:],axis=0)*h for i in range(N)])
#for aa in range(N):
#    print np.cross(alpha[aa,:],omega[aa,:])
deltaalpha = np.sum([np.cross(alpha[i,:],omega[i,:]) for i in range(N)],axis=0)*h
print deltaalpha

rho = alpha[-1,:] + deltaalpha
print rho