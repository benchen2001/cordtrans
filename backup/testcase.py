# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 11:55:19 2014

@author: user
"""

import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d as p3

def CK(rotvec): #passive, r.h. rule
    amp = np.sqrt(np.dot(rotvec,rotvec))
#    print(amp)
    axis = -rotvec/amp
    phi = amp % (2*np.pi)
    a = np.cos(phi/2)
    b,c,d = axis*np.sin(phi/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])


A = np.array([1,0,0])
B = np.array([0,1,0])
C = np.array([0,0,1])
print(np.dot(CK(A+B),C))
k=np.dot(CK(A),CK(B))
print(np.dot(k,C))


'''
omega1 = np.array([1,0,0])*np.pi/180*1

testmat = CK(omega1)
print(testmat)
out = np.dot(testmat,np.array([0,1,1]))

fig=plt.figure()
ax = p3.Axes3D(fig)
ax.view_init(elev=20, azim=15)

lineL=ax.plot([0,out[0]],
              [0,out[1]],
              [0,out[2]],'k-')

plt.show()
'''