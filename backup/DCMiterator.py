# -*- coding: utf-8 -*-
"""
An integrator for orientation's direction cosine matrix. 

Theary is available in Chinese at whymranderson.blogspot.tw/2014/03/4th-runge-kutta.html

DCMiterator.DCMiter(wb,A0) takes a vector array wb(t1,t2,tN) (angular velocity along
body) whose shape = [N+1,3], eg. wb[t1,:] is wb(t1). And a initial rotation vector
A0 where np.dot(CK(A0),np.eye(3)) = initial orientation of A(t0) in
the space frame. CK(vector) is the CK rotation matrix originated from Rodriguez 
rotational formula. Returns subsequent orientation of A(t1,t2,...tN) where N = 
len(wb).
Created on Thu Jun 19 12:30:40 2014

@author: user
"""

import numpy as np

## Important for CK, quarternion is with its couterclockwise active sense here
## euler para rotational formula in counterclockwise active sense
def CK(rotvec): 
    amp = np.sqrt(np.dot(rotvec,rotvec))
    axis = rotvec/amp
    phi = amp % (2*np.pi)
    print(phi)
    a = np.cos(phi/2)
    b,c,d = axis*np.sin(phi/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def topEOM(wi,torquei):
    fderiv = np.dot(np.array([[0,b,0],[-b,0,0],[0,0,0]]),wi)+np.dot(
        np.array([[1/Ix,0,0],[0,1/Iy,0],[0,0,1/Iz]]),
        torquei)
    return fderiv
def topRK(wi,torquei):
    K1 = h*topEOM(wi,torquei)
    K2 = h*topEOM(wi+K1/2,torquei)
    K3 = h*topEOM(wi+K2/2,torquei)
    K4 = h*topEOM(wi+K3,torquei)
    nextw = wi + (K1+2*K2+2*K3+K4)/6
    return nextw

def DCMiter(wb,A0): # w_b are presolved from euler equations
# wb shape = [N+1,3 (xyz)]
    orien = np.array(A0)#starting orientation rotation vector from space xyz
    TrackMat_t0 = CK(orien)
    TrackMat_temp = TrackMat_t0  # shallow copy by reference
    N = len(wb[:,0])
    cordvec = np.zeros((N,3,3))     #(t(i), 012xyz components, 012xyz axes + other vertices)
    cordvec[0,:,2] = np.dot(CK(orien),np.array([0,0,1]))  #initial z_s(t0)
    cordvec[0,:,1] = np.dot(CK(orien),np.array([0,1,0]))  #initial y_s(t0)
    cordvec[0,:,0] = np.array([1,0,0])

    for i in range(1,N):
        TrackMat_temp = np.dot(TrackMat_temp,CK(wb[i,:]))        #wi here, 
        for j in range(3):
            cordvec[i,:,j]=np.dot(TrackMat_temp,np.eye(3)[j,:])
    return cordvec
    
def EulerDCMiter(w,Tau_body_temp,arm,h,w_option,cordvec,w_b,L_s,F,euler_angles_hasbun,N,TrackMat_temp): # w_b are solved simultaneously with euler eqs
    for i in range(1,N+1):
        w[i,:]=topRK(w[i-1,:],Tau_body_temp)
        #    print(TrackMat_temp)
        TrackMat_temp = np.dot(TrackMat_temp,CK((w_option[i,:])*h))        # +w_body_hasbun[i,:] wi here, _body_hasbun
        # TrackMat_temp is the transformation matrix for A_s(ti-1) to A_lab(ti)
        #    print(TrackMat_temp)
        for j in range(3):
            cordvec[i,:,j]=np.dot(TrackMat_temp,np.eye(3)[j,:])
        w_b[i-1,:] = np.dot(TrackMat_temp,w[i-1,:])
# if all rotation vectors are small enough, this block can replace the above--
#        transmatrix = np.dot(CK(orien),CK(h*np.sum(w[1:i+1,:],0)))      
#        cordvec[i,:,j] = np.dot(
#        np.transpose(transmatrix),
#        np.eye(3)[j,:])
#    print(CK(orien + h*np.sum(w[1:i+1,:])))
        Tau_s_temp = np.cross(arm*cordvec[i,:,2],F)             # change from cordvec(i-1) to i
        Tau_body_temp = np.dot(np.transpose(TrackMat_temp),Tau_s_temp)
        L_s[i,:]=Tau_s_temp*h + L_s[i-1,:]
#    time.sleep(5)
    for i in range(1,N):
        for j in range(3):  # calculate Hasbun's xyz axes in space
            cordvec[i,:,j+7]=np.dot(np.transpose(euler2space(euler_angles_hasbun[i,:])),np.eye(3)[j,:])
    return w_b,cordvec,L_s,w
