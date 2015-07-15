import matplotlib.pyplot as plt
import numpy as np

M = 1   #kg
R = 0.025#m
L = 0.01     #m width of disk
arm = 0.1   #m
Iy= 0.25*M*R**2 + 1/12*M*L**2 +M*arm**2
Ix= 0.5*M*R**2
Iz= Iy
g = 9.8 #m/s^2

Nx = 0
tileangle = 45  #in radian
Ny = M*g*np.cos(np.radians(45))*arm
Nz = 0

a1 = Ny/Iy
a2 = 0
wx = 2*np.pi*25
b = (Iy/Ix - 1)*Ix*wx/Iy
print(a1,a2,wx,b)
N = 3000
tn=1
t0=0

def augfun(vector):
    fout = np.dot(np.array([[0,b],[-b,0]]),vector)+ np.array([a1,a2]).reshape((2,1))
    return fout


h = (tn-t0)/N
xt0 = 0     #wy intial value
yt0 = 0     #wz initial value

vector = np.array([xt0,yt0]).reshape((2,1))
t = t0
veclist = np.zeros((2,N))
tlist = np.zeros((1,N))
wxvec = np.zeros(N) + wx
for i in range(1,N):
    K1 = h*augfun(vector)
    K2 = h*augfun(vector+K1/2)
    K3 = h*augfun(vector+K2/2)
    K4 = h*augfun(vector+K3)

    vector = vector + (K1+2*K2+2*K3+K4)/6
    t = t0+i*h
    veclist[:,i] = vector.reshape(1,2)
    tlist[:,i]= t

Ltotalinbody = (Ix*wxvec)**2 + (Iy*veclist[0,:])**2+(Iz*veclist[1,:])**2
print(max(Ltotalinbody)-min(Ltotalinbody))
plt.figure()
plt.plot(tlist[0,:],veclist[0,:], 'r-')
plt.plot(tlist[0,:],veclist[1,:], 'b-')
plt.plot(tlist[0,:],wxvec, 'k-')
plt.plot(tlist[0,:],Ltotalinbody, 'k-')
plt.show()

np.savez('workfile', wp=veclist, wx=wxvec, deltat=h, Ix=Ix,Iy=Iy,Iz=Iz)

