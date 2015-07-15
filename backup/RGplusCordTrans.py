import matplotlib.pyplot as plt
import numpy as np
import numpy as np
from matplotlib import pyplot
import matplotlib.gridspec as gridspec
import pylab
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

def topEOM(wi,torquei):
    fout = np.dot(np.array([[0,b,0],[-b,0,0],[0,0,0]]),wi)+np.dot(
        np.array([[1/Ix,0,0],[0,1/Iy,0],[0,0,1/Iz]]),
        torquei)
    return fout
def topRK(wi,torquei):
    K1 = h*topEOM(wi,torquei)
    K2 = h*topEOM(wi+K1/2,torquei)
    K3 = h*topEOM(wi+K2/2,torquei)
    K4 = h*topEOM(wi+K3,torquei)
    nextw = wi + (K1+2*K2+2*K3+K4)/6
    return nextw
##_______________________ Important, quarternion is with its couterclockwise active sense here
def CK(rotvec): ##euler para rotational formula in counterclockwise active sense
    amp = np.sqrt(np.dot(rotvec,rotvec))
    axis = rotvec/amp
    phi = amp*h
    a = np.cos(phi/2)
    b,c,d = axis*np.sin(phi/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

##_______________________ Parameters of top
M = 1   #kg
R = 0.025#m
L = 0.01     #m width of disk
arm = 0.1   #m
Iy= 0.25*M*R**2 + 1/12*M*L**2 +M*arm**2
Iz= 0.5*M*R**2
Ix= Iy

g = 2 #m/s^2
F = np.array([0,0,-M*g])
freq = 50 #hertz
N = 30000 #number of steps
tn=5    #in seconds
t0=0
h = (tn-t0)/N
##_______________________ Initial condition
orien = np.array([-np.pi/4,0,0])     #starting orientation turn of top from space xyz
w = np.zeros([N+1,3])
totalKE = np.zeros([N+1])
w[:,2]= w[:,2]-2*np.pi*freq                #wz is the symetry axis
#w[0,0:2]= np.zeros(2)              #initial w(t0)
b=-(Iz-Iy)/Ix*w[0,2]
cordvec = np.zeros((N+1,3,3))     #(i element, 012xyz components, 012xyz cord vectors)
cordvec[0,:,2] = np.dot(CK(orien/h),np.array([0,0,1]))     #initial z_s(t0)
cordvec[0,:,1] = np.dot(CK(orien/h),np.array([0,1,0]))     #initial y_s(t0)
cordvec[0,:,0] = np.array([1,0,0])                       #initial x_s(t0)
Tau_s_t0 = np.cross(arm*cordvec[0,:,2],F)
TrackMat_t0 = CK(orien/h)
Ltotalinbody,Ltotalinspace = np.zeros((N,3)),np.zeros((N,3))

print(np.linalg.norm(Tau_s_t0)/np.linalg.norm(np.array([Ix,Iy,Iz]))/2/np.pi)

##____________________ main calculation
Tau_s_temp = Tau_s_t0
TrackMat_temp = TrackMat_t0
Tau_body_t0 = np.dot(np.transpose(TrackMat_temp),Tau_s_t0)
Tau_body_temp = Tau_body_t0

for i in range(1,N+1):
    w[i,:]=topRK(w[i-1,:],Tau_body_temp)
    TrackMat_temp = np.dot(TrackMat_temp,CK(w[i,:]))
    cordvec[i,:,0]=np.dot(TrackMat_temp,np.array([1,0,0]))
    cordvec[i,:,1]=np.dot(TrackMat_temp,np.array([0,1,0]))
    cordvec[i,:,2]=np.dot(TrackMat_temp,np.array([0,0,1]))
    Tau_s_temp = np.cross(arm*cordvec[i-1,:,2],F)
    Tau_body_temp = np.dot(np.transpose(TrackMat_temp),Tau_s_temp)
    totalKE[i] = np.dot(np.array([Ix,Iy,Iz]),np.square(w[i,:]))

gs = gridspec.GridSpec(6, 2)
plt.figure()
plt.subplot(gs[0,:])
plt.plot(w[:,0], 'r-')
plt.plot(w[:,1], 'b-')
plt.plot(w[:,2], 'y-')
plt.plot(totalKE,'k-')
plt.subplot(gs[1,:])
plt.plot(cordvec[:,1,2],cordvec[:,2,2],'k.')
plt.subplot(gs[2:6,:])
plt.plot(cordvec[:,0,2],cordvec[:,1,2],'k.')
plt.show()

##____________________ Plot figure
fig = pylab.figure()
ax = p3.Axes3D(fig)
ax.view_init(elev=20, azim=45) 

def update_line(x, linex,liney,linez):
    for linexi,lineyi,linezi in zip(linex,liney,linez): #must keep the angles part for zip to work and rip the first element
        linexi.set_data([0,cordvec[x,0,0]],[0,cordvec[x,1,0]])
        linexi.set_3d_properties([0,cordvec[x,2,0]])
        lineyi.set_data([0,cordvec[x,0,1]],[0,cordvec[x,1,1]])
        lineyi.set_3d_properties([0,cordvec[x,2,1]])
        linezi.set_data([0,cordvec[x,0,2]],[0,cordvec[x,1,2]])
        linezi.set_3d_properties([0,cordvec[x,2,2]])
#       lineLi.set_data([0,Ltotalinspace[x,0]],[0,Ltotalinspace[x,1]])
#       lineLi.set_3d_properties([0,Ltotalinspace[x,2]])
    return linex,liney,linez

linex=ax.plot([0,cordvec[0,0,0]],[0,cordvec[0,1,0]],[0,cordvec[0,2,0]],'r')
liney=ax.plot([0,cordvec[0,0,1]],[0,cordvec[0,1,1]],[0,cordvec[0,2,1]],'b')
linez=ax.plot([0,cordvec[0,0,2]],[0,cordvec[0,1,2]],[0,cordvec[0,2,2]],'y')
linetrace=ax.plot(cordvec[:,0,2],cordvec[:,1,2],cordvec[:,2,2],'k.')
#lineL=ax.plot([0,Ltotalinspace[i,0]],
#              [0,Ltotalinspace[i,1]],
#              [0,Ltotalinspace[i,2]],'k')

ax.set_xlim3d([-1, 1])
ax.set_ylim3d([-1, 1])
ax.set_zlim3d([0, 1.5])
##ax.axis('equal')
line_ani = animation.FuncAnimation(fig, update_line, list(range(0,N,152)),
                                   fargs=(linex,liney,linez),
                                   interval=100, blit=False,repeat=None)
#line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')
pyplot.show()


