## Simulation of symetric top with direct intergration of euler equations

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d as p3
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

## Important for CK, quarternion is with its couterclockwise active sense here
## euler para rotational formula in counterclockwise active sense
def CK(rotvec): 
    amp = np.sqrt(np.dot(rotvec,rotvec))
    axis = rotvec/amp
    phi = amp*h % (2*np.pi)
    a = np.cos(phi/2)
    b,c,d = axis*np.sin(phi/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

##_______________________ Parameters of top
M = 1       #kg
R = 0.025   #m
L = 0.01    #m width of disk
arm = 0.04  #m
Iy= 0.002   #0.25*M*R**2 + 1/12*M*L**2 +M*arm**2
Iz= 0.0008  #0.5*M*R**2
Ix= 0.002

g = 9.8                     #m/s^2
F = np.array([0,0,-M*g])    #gravity
freq = 20                   #hertz
N = 1500                     #number of steps
tn=1.6                      #in seconds
t0=0
h = (tn-t0)/N

##_______________________ Initial condition
orien = np.array([-np.pi/4,0,0])    #starting orientation turn of top from space xyz
w = np.zeros([N+1,3])
w[:,2]= w[:,2]+2*np.pi*freq         #wz is the symetry axis
w[0,0] = 10.04597
w[0,1]= -46.25948                           #initial w(t0)
b=-(Iz-Iy)/Ix*w[0,2]
cordvec = np.zeros((N+1,3,3))       #(i element, 012xyz components, 012xyz cord vectors)
cordvec[0,:,2] = np.dot(CK(orien/h),np.array([0,0,1]))  #initial z_s(t0)
cordvec[0,:,1] = np.dot(CK(orien/h),np.array([0,1,0]))  #initial y_s(t0)
cordvec[0,:,0] = np.array([1,0,0])                      #initial x_s(t0)
Tau_s_t0 = np.cross(arm*cordvec[0,:,2],F)
TrackMat_t0 = CK(orien/h)
L_b_t0= np.dot(np.array([[Ix,0,0],[0,Iy,0],[0,0,Iz]]),w[0,:])
L_s = np.zeros((N+1,3))
L_s[0,:]=np.dot(TrackMat_t0,L_b_t0)

##____________________ main calculation
Tau_s_temp = Tau_s_t0           # shallow copy by reference
TrackMat_temp = TrackMat_t0
Tau_body_t0 = np.dot(np.transpose(TrackMat_temp),Tau_s_t0)
Tau_body_temp = Tau_body_t0

for i in range(1,N+1):
    w[i,:]=topRK(w[i-1,:],Tau_body_temp)
    TrackMat_temp = np.dot(TrackMat_temp,CK(w[i,:]))        #wi here, 
    cordvec[i,:,0]=np.dot(TrackMat_temp,np.array([1,0,0]))
    cordvec[i,:,1]=np.dot(TrackMat_temp,np.array([0,1,0]))
    cordvec[i,:,2]=np.dot(TrackMat_temp,np.array([0,0,1]))
    Tau_s_temp = np.cross(arm*cordvec[i,:,2],F)             # change from cordvec(i-1) to i
    Tau_body_temp = np.dot(np.transpose(TrackMat_temp),Tau_s_temp)
    L_s[i,:]=Tau_s_temp*h + L_s[i-1,:]

L_plot=np.array([row[:]/np.linalg.norm(row) for row in L_s])# plotting L_s norm     

##____________________ Plot figure 1
gs = gridspec.GridSpec(4, 8)
plt.figure(figsize=(9, 6))
axw=plt.subplot(gs[:2,:4])
gs.update(left=0.1, bottom=None, right=0.95, top=None, wspace=0.5, hspace=0.5)
plt.plot(w[:,0], 'r-')
plt.plot(w[:,1], 'b-')
axw.set_title('angular velocity wx wy')
axyz=plt.subplot(gs[2:4,:4])
plt.plot(cordvec[:,1,2],cordvec[:,2,2],'k.')
axyz.set_title('position of the tip of z axis in y-z plane')
axxz=plt.subplot(gs[:4,4:8])
plt.plot(cordvec[:,0,2],cordvec[:,1,2],'b.')
axxz.set_title('position of the tip of z axis in x-z plane')
string = 'f='+str(
    freq)+'Hz,'+'g='+str(g)+'m/s^2,'+'m='+str(M)+'kg,'+'arm='+str(
        arm)+'m,'
string2 = '[wx,wy]='+str([w[0,0],w[0,1]])+'rad/s'
string3 = str(N)+'steps in '+str(tn) + ' seconds'
axxz.text(0.02, 0.05, string, fontsize=12, transform = axxz.transAxes)
axxz.text(0.02, 0.1, string2, fontsize=12, transform = axxz.transAxes)
axxz.text(0.02, 0.15, string3, fontsize=12, transform = axxz.transAxes)

plt.show()
plt.close()

## check for the z positions of the top of nutation (when wx=wy=0) consistency
print((np.max(cordvec[:,2,2])-cordvec[0,2,2])/cordvec[0,2,2]*100)

##____________________ Plot figure 2

fig2 = plt.figure(2)
ax = p3.Axes3D(fig2)
ax.view_init(elev=20, azim=135) 

def update_line(x, linex,liney,linez,lineL):
    for linexi,lineyi,linezi,lineLi in zip(linex,liney,linez,lineL): #must keep the angles part for zip to work and rip the first element
        linexi.set_data([0,cordvec[x,0,0]],[0,cordvec[x,1,0]])
        linexi.set_3d_properties([0,cordvec[x,2,0]])
        lineyi.set_data([0,cordvec[x,0,1]],[0,cordvec[x,1,1]])
        lineyi.set_3d_properties([0,cordvec[x,2,1]])
        linezi.set_data([0,cordvec[x,0,2]],[0,cordvec[x,1,2]])
        linezi.set_3d_properties([0,cordvec[x,2,2]])
        lineLi.set_data([0,L_plot[x,0]],[0,L_plot[x,1]])
        lineLi.set_3d_properties([0,L_plot[x,2]])
    return linex,liney,linez

vtx = [cordvec[0,:,0],cordvec[0,:,1],[0,0,0]]
poly = p3.art3d.Poly3DCollection([vtx])
#ax.add_collection3d(poly)

linex=ax.plot([0,cordvec[0,0,0]],[0,cordvec[0,1,0]],[0,cordvec[0,2,0]],'r')
liney=ax.plot([0,cordvec[0,0,1]],[0,cordvec[0,1,1]],[0,cordvec[0,2,1]],'b')
linez=ax.plot([0,cordvec[0,0,2]],[0,cordvec[0,1,2]],[0,cordvec[0,2,2]],'y',linewidth=2.0)
linetrace=ax.plot(cordvec[:,0,2],cordvec[:,1,2],cordvec[:,2,2],'g.',markersize=2)
lineLtrance=ax.plot(L_plot[:,0],L_plot[:,1],L_plot[:,2],'k.',markersize=2)
lineL=ax.plot([0,L_plot[0,0]],
              [0,L_plot[0,1]],
              [0,L_plot[0,2]],'k-')

ax.set_xlim3d([-1, 1])
ax.set_ylim3d([-1, 1])
ax.set_zlim3d([0, 1.5])
line_ani = animation.FuncAnimation(fig2, update_line, list(range(0,N,round(5))),
                                   fargs=(linex,liney,linez,lineL),
                                   interval=10, blit=False,repeat=None)
#line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')
plt.show()


