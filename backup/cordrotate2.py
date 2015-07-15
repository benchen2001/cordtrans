import numpy as np
from matplotlib import pyplot
import pylab
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

N = 2000
npzfile = np.load('workfile.npz')
wp = npzfile['wp']
wz = npzfile['wx']      #swap wx to wz
wx = wp[0,:]
wy = wp[1,:]
Iz = npzfile['Ix']
Ix,Iy = npzfile['Iy'],npzfile['Iy']
#print(wx[:30])
#print(wy[:30])
#print(wz[:30])

deltat = npzfile['deltat']
                         #400 steps
def rotation_matrix(rotvec): #rotational formula in euler para or ck para
    #rotvec = np.array([x,y,z])
    amp = np.sqrt(np.dot(rotvec,rotvec))
    axis = rotvec/amp
    phi = amp*deltat
    a = np.cos(phi/2)
    b,c,d = axis*np.sin(phi/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

cordvec = np.zeros((N,3,3))     #(i element, 012xyz components, 012xyz cord vectors)
cordvec[0,0,0] = 1            # initial cord vector
cordvec[0,1,1] = 1  
cordvec[0,2,2] = 1  
Ltotalinbody,Ltotalinspace = np.zeros((N,3)),np.zeros((N,3))

backtoprime = np.array([0,1,0])*np.radians(45)/deltat

k = np.zeros((N,3,3))
k[1,:,:]=np.identity(3)
cordvecprime = cordvec

for i in range(1,N):          # i = 1~N
    #print(i)
    if i >=2:
        rotvecim1 = np.array([wx[i-1],wy[i-1],wz[i-1]])
        #print(np.transpose(rotation_matrix(rotvecim1)))
        k[i,:,:]=np.dot(k[i-1,:,:],
                    np.transpose(rotation_matrix(rotvecim1)))
    rotveci = np.array([wx[i],wy[i],wz[i]])
    for j in range(3):
        cordvec[i,:,j]=np.dot(k[i,:,:],
                              np.dot(rotation_matrix(rotveci),
                                     cordvec[0,:,j]))

        cordvec[i,:,j] = np.dot(np.transpose(rotation_matrix(backtoprime)),
                                cordvec[i,:,j])
    #Ltotalinbody[i,:]= np.array([wx[i],wy[i],wz[i]])
    Ltotalinspace[i,:] = np.dot(k[i,:,:],rotveci)
    Ltotalinspace[i,:] = np.dot(np.transpose(rotation_matrix(backtoprime)),
                      Ltotalinspace[i,:])/50
    
fig = pylab.figure()
ax = p3.Axes3D(fig)
ax.view_init(elev=45, azim=135) 

def update_line(x, linex,liney,linez,lineL):
    for linexi,lineyi,linezi,lineLi in zip(linex,liney,linez,lineL): #must keep the angles part for zip to work and rip the first element
        linexi.set_data([0,cordvec[x,0,0]],[0,cordvec[x,1,0]])
        linexi.set_3d_properties([0,cordvec[x,2,0]])
        lineyi.set_data([0,cordvec[x,0,1]],[0,cordvec[x,1,1]])
        lineyi.set_3d_properties([0,cordvec[x,2,1]])
        linezi.set_data([0,cordvec[x,0,2]],[0,cordvec[x,1,2]])
        linezi.set_3d_properties([0,cordvec[x,2,2]])
        lineLi.set_data([0,Ltotalinspace[x,0]],[0,Ltotalinspace[x,1]])
        lineLi.set_3d_properties([0,Ltotalinspace[x,2]])
    return linex,liney,linez

linex=ax.plot([0,cordvec[0,0,0]],[0,cordvec[0,1,0]],[0,cordvec[0,2,0]],'b')
liney=ax.plot([0,cordvec[0,0,1]],[0,cordvec[0,1,1]],[0,cordvec[0,2,1]],'r')
linez=ax.plot([0,cordvec[0,0,2]],[0,cordvec[0,1,2]],[0,cordvec[0,2,2]],'y')
lineL=ax.plot([0,Ltotalinspace[i,0]],
              [0,Ltotalinspace[i,1]],
              [0,Ltotalinspace[i,2]],'k')

ax.set_xlim3d([-1, 1])
ax.set_ylim3d([-1, 1])
ax.set_zlim3d([0, 1.5]) 

line_ani = animation.FuncAnimation(fig, update_line, list(range(1,N,10)),
                                   fargs=(linex,liney,linez,lineL),
                                   interval=1, blit=False,repeat=None)
#line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')
pyplot.show()

#    ytemp = np.dot(rotation_matrix(cordvec[i,:,0],deltatheta[i,0]),cordvec[i-1,:,1])
#    ztemp = np.dot(rotation_matrix(cordvec[i,:,0],deltatheta[i,0]),cordvec[i-1,:,2])
#    cordvec[i,:,2] = np.dot(rotation_matrix(ytemp,deltatheta[i,1]),ztemp)
#    xtemp = np.dot(rotation_matrix(ytemp,deltatheta[i,1]),cordvec[i-1,:,0])
#    cordvec[i,:,0] = np.dot(rotation_matrix(cordvec[i,:,2],deltatheta[i,2]),xtemp)
#    cordvec[i,:,1] = np.dot(rotation_matrix(cordvec[i,:,2],deltatheta[i,2]),ytemp)
#linex=ax.plot([0,ix[0]],[0,ix[1]],[0,ix[2]],'b')
#lineiy=ax.plot([0,iy[0]],[0,iy[1]],[0,iy[2]],'r')
#lineiiz=ax.plot([0,iiz[0]],[0,iiz[1]],[0,iiz[2]],'y')
