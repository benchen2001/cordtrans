import numpy as np
from matplotlib import pyplot
import pylab
import mpl_toolkits.mplot3d.axes3d as p3

def rotation_matrix(axis,theta): #theta in radian
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    b,c,d = -axis*np.sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

x = np.array([1,0,0])
y = np.array([0,1,0])
z = np.array([0,0,1])

theta = 10/180*np.pi
firstrot= np.transpose(rotation_matrix(x,theta))

iy = np.dot(firstrot,y)
iz = np.dot(firstrot,z)

secrot= np.transpose(rotation_matrix(iy,theta))

ix = np.dot(secrot,x)
iiz = np.dot(secrot,iz)

fig = pylab.figure()
ax = p3.Axes3D(fig)
ax.view_init(elev=15, azim=45) 

linex=ax.plot([0,1],[0,0],[0,0],'b')
liney=ax.plot([0,iy[0]],[0,iy[1]],[0,iy[2]],'r')
linez=ax.plot([0,iz[0]],[0,iz[1]],[0,iz[2]],'y')

linex=ax.plot([0,ix[0]],[0,ix[1]],[0,ix[2]],'b')
lineiy=ax.plot([0,iy[0]],[0,iy[1]],[0,iy[2]],'r')
lineiiz=ax.plot([0,iiz[0]],[0,iiz[1]],[0,iiz[2]],'y')

ax.set_xlim3d([0, 1])
ax.set_ylim3d([0.0, 1])
ax.set_zlim3d([0.0, 1]) 

pyplot.show()
