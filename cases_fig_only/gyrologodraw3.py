import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3

a0 = 1.0
xx = np.linspace(0.001 ,25, 50)

def hydrogenN3l2(x):
    f = 1.0/81/np.sqrt(6*np.pi)/np.power(a0,1.5)*x*x/a0/a0*np.exp(-x/3/a0)   
    return f
def gyroshape():
    y = 200*np.pi*xx*xx*hydrogenN3l2(xx)*hydrogenN3l2(xx)
    vec = np.array([np.zeros(50),xx,y])
    vecm = np.array([np.zeros(50),xx[::-1],-y[::-1]])
    vect = np.concatenate((vec,vecm),axis=1)
  #  print(np.shape(vect))
    return vect
def rotation_matrix(axis,theta):
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    b,c,d = axis*np.sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
def circle_arc(axis,start_v,end_v,num_points):
    #start_v,end_v have to be perpendicular to zxis
    axis = axis/np.linalg.norm(axis)
    start_v = start_v/np.linalg.norm(start_v)
   # print(start_v)
    end_v = end_v/np.linalg.norm(end_v)
    theta = np.arccos(np.dot(start_v,end_v))
#    print(np.dot(start_v,end_v))
 #   print(theta)
    theta_s = list(np.arange(0.0, theta + theta/num_points, theta/num_points))
    circle_vecs = np.zeros([len(theta_s),3])
    for i,thetai in enumerate(theta_s):
        makecir = rotation_matrix(axis,thetai)
        circle_vecs[i,:] = np.dot(makecir,start_v)
    return circle_vecs

gyrovec = gyroshape()
gyrovec2 = gyroshape()
mmg,nng = np.shape(gyrovec)
for nni in range(nng):
    gyrovec[:,nni] = np.dot(rotation_matrix([-1,0,0],np.radians(100)),
                            gyrovec[:,nni])
    gyrovec2[:,nni] = np.dot(rotation_matrix([0,-1.0*np.sin(np.radians(10)),
                            -1.0*np.cos(np.radians(10))],np.radians(100)),
                            gyrovec[:,nni])




#fig2=plt.figure(2)
#plt.text(0,0,'\n Technology',fontsize=28)
#plt.plot(xx,200*np.pi*xx*xx*hydrogenN3l2(xx)*hydrogenN3l2(xx))

fig3 = plt.figure(3,figsize=(1.5,2),dpi=200)
ax3 = p3.Axes3D(fig3)
ax3.view_init(elev=5, azim=2)
#ax3.set_color_cycle('b')

#print(fig3.get_dpi())
#print(fig3.get_size_inches())

markersizearray =  2600*np.pi*xx*xx*hydrogenN3l2(xx)*hydrogenN3l2(xx)

grx = np.array([8.5,12.0])
grr = 200*np.pi*grx*grx*hydrogenN3l2(grx)*hydrogenN3l2(grx)
gyroZvec = np.dot(rotation_matrix((-1,0,0),np.radians(100)),[0,1,0])
arc2unrotated = grr[0]*circle_arc([0,1,0],[0,0,grr[0]],[0.001,0,-grr[0]],60) + grx[0]*np.array([0,1,0])
print np.shape(np.transpose(arc2unrotated))
arc2 = np.dot(rotation_matrix((-1,0,0),np.radians(100)),np.transpose(arc2unrotated)) 
print np.shape(arc2)
ax3.scatter(arc2[0,:],arc2[1,:],arc2[2,:],alpha=0.2,s = markersizearray[::-1]*0.65)
arc3 = grr[1]*circle_arc([0,1,np.sqrt(3)],[0,-np.sqrt(3),1],[0.01,np.sqrt(3),-1],20) + grx[1]*np.array([0,1,np.sqrt(3)])/2
#larc3, = ax3.plot(arc3[:,0],arc3[:,1],arc3[:,2],'k',lw=4)
#ax3.text(0,15,17,r'',fontsize=28,verticalalignment='top')

#ax3.plot(*gyrovec2[:,:50],linewidth=7,color='k',alpha=0.4)
#ax3.plot(*gyrovec,linewidth=6,color='k')

#plt.figure(4)
#plt.plot(markersizearray)
ax3.scatter(*gyrovec2[:,:50],s=markersizearray*0.8,color='r',alpha=0.2)
ax3.scatter(*gyrovec,s=markersizearray,color='k')

ax3.set_axis_off()  #-> this can turn off the background curtain

ax3.tick_params(labelbottom='off', labeltop='off', labelleft='off', labelright='off')
ax3.set_xlim(-6,3)
ax3.set_ylim(-6,3)
ax3.set_zlim(-20,-3)

print(fig3.get_size_inches())

plt.show()
plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\gyrologo3.png')
