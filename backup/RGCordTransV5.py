## Simulation of symetric top with direct intergration of euler equations

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import time

tnow = time.time()

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

## Important for CK, quarternion is with its couterclockwise active sense here
## euler para rotational formula in counterclockwise active sense
def CK(rotvec): 
    amp = np.sqrt(np.dot(rotvec,rotvec))
    axis = rotvec/amp
#    print(amp)
    phi = amp % (2*np.pi)
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
freq = 17                   #hertz
N = 2000*2                    #number of steps
tn=1.2*2                        #in seconds
t0=0
h = (tn-t0)/N

##_______________________ Initial condition
orien = np.array([-.9525,0,0])    #starting orientation turn of top from space xyz
w = np.zeros([N+1,3])
w[0,2]= 2*np.pi*freq         #wz is the symetry axis
w[0,0] = 0
w[0,1]=  0                           #initial w(t0)
b=-(Iz-Iy)/Ix*w[0,2]
cordvec = np.zeros((N+1,3,3+4))     #(t(i), 012xyz components, 012xyz axes + vertices)
cordvec[0,:,2] = np.dot(CK(orien),np.array([0,0,1]))  #initial z_s(t0)
cordvec[0,:,1] = np.dot(CK(orien),np.array([0,1,0]))  #initial y_s(t0)
cordvec[0,:,0] = np.array([1,0,0])                      #initial x_s(t0)
Tau_s_t0 = np.cross(arm*cordvec[0,:,2],F)
TrackMat_t0 = CK(orien)
L_b_t0= np.dot(np.array([[Ix,0,0],[0,Iy,0],[0,0,Iz]]),w[0,:])
L_s = np.zeros((N+1,3))
L_s[0,:]=np.dot(TrackMat_t0,L_b_t0)

##____________________ main calculation
Tau_s_temp = Tau_s_t0           # shallow copy by reference
TrackMat_temp = TrackMat_t0
Tau_body_t0 = np.dot(np.transpose(TrackMat_temp),Tau_s_t0)
Tau_body_temp = Tau_body_t0


#method of trackMat start
for i in range(1,N+1):
    w[i,:]=topRK(w[i-1,:],Tau_body_temp)
#    print(TrackMat_temp)
    TrackMat_temp = np.dot(TrackMat_temp,CK(w[i,:]*h))        #wi here, 
#    print(TrackMat_temp)
    for j in range(3):
        cordvec[i,:,j]=np.dot(TrackMat_temp,np.eye(3)[j,:])
# if rotation vector is small enough, the following should work
#        transmatrix = np.dot(CK(orien),CK(h*np.sum(w[1:i+1,:],0)))      
#        cordvec[i,:,j] = np.dot(
#        np.transpose(transmatrix),
#        np.eye(3)[j,:])
#    print(CK(orien + h*np.sum(w[1:i+1,:])))
    Tau_s_temp = np.cross(arm*cordvec[i,:,2],F)             # change from cordvec(i-1) to i
    Tau_body_temp = np.dot(np.transpose(TrackMat_temp),Tau_s_temp)
    L_s[i,:]=Tau_s_temp*h + L_s[i-1,:]
#    time.sleep(5)
# method of trackMat end

'''# method of omega sum start
for i in range(1,N+1):
    w[i,:]=topRK(w[i-1,:],Tau_body_temp)   
    cordvec[i,:,:3] = np.dot(
    np.transpose(CK(orien + h*np.sum(w[1:i,:],0))),
    np.eye(3))
    Tau_s_temp = np.cross(arm*cordvec[i,:,2],F)             # change from cordvec(i-1) to i
    Tau_body_temp = np.dot(cordvec[i,:,:3],Tau_s_temp)
    L_s[i,:]=Tau_s_temp*h + L_s[i-1,:]
'''# method of omega sum end


L_plot=np.array([row[:]/np.linalg.norm(row) for row in L_s])# plotting L_s norm


'''##____________________ Plot figure 1
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
plt.pause(2)
plt.close()
print(time.time()-tnow)
## check for the z positions of the top of nutation (when wx=wy=0) consistency
#print((np.min(cordvec[:,2,2])-cordvec[0,2,2])/cordvec[0,2,2]*100)
'''

##____________________ Plot figure 2

fig2 = plt.figure(2,figsize=(10,5))
ax = p3.Axes3D(fig2)
ax.view_init(elev=20, azim=135)

## shape of a rectangulor top
for i in range(0,N+1):
    for ind,val in enumerate([0,1,0,-1]):
        sign = (-1 if ind > 1 else 1)
        cordvec[i,:,ind+3] = ( sign*cordvec[i,:,sign*val]/3 + cordvec[i,:,2] )/2
		
def update_line(x, blines,lineL):
    for ind,bline in enumerate(blines):
        if ind < 3:	# update baxes body axes
            for blinei,lineLi in zip(bline,lineL): 
                if ind == 2:                
                    blinei.set_data([0,(cordvec[x,0,ind]+cordvec[x,0,2])/2],
                                 [0,(cordvec[x,1,ind]+cordvec[x,1,2])/2])                
                    blinei.set_3d_properties([0,(cordvec[x,2,ind]+cordvec[x,2,2])/2])
                else:
                    blinei.set_data(
                    [0+cordvec[x,0,2]/2 , (cordvec[x,0,ind]/2+cordvec[x,0,2])/2],
                    [0+cordvec[x,1,2]/2 , (cordvec[x,1,ind]/2+cordvec[x,1,2])/2]) 
                    blinei.set_3d_properties(
                    [0+cordvec[x,2,2]/2 , (cordvec[x,2,ind]/2+cordvec[x,2,2])/2])
##            baxisi.set_data([0,cordvec[x,0,0]],[0,cordvec[x,1,0]])
##            baxisi.set_3d_properties([0,cordvec[x,2,0]])
                lineLi.set_data([0,L_plot[x,0]],[0,L_plot[x,1]])
                lineLi.set_3d_properties([0,L_plot[x,2]])
        else:           # update rectangle
            for blinei,lineLi in zip(bline,lineL):
                blinei.set_data([cordvec[x,0,polynum[ind-3]],cordvec[x,0,polynum[ind+1-3]]]
                                ,[cordvec[x,1,polynum[ind-3]],cordvec[x,1,polynum[ind+1-3]]])
                blinei.set_3d_properties([cordvec[x,2,polynum[ind-3]]
                                          ,cordvec[x,2,polynum[ind+1-3]]])	
    return blines,lineL


baxes = [ax.plot(*[ [0+cordvec[0,j,2]/2,(cordvec[0,j,i]+cordvec[0,j,2])/2] for j in range(3) ]) for i in range(3)]
## unpacking list to tuple and list comprehension, see below
## baxes[1]=ax.plot([0,cordvec[0,0,0]],[0,cordvec[0,1,0]],[0,cordvec[0,2,0]]), plot x axis line
polynum = [3,4,5,6,3]   # indexing vertices of a closed rectangle
polylines = [ ax.plot(*[ [cordvec[0,j,polynum[i]],
                          cordvec[0,j,polynum[i+1]]] for j in range(3) ]
                      ) for i in range(4) ]
setp(polylines,lw=5)
linetrace=ax.plot(cordvec[:,0,2],cordvec[:,1,2],cordvec[:,2,2],'b.',markersize=1)
#lineLtrance=ax.plot(L_plot[:,0],L_plot[:,1],L_plot[:,2],'k.',markersize=2)
lineL=ax.plot([0,L_plot[0,0]],
              [0,L_plot[0,1]],
              [0,L_plot[0,2]],'k-')

ax.set_xlim3d([-0.8, 0.8])
ax.set_ylim3d([-0.8, 0.8])
ax.set_zlim3d([0, 0.8])
# Keep "line_ani =" part, otherwise won't work, timedanimation need an return to clear the figure
line_ani = animation.FuncAnimation(fig2, update_line, list(range(0,N,3)),
                                   fargs=(baxes+polylines,lineL),
                                   interval=1, blit=False,repeat=True)
#line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')
plt.show()

f = open('workfilex.txt', 'w')
cordvec[:,0,2].tofile(f,sep=" ",format="%s")
f.close()
f1 = open('workfiley.txt', 'w')
cordvec[:,1,2].tofile(f1,sep=" ",format="%s")
f1.close()
f2 = open('workfilez.txt', 'w')
cordvec[:,2,2].tofile(f2,sep=" ",format="%s")
f2.close()
