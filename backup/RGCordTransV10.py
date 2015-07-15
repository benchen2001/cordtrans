## Simulation of symetric top with direct intergration of euler equations

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import time
#from MATLAB.chapter12 import top
import sys
sys.path.insert(0, 'C:/Documents and Settings/user/My Documents/MATLAB/chapter12/')
import top

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
M = 1       #mass in kg
#R = 0.025  #radius of disk in m
#L = 0.01   #width of disk in m
arm = 0.04  #level arm of Center of mass in m
Iy= 0.002   #moment of inertial, substitute disk formula if needed Iy=0.25*M*R**2 + 1/12*M*L**2 +M*arm**2
Iz= 0.0008  #Iz = 0.5*M*R**2
Ix= 0.002

g = 9.8                     #gravity m/s^2
F = np.array([0,0,-M*g])    #gravity f
freq = 15                 #top revolution speed in hertz, along symmetric axis
tn=1                      # end of simulation time, take longer time if tn > 10 seconds
t0=0.0                      # start of time zero
samplerate = 2000           # rate of iteration in Hz
N = int(round((tn-t0)*samplerate))    #number of steps
h = (tn-t0)/N               #time step interval

#%%_________________________________________ Initial condition
orien = np.array([-np.pi/3,0,0])            #starting orientation vector of top from lab xyz
#orien = np.array([-np.radians(5),0,0])     #starting orientation turn of top from space xyz
w = np.zeros([N+1,3])                       # w_s
w_b = np.zeros([N+1,3])                     # w_lab
w[0,2]= 2*np.pi*freq                        #wz is the symetry axis



#%%initial w(t0), selection of four typical gyroscopic motions: 1,2,3,4
case = 3
if case == 1: #sharp point
    w[0,0] = 0
    w[0,1]=  0
elif case == 2: #round point
    w[0,0] = 4
    w[0,1]=  0
elif case == 3: #wave like , b/a > u2, f=20, w[0,1]=-4.5
    w[0,0] = 0
    w[0,1]=  -5
elif case == 4: #circular , w[0,0]=0,w[0,1]=-5,f=15
    print('check theta0 sign')
    import uniformprecess
    w11,w21,w12,w22 = uniformprecess.uniformprecesscal(M,arm,g,Iz,Ix,Iy,
                                                       w[0,2],-np.abs(np.degrees(orien[0])))
    possible = 2    
    if possible == 1:
        w[0,0] = w11
        w[0,1]=  w21
    else:
        w[0,0] = w22
        w[0,1]=  w12        
else:
    w[0,0] = 0
    w[0,1]=  0
    print('enter Initial body angular velocity')
    
                        
b=-(Iz-Iy)/Ix*w[0,2]
cordvec = np.zeros((N+1,3,3+4+3))     
#(t(i), 012xyz components, 012xyz axes + four rectangle vertices + three hasbun's xyz 789)
cordvec[0,:,2] = np.dot(CK(orien),np.array([0,0,1]))  #initial z_s(t0)
cordvec[0,:,1] = np.dot(CK(orien),np.array([0,1,0]))  #initial y_s(t0)
cordvec[0,:,0] = np.array([1,0,0])                      #initial x_s(t0)
cordvec[0,:,9] = np.dot(CK(orien),np.array([0,0,1]))  #initial z_s(t0) hasbun
cordvec[0,:,8] = np.dot(CK(orien),np.array([0,1,0]))  #initial y_s(t0) hasbun
cordvec[0,:,7] = np.array([1,0,0])                      #initial x_s(t0) hasbun

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



#%% Hasbun's code-----------------------
##---------Euler angle conversion for Hasbun's code
###--------
phi0 = 0.0
theta0 = orien[0]
psi0 = 0.0
eulerMat = np.array([[np.sin(theta0)*np.sin(psi0), np.cos(psi0),0],
                     [np.sin(theta0)*np.cos(psi0), -np.sin(psi0),0],
                     [np.cos(theta0),          0 ,            1   ]])   #phi theta shi
phi0d,theta0d,psi0d = np.linalg.solve(eulerMat,w[0,:])
# Calculate initial derivative values
print(phi0d,theta0d,psi0d,w[0,:])

###EOM
Lx,Ly,Lz,Lzppx,Lzppy,Lzppz,euler_angles_hasbun = top.HasbunMatlabTop(Iy,Iz,g,M,arm,M*g*arm,tn,h,
                                                 phi0,theta0,psi0,
                                                 phi0d,theta0d,psi0d)
                                                 
#w(0):phi, w(1):phi_dot, w(2):theta, w(3):theta_dot, w(4):psi

def eulerW2bodyW(eulerW):
    wbx = eulerW[1]*np.sin(eulerW[2])*np.sin(
            eulerW[4]) + eulerW[3]*np.cos(eulerW[4])
    wby = eulerW[1]*np.sin(eulerW[2])*np.cos(
            eulerW[4]) - eulerW[3]*np.sin(eulerW[4])
    wbz = w[0,2]#eulerW[1]*np.cos(eulerW[2]) + eulerW[4]# eulerW[0]*np.cos(eulerW[1])
    return wbx,wby,wbz
def euler2space(eulerW):
    # passive right-hand matrix    #w(0):phi, w(1):phi_dot, w(2):theta, w(3):theta_dot, w(4):psi
    euler2spaceMatrix = np.array([[np.cos(eulerW[4])*np.cos(eulerW[0])-np.cos(eulerW[2])*np.sin(eulerW[0])*np.sin(eulerW[4]),
                                   np.cos(eulerW[4])*np.sin(eulerW[0])+np.cos(eulerW[2])*np.cos(eulerW[0])*np.sin(eulerW[4]),
                                   np.sin(eulerW[4])*np.sin(eulerW[2])],
                                  [-np.sin(eulerW[4])*np.cos(eulerW[0])-np.cos(eulerW[2])*np.sin(eulerW[0])*np.cos(eulerW[4]),
                                   -np.sin(eulerW[4])*np.sin(eulerW[0])+np.cos(eulerW[2])*np.cos(eulerW[0])*np.cos(eulerW[4]),
                                   np.cos(eulerW[4])*np.sin(eulerW[2])],
                                  [np.sin(eulerW[2])*np.sin(eulerW[0]),          
                                   -np.sin(eulerW[2])*np.cos(eulerW[0]),
                                   np.cos(eulerW[2])]
                                 ])
    return euler2spaceMatrix

mm,nn = shape(euler_angles_hasbun)
# euler_angles_hasbun has size N+1 -1
#print(mm)
w_body_hasbun = np.zeros([N+1,3])

for ii in range(mm):
    w_body_hasbun[ii,:] = eulerW2bodyW(euler_angles_hasbun[ii,:])

#print(w_body_hasbun[1000,:])
#print(w_body_hasbun[N-1,1])
Lzpp = np.array([Lzppx,
                 Lzppy,
                 Lzppz])
Lzppt = np.transpose(Lzpp)
#print(shape(Lzppt))
#print(shape(cordvec))
Lzppnorm = np.array([row[:]/np.linalg.norm(row) for row in Lzppt])
###-------
##--------
#




#%% main DCM interation
##------------
###-----------
for i in range(1,N+1):
    w[i,:]=topRK(w[i-1,:],Tau_body_temp)
#    print(TrackMat_temp)
    TrackMat_temp = np.dot(TrackMat_temp,CK((w[i,:])*h))        # +w_body_hasbun[i,:] wi here, _body_hasbun
    # TrackMat_temp is the transformation matrix for A_s(ti-1) to A_lab(ti)
#    print(TrackMat_temp)
    for j in range(3):
       cordvec[i,:,j]=np.dot(TrackMat_temp,np.eye(3)[j,:])
    w_b[i-1,:] = np.dot(TrackMat_temp,w[i,:])
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
    
###-----------
##------------
#


'''
#%% method of omega sum start
##---------
###--------
for i in range(1,N+1):
    w[i,:]=topRK(w[i-1,:],Tau_body_temp)   
    cordvec[i,:,:3] = np.dot(
    np.transpose(CK(orien + h*np.sum(w[1:i,:],0))),
    np.eye(3))
    Tau_s_temp = np.cross(arm*cordvec[i,:,2],F)             # change from cordvec(i-1) to i
    Tau_body_temp = np.dot(cordvec[i,:,:3],Tau_s_temp)
    L_s[i,:]=Tau_s_temp*h + L_s[i-1,:]
###--------
##---------
# method of omega sum end
'''

plt.figure(5)#,figsize=(5,4))
plt.plot(w_body_hasbun[:,0]-w[:,0])
#print(w_body_hasbun[0,1])
plt.plot()

L_plot=np.array([row[:]/np.linalg.norm(row) for row in L_s])# plotting L_s norm





'''
#%%
##____________________ Plot figure 1
###
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
###
##
#'''

'''
#%%4down Plot of differences between my code and Hasbun's in solid angle in degrees
##------------
###-----------
plt.figure(4,figsize=(2,2))
diffindeg = [np.arccos(np.dot(Lzppnorm[i,:],
                              cordvec[i,:,2]/np.linalg.norm(cordvec[i,:,2]))) for i in range(0,N,1)]
#print(diffindeg[1000])
plt.plot(np.linspace(0.0,tn,N/1),np.degrees(diffindeg))
#plt.xlabel('xlabel', fontsize=6)
#plt.ylabel('ylabel', fontsize=6)
#plt.xlabel('time in seconds',fontsize=6)
#plt.ylabel('in degrees',fontsize=6)
plt.title('case='+str(case)+',SampleRate='+str(samplerate)+'Hz,',fontsize=6)
#title(r'angles between two symmetric axes in Hasbun and my code', fontsize=6)
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\compare_smallfig_case2.pgf')
###-----------
##------------
#4up
'''

'''
#%% 6 down Plot the DCM xyz axes and Hasbun xyz axes for comparison
##------------
###-----------
fig3=plt.figure(6)#,figsize=(2,2))
ax3 = p3.Axes3D(fig3)
ax3.view_init(elev=20, azim=135)

#DCM xyz
baxes = [ax3.plot(*[[0,cordvec[1000,j,i]] for j in range(3) ]) for i in range(3)]
plt.setp(baxes[2],marker='o') #set z axis marker

#Hasbun
baxes_hasbun = [ax3.plot(*[[0,cordvec[1000,j,i]] for j in range(3) ]) for i in range(7,10)]
plt.setp(baxes_hasbun,lw=3) #set z axis marker
plt.setp(baxes_hasbun[2],marker='o') #set z axis marker


plt.show()
#plt.xlabel('xlabel', fontsize=6)
#plt.ylabel('ylabel', fontsize=6)
#plt.xlabel('time in seconds',fontsize=6)
#plt.ylabel('in degrees',fontsize=6)
#title(r'angles between two symmetric axes in Hasbun and my code', fontsize=6)
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\compare_smallfig_case2.pgf')
###-----------
##------------
# 6 up
'''


#%%____________________ Plot figure top animation
##
###
fig2 = plt.figure(2)#,figsize=(6,4))
ax = p3.Axes3D(fig2)
ax.view_init(elev=20, azim=135)

## shape of a rectangulor top
for i in range(0,N+1):
    for ind,val in enumerate([0,1,0,-1]):
        sign = (-1 if ind > 1 else 1)
        cordvec[i,:,ind+3] = ( sign*cordvec[i,:,sign*val]/3 + cordvec[i,:,2] )/2
		
def update_line(x, blines,lineL,lineW):
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
        if 3<=ind<7:           # update rectangle
            for blinei,lineLi in zip(bline,lineL):
                blinei.set_data([cordvec[x,0,polynum[ind-3]],cordvec[x,0,polynum[ind+1-3]]]
                                ,[cordvec[x,1,polynum[ind-3]],cordvec[x,1,polynum[ind+1-3]]])
                blinei.set_3d_properties([cordvec[x,2,polynum[ind-3]]
                                          ,cordvec[x,2,polynum[ind+1-3]]])	
        if ind>=7:     # update hasbun's body axes    
            for blinei,lineLi in zip(bline,lineL): 
                if ind == 9:                
                    blinei.set_data([0,(cordvec[x,0,ind]+cordvec[x,0,9])/2],
                                 [0,(cordvec[x,1,ind]+cordvec[x,1,9])/2])                
                    blinei.set_3d_properties([0,(cordvec[x,2,ind]+cordvec[x,2,9])/2])
                else:
                    blinei.set_data(
                    [0+cordvec[x,0,9]/2 , (cordvec[x,0,ind]/2+cordvec[x,0,9])/2],
                    [0+cordvec[x,1,9]/2 , (cordvec[x,1,ind]/2+cordvec[x,1,9])/2]) 
                    blinei.set_3d_properties(
                    [0+cordvec[x,2,9]/2 , (cordvec[x,2,ind]/2+cordvec[x,2,9])/2])
    #plotting angular velocity vector norm to w[0,2]
    for lineWi,lineLi in zip(lineW,lineL):
        lineWi.set_data([0,w_b_norm[x,0]],[0,w_b_norm[x,1]])
        lineWi.set_3d_properties([0,w_b_norm[x,2]])
    return blines,lineL,lineWi

# plot the initial body xyz axes DCM method
baxes = [ax.plot(*[ [0+cordvec[0,j,2]/2,(cordvec[0,j,i]+cordvec[0,j,2])/2] for j in range(3) ]) for i in range(3)]
plt.setp(baxes[2],marker='o') #set z axis marker
## unpacking list to tuple and list comprehension, see below
## baxes[1]=ax.plot([0,cordvec[0,0,0]],[0,cordvec[0,1,0]],[0,cordvec[0,2,0]]), plot x axis line

# plot the initial body xyz axes from Hasbun's euler angle method
baxes_hasbun = [ax.plot(*[ [0+cordvec[0,j,9]/2,(cordvec[0,j,i]+cordvec[0,j,9])/2] for j in range(3) ]) for i in range(7,10)]
plt.setp(baxes_hasbun,linestyle=':') #set z axis marker
plt.setp(baxes_hasbun,lw=5)

polynum = [3,4,5,6,3]   # indexing vertices of a closed rectangle

# Color of the rectangle
##
import matplotlib.cm as mplcm
import matplotlib.colors as colors
cm = plt.get_cmap('Oranges')
cNorm  = colors.Normalize(vmin=-1, vmax=3)
scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(4)])
##
#

# plot the initial square box
polylines = [ ax.plot(*[ [cordvec[0,j,polynum[i]],
                          cordvec[0,j,polynum[i+1]]] for j in range(3) ]
                      ) for i in range(4) ]
plt.setp(polylines,lw=5)

#change previous w(t_i-1) for w(t i-1) comparison
linetrace=ax.plot(cordvec[:,0,2],cordvec[:,1,2],cordvec[:,2,2],'b-',
                  markersize=1,label='DCM with $\omega(t_i)$')
#
#lineLtrance=ax.plot(L_plot[:,0],L_plot[:,1],L_plot[:,2],'k.',markersize=2)
HasbunLtrance=ax.plot(Lzppnorm[:,0],Lzppnorm[:,1],Lzppnorm[:,2],'k-',
                      markersize=2,label='$\omega(t_{i+1})$ and Hasbun\'s result')
#ax.legend()

lineL=ax.plot([0,L_plot[0,0]],
              [0,L_plot[0,1]],
              [0,L_plot[0,2]],'k-')

#HasbunL=ax.plot([0,L_plot[0,0]],
#               [0,L_plot[0,1]],
#               [0,L_plot[0,2]],'k-')

# plot angular velocity vector normalized to w(t0)--
w_b_norm = w_b/w[0,2]
#print(w[0,2],np.linalg.norm(w_b[0,:]))
lineW=ax.plot([0,w_b_norm[0,0]],
              [0,w_b_norm[0,1]],
              [0,w_b_norm[0,2]],'g-')
# ----



ax.set_xlim3d([-0.8, 0.8])
ax.set_ylim3d([-0.8, 0.8])
ax.set_zlim3d([0, 0.8])
# Keep "line_ani =" part, otherwise won't work, timedanimation need an return to clear the figure
line_ani = animation.FuncAnimation(fig2, update_line, list(range(1,N,5)),
                                  fargs=(baxes+polylines+baxes_hasbun,lineL,lineW),
                                  interval=10, blit=False,repeat=False)
#line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')
plt.show()
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\wti_wtiplus1.pgf')
###
##
#
