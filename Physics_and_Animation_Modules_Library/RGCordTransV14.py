"""
Module of the Rigid Body class and functions for spinning top's equation of motions. 

Created on Tue Nov 11 16:05:05 2014

@author: whymranderson.blogspot.tw

There are two global functions CK and euler2space.
"""


import numpy as np

import time
#from MATLAB.chapter12 import top
import sys
sys.path.insert(0, 'C:/Documents and Settings/user/My Documents/MATLAB/chapter12/')
import top

tnow = time.time()

print top.__doc__


def CK(rotvec):
    '''Cayley-Klein parameters. Important: the built rotation matrix is with its couterclockwise active sense.
    This is basically the Rodriguez rotation formula in a matric form.
    '''
    amp = np.sqrt(np.dot(rotvec,rotvec))
    if amp == 0:
        ret = np.eye(3)
    else:
        axis = rotvec/amp
        phi = amp % (2*np.pi)
        a = np.cos(phi/2)
        b,c,d = axis*np.sin(phi/2)
        ret =  np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
    return ret
    
def euler2space(eulerW):
    '''Return Euler angles transformation matrix from lab frame to time t frame. 
    
    Matrix takes on passive right-hand sense. The list of input argument w is w(0):phi, w(1):phi_dot, w(2):theta, w(3):theta_dot, w(4):psi
    '''    
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

class RigidBodyObject:
    '''Rigid body class'''
    def __init__(self):
        pass

    # Parameters of top
    # Independent variables    
    M = 0.3      #: mass of top in kg
    R = 0.025  #: top is a disk with radius R in meters
    L = 0.005   #: width of disk in m
    arm = 0.01  #: location of center of mass of top from origin in meters
    counter_weight = 0.05 #: mass doesn't spin along symmetry axis
    counter_weight_location_from_origin = 0.1 #: location from origin
    #Iy= 0.002   #moment of inertial, substitute disk formula if needed 
    Iy=0.25*M*R**2 + 1/12*M*L**2 +M*arm**2 + counter_weight*counter_weight_location_from_origin**2
    '''y moment of inertia 0.25*M*R**2 + 1/12*M*L**2 +M*arm**2 + counter_weight*counter_weight_location_from_origin**2'''    
    #Iz= #0.0008  #
    Iz = 0.5*M*R**2
    '''z moment of inertia 0.5*M*R**2'''    
    Ix = Iy #0.002
    '''x moment of inertia 0.5*M*R**2'''    
    g = 9.8     #: gravity m/s^2
    freq = 20.0 #: top revolution speed in hertz, along symmetric axis
    tn=1.3        #: end of simulation time, take longer time if tn > 10 seconds
    t0=0.0      #: start of time zero
    samplerate = 2000  #: rate of iteration in Hz
    classical_case = 1 #: selection of four typical gyroscopic motions: 1,2,3,4
    # Initial condition
    orien = np.array([-np.pi/3,0,0])            #: starting orientation vector of z' of top from lab xyz
    #orien = np.array([-np.radians(5),0,0])     #starting orientation turn of top from space xyz
    view_azim_angle = 185 #: 3D view azimuthal angle 
    UsePY_ODE = 1    
    use_Jcycle=0
    
    def GenerateDependentVariables(self):
        '''Build and initialize needed arrays and parameters'''
        self.F = np.array([0,0,-self.M*self.g])    #gravity f
        self.N = int(round((self.tn-self.t0)*self.samplerate))    #number of steps
        self.h = (self.tn-self.t0)/self.N               #time step interval
        self.tlist = np.linspace(self.t0,self.tn,self.N+1)
        
        self.w = np.zeros([self.N+1,3])                       # w_s
        self.w_lab = np.zeros([self.N+1,3])                     # w_lab
        self.w[0,2]= 2*np.pi*self.freq                        #wz is the symetry axis
        self.theta0 = self.orien[0]             #euler angle's theta[t0] in rad
        self.b=-(self.Iz-self.Iy)/self.Ix*self.w[0,2]
        self.cordvec = np.zeros((self.N+1,3,3+4+3+3))     
        #(t(i), 012xyz components, 012xyz axes + four rectangle vertices + three hasbun's xyz 789)
        self.cordvec[0,:,2] = np.dot(CK(self.orien),np.array([0,0,1]))  #initial z_s(t0)
        self.cordvec[0,:,1] = np.dot(CK(self.orien),np.array([0,1,0]))  #initial y_s(t0)
        self.cordvec[0,:,0] = np.array([1,0,0])                      #initial x_s(t0)
        self.cordvec[0,:,9] = np.dot(CK(self.orien),np.array([0,0,1]))  #initial z_s(t0) hasbun
        self.cordvec[0,:,8] = np.dot(CK(self.orien),np.array([0,1,0]))  #initial y_s(t0) hasbun
        self.cordvec[0,:,7] = np.array([1,0,0])                      #initial x_s(t0) hasbun
        self.cordvec[0,:,12] = np.dot(CK(self.orien),np.array([0,0,1]))  #initial z_s(t0) hasbun
        self.cordvec[0,:,11] = np.dot(CK(self.orien),np.array([0,1,0]))  #initial y_s(t0) hasbun
        self.cordvec[0,:,10] = np.array([1,0,0])                      #initial x_s(t0) hasbun

        self.Tau_lab_t0 = np.cross(self.arm*self.cordvec[0,:,2],self.F)
        self.TrackMat_t0 = CK(self.orien)
        self.L_b_t0= np.dot(np.array([[self.Ix,0,0],[0,self.Iy,0],[0,0,self.Iz]]),self.w[0,:])
        self.L_lab = np.zeros((self.N+1,3))
        self.L_lab[0,:]=np.dot(self.TrackMat_t0,self.L_b_t0)
        self.Tau_lab_temp = self.Tau_lab_t0           # shallow copy by reference?
#        print(np.dot(self.Tau_lab_temp,self.L_lab[0,:]))
#        print(self.L_lab[0,:])
        self.TrackMat_temp = self.TrackMat_t0
        self.TrackMat_temp_hasbun = self.TrackMat_t0
        self.Tau_body_t0 = np.dot(np.transpose(self.TrackMat_temp),self.Tau_lab_t0)
        self.Tau_body_temp = self.Tau_body_t0
        self.DrawOption = {'A_axes': False,
                           'B_axes': False,
                           'C_axes': False,
                           'A_z_axis_trace': False,
                           'B_z_axis_trace': False,
                           'C_z_axis_trace': False,
                           'A_Angular Momentum Trace' : False,
                           'B_Angular Momentum Trace' : False,
                           'Angular Velocity Trace' : False,                           
                           'A_Angular Momentum Vec' : False,
                           'B_Angular Momentum Vec' : False,
                           'A_Angular Velocity Vec' : False,
                           'A_square': False}
        '''Available drawing options'''
        self.linesarg = ()

    def setcase(self,case,*arg):
        '''Choose from four classical nutation and precession motions or set customary initial condition'''
        self.classical_case = case
        if case == 1: #sharp point
            self.w[0,0] = 0
            self.w[0,1]=  0
        elif case == 2: #round point
            self.w[0,0] = 4
            self.w[0,1]=  0
        elif case == 3: #wave like
            self.w[0,0] = 0
            self.w[0,1]=  -4
        elif case == 4: #circular
            w11,w21,w12,w22 = self.uniformprecesscal()
            possible = 2    
            if possible == 1:
                self.w[0,0] = w11
                self.w[0,1]=  w21
            else:
                self.w[0,0] = w12
                self.w[0,1]=  w22        
        else:
            self.w[0,0] = arg[0]
            self.w[0,1]=  arg[1]
            print('enter Initial body angular velocity')
            print(arg[0],arg[1])
        self.L_b_t0= np.dot(np.array([[self.Ix,0,0],[0,self.Iy,0],[0,0,self.Iz]]),self.w[0,:])
        self.L_lab[0,:]=np.dot(self.TrackMat_t0,self.L_b_t0)

            
    def topEOM(self,wi,tlist,torquei):
        '''EOM of a symmetric top'''
        fderiv = np.dot(np.array([[0,self.b,0],[-self.b,0,0],[0,0,0]]),wi)+np.dot(
                    np.array([[1/self.Ix,0,0],[0,1/self.Iy,0],[0,0,1/self.Iz]]),
                    torquei)
        return fderiv
    def topRK(self,wi,torquei):
        '''4th order Rugge Kutta numeric solver'''
        tnouse=[]
        K1 = self.h*self.topEOM(wi,tnouse,torquei)
        K2 = self.h*self.topEOM(wi+K1/2,tnouse,torquei)
        K3 = self.h*self.topEOM(wi+K2/2,tnouse,torquei)
        K4 = self.h*self.topEOM(wi+K3,tnouse,torquei)
        nextw = wi + (K1+2*K2+2*K3+K4)/6
        return nextw


    #%% Hasbun's code-----------------------
    ##---------Euler angle conversion for Hasbun's code
    ###--------
    phi0 = 0.0
    #theta0 = orien[0] # generated by GenerateDependentVariables method
    psi0 = 0.0
    def HasbunEulerEquationODEsolve(self):
        '''Numerically solve  Newton-Euler equations with Euler angles using Python ode library.
        
        This is rewritten from Prof Hasbun's matlab code in his book **Classical Mechanics with MATLAB Applications**
        This method is essentially the same with the Lagrange method. See text.
        The list of Euler_angles_hasbun = w(0):phi, w(1):phi_dot, w(2):theta, w(3):theta_dot, w(4):psi
        '''        
        eulerMat = np.array([[np.sin(self.theta0)*np.sin(self.psi0), np.cos(self.psi0),0],
                          [np.sin(self.theta0)*np.cos(self.psi0), -np.sin(self.psi0),0],
                           [np.cos(self.theta0),          0 ,            1   ]])   #phi theta shi
        self.phi0d,self.theta0d,self.psi0d = np.linalg.solve(eulerMat,self.w[0,:]) # Calculate initial derivative values
        # ODE solve EOM
        [self.Lx,
         self.Ly,
         self.Lz,
         self.Lzppx,
         self.Lzppy,
         self.Lzppz,
         self.euler_angles_hasbun] = top.HasbunMatlabTop(self.Iy,
                                                         self.Iz,
                                                         self.g,
                                                         self.M,
                                                         self.arm,
                                                         self.M*self.g*self.arm,
                                                         self.tn,
                                                         self.N+1,
                                                         self.phi0,
                                                         self.theta0,
                                                         self.psi0,
                                                         self.phi0d,
                                                         self.theta0d,
                                                         self.psi0d)
        
        # For plotting Hasbun's Lz on unit sphere        
        Lzpp = np.array([self.Lzppx,
                         self.Lzppy,
                         self.Lzppz])
        Lzppt = np.transpose(Lzpp)
        self.Lzppnorm = np.array([row[:]/np.linalg.norm(row) for row in Lzppt])
        L_Bmethod = np.array([self.Lx,
                         self.Ly,
                         self.Lz])
        L_Bt = np.transpose(L_Bmethod)
        self.L_Bnorm = np.array([row[:]/np.linalg.norm(L_Bt[0,:]) for row in L_Bt])        
        for i in range(1,self.N+1):
            for j in range(3):
                self.cordvec[i,:,j+7]=np.dot(np.transpose(euler2space(self.euler_angles_hasbun[i,:])),np.eye(3)[j,:]) # calculate Hasbun's xyz axes in space
#        self.DrawOption['B_axes']=True    
        self.w_body_hasbun = np.zeros([self.N+1,3])
        self.w_body_hasbun = self.eulerW2bodyW(self.euler_angles_hasbun)

    def eulerW2bodyW(self, euler_angles_hasbun):
        '''Convert Hasbun's euler angles solution to body angular velocity
        
        The euler_angles_hasbun array has size N+1-1'''
        mm,nn = np.shape(euler_angles_hasbun)
        w_body_hasbun_temp = np.zeros([mm,3])
        for ii in range(mm):
            w_body_hasbun_temp[ii,0] = euler_angles_hasbun[ii,1]*np.sin(euler_angles_hasbun[ii,2])*np.sin(
                euler_angles_hasbun[ii,4]) + euler_angles_hasbun[ii,3]*np.cos(euler_angles_hasbun[ii,4])
            w_body_hasbun_temp[ii,1] = euler_angles_hasbun[ii,1]*np.sin(euler_angles_hasbun[ii,2])*np.cos(
                euler_angles_hasbun[ii,4]) - euler_angles_hasbun[ii,3]*np.sin(euler_angles_hasbun[ii,4])
            w_body_hasbun_temp[ii,2] = self.w[0,2]# eulerW[0]*np.cos(eulerW[1])
        return w_body_hasbun_temp
    ###-------
    ##--------
    #

    def EulerDCMiter(self): 
        '''Integration of body angular velocity with a DCM method and advances with Newton-Euler equation. 
        
        This is the A method. Also this method use w(t_{i+1})dt in rotation 
        vector approximation by default. One can also choose to use the J-cycle rotation 
        approximation with the use_Jcycle=1 attribute.'''
        for i in range(1,self.N+1):
            #%%'''
            
            if self.UsePY_ODE == 1:  #Use python odeint ODE solver
                from scipy.integrate import odeint            
                ictemp = self.w[i-1,:]
                
                if self.use_Jcycle == 1: #Use J-cycle to better approximate coning rotation vector
                    inserted_tlist = np.linspace(self.tlist[i-1],self.tlist[i],4)
                    self.w[i-1,:],wi0,wi1,self.w[i,:]= odeint(self.topEOM,ictemp,inserted_tlist,(self.Tau_body_temp,))
                    J_dt = inserted_tlist[1]-inserted_tlist[0]               
                    J_rho_i = self.CalculateRhoInJCycle(self.w[i-1,:],wi0,wi1,J_dt)
                    rotation_vec_approximation = J_rho_i
                    
                else:
                    self.w[i-1,:],self.w[i,:]= odeint(self.topEOM,ictemp,self.tlist[i-1:i+1],(self.Tau_body_temp,))
                    rotation_vec_approximation = (self.w[i,:])*self.h
            else:
                self.w[i,:]=self.topRK(self.w[i-1,:],self.Tau_body_temp)
                rotation_vec_approximation = (self.w[i,:])*self.h
            #%%'''
            self.w_lab[i-1,:] = np.dot(self.TrackMat_temp,self.w[i-1,:])
            self.L_lab[i,:]=self.Tau_lab_temp*self.h + self.L_lab[i-1,:]
            #print np.dot(self.cordvec[i-1,:2,2],self.L_lab[i-1,:2])      
            self.TrackMat_temp = np.dot(self.TrackMat_temp,CK(rotation_vec_approximation))        # +w_body_hasbun[i,:] wi here, _body_hasbun
            # TrackMat_temp is the transformation matrix for A_s(ti-1) to A_lab(ti)
            for j in range(3):
                self.cordvec[i,:,j]=np.dot(self.TrackMat_temp,np.eye(3)[j,:])
            self.Tau_lab_temp = np.cross(self.arm*self.cordvec[i,:,2],self.F)             # change from cordvec(i-1) to i
            self.Tau_body_temp = np.dot(np.transpose(self.TrackMat_temp),self.Tau_lab_temp)

        self.L_plot=np.array([row[:]/np.linalg.norm(self.L_lab[0,:]) for row in self.L_lab])# plotting L_lab norm to initial value

        ## construct the shape of a rectangulor representing the top
        for i in range(0,self.N+1):
            for ind,val in enumerate([0,1,0,-1]):
                sign = (-1 if ind > 1 else 1)
                self.cordvec[i,:,ind+3] = ( sign*self.cordvec[i,:,sign*val]/5 + self.cordvec[i,:,2] )/2
#        self.DrawOption['A_axes']=True

    def EulerDCMiter_wt_i(self): 
        ''''This is A method with \\omega(t_{i}) approximation.'''
        for i in range(1,self.N+1):
            self.w[i,:]=self.topRK(self.w[i-1,:],self.Tau_body_temp)
            self.w_lab[i-1,:] = np.dot(self.TrackMat_temp,self.w[i-1,:])
            self.L_lab[i,:]=self.Tau_lab_temp*self.h + self.L_lab[i-1,:]
            #print np.dot(self.cordvec[i-1,:2,2],self.L_lab[i-1,:2])      
            self.TrackMat_temp = np.dot(self.TrackMat_temp,CK((self.w[i-1,:])*self.h))        # +w_body_hasbun[i,:] wi here, _body_hasbun
            # TrackMat_temp is the transformation matrix for A_s(ti-1) to A_lab(ti)
            for j in range(3):
                self.cordvec[i,:,j]=np.dot(self.TrackMat_temp,np.eye(3)[j,:])
            self.Tau_lab_temp = np.cross(self.arm*self.cordvec[i,:,2],self.F)             # change from cordvec(i-1) to i
            self.Tau_body_temp = np.dot(np.transpose(self.TrackMat_temp),self.Tau_lab_temp)

        self.L_plot=np.array([row[:]/np.linalg.norm(self.L_lab[0,:]) for row in self.L_lab])# plotting L_lab norm to initial value

        ## construct the shape of a rectangulor representing the top
        for i in range(0,self.N+1):
            for ind,val in enumerate([0,1,0,-1]):
                sign = (-1 if ind > 1 else 1)
                self.cordvec[i,:,ind+3] = ( sign*self.cordvec[i,:,sign*val]/3 + self.cordvec[i,:,2] )/2

        
    def directDCMiter(self): 
        '''direct DCM iteration from body anguler velocity. This is method C.
        
        w_lab are calculated from method B'''
        C_UseJCycle = 1
        if C_UseJCycle == 1:
            print 'use J-cycle for rotation approximation'
            [self.Lx,
             self.JLy,
             self.JLz,
             self.JLzppx,
             self.JLzppy,
             self.JLzppz,
             self.J_euler_angles_hasbun] = top.HasbunMatlabTop(self.Iy,
                                                             self.Iz,
                                                             self.g,
                                                             self.M,
                                                             self.arm,
                                                             self.M*self.g*self.arm,
                                                             self.tn,
                                                             (self.N +1)*3,
                                                             self.phi0,
                                                             self.theta0,
                                                             self.psi0,
                                                             self.phi0d,
                                                             self.theta0d,
                                                             self.psi0d)
            self.J_w_body_hasbun = np.zeros([(self.N + 1)*3,3])
            self.J_w_body_hasbun = self.eulerW2bodyW(self.J_euler_angles_hasbun)
            print np.shape(self.J_w_body_hasbun)

            for i in range(1,self.N +1):
                C_rotation_vect_i = self.CalculateRhoInJCycle(
                                                         self.J_w_body_hasbun[3*(i-1),:],
                                                         self.J_w_body_hasbun[3*(i-1)+1,:],
                                                         self.J_w_body_hasbun[3*(i-1)+2,:],self.h/3)
                self.TrackMat_temp_hasbun = np.dot(self.TrackMat_temp_hasbun,
                                                   CK(C_rotation_vect_i))        # +w_body_hasbun[i,:] wi here, _body_hasbun
                # TrackMat_temp is the transformation matrix for A_s(ti-1) to A_lab(ti)
                for j in range(3):
                    self.cordvec[i,:,j+10]=np.dot(self.TrackMat_temp_hasbun,np.eye(3)[j,:])   

        
        else:        
            for i in range(1,self.N+1):
                C_rotation_vect_i = self.w_body_hasbun[i-1,:]*self.h
                self.TrackMat_temp_hasbun = np.dot(self.TrackMat_temp_hasbun,CK(C_rotation_vect_i))        # +w_body_hasbun[i,:] wi here, _body_hasbun
                # TrackMat_temp is the transformation matrix for A_s(ti-1) to A_lab(ti)
                for j in range(3):
                    self.cordvec[i,:,j+10]=np.dot(self.TrackMat_temp_hasbun,np.eye(3)[j,:])   
        #self.DrawOption['C_axes'] = True
    
    def IncludeNoiseInOmega(self,NoiseAmplitudeFactor):
        '''Include noise in body angular velocity in C method. 
        
        Return a noisy w_body_hasbun. Take one argument as noise amplitude w.r.t original velocity amplitude.'''
        import random as RandomNumber
        for i in range(0,self.N+1):        
            self.w_body_hasbun[i,:] = np.array([self.w_body_hasbun[i,0]*(1.0)+(NoiseAmplitudeFactor*RandomNumber.uniform(-1,1)),
                                                self.w_body_hasbun[i,1]*(1.0)+(NoiseAmplitudeFactor*RandomNumber.uniform(-1,1)),
                                                self.w_body_hasbun[i,2]*(1.0)+(NoiseAmplitudeFactor*RandomNumber.uniform(-1,1))])
    
    def uniformprecesscal(self):
        '''return the initial condition of body angular velocity required to generate a uniform precession top.
        
        Calculation base on the parameters and the setup of the top. The formula 
        is in Goldstein's book. The formula has a singular point at $\theta$ = 
        90 degree, but a unique way is developed to deal with this. So the program
        will work at any angle.'''
        M,L,g,I3,I1,I2,w3,theta0,psi0 = self.M,self.arm,self.g,self.Iz,self.Ix,self.Iy,self.w[0,2],self.theta0,self.psi0
        a = I1*np.cos(theta0)
        b = -I3*w3
        c = M*g*L
        d = b**2-4*a*c # discriminant
        if np.abs(np.degrees(theta0)) == 90:
            #sys.exit("*** theta0 can not be 90 degrees in uniform precession case ***")
            x = c/-b
            
            wx,wy = x*np.sin(theta0)*np.sin(psi0),x*np.sin(theta0)*np.cos(psi0)
            print 'theta0 = 90 degrees'            
            return wx,wy,wx,wy
            
        elif d < 0:
            print('This equation has no real solution')
        elif d == 0:
            x = (-b+np.sqrt(b**2-4*a*c))/2/a
            print('This equation has one solutions: ' + str(x))
        else:
            x1 = (-b+np.sqrt(b**2-4*a*c))/2/a
            x2 = (-b-np.sqrt(b**2-4*a*c))/2/a
            print('This equation has two solutions: '+ str(x1) + 'and' + str(x2))
    
    #    psi1,psi2=w3-x1*np.cos(theta0),w3-x2*np.cos(theta0)
    
            w11,w12 =x1*np.sin(theta0)*np.sin(psi0),x2*np.sin(theta0)*np.sin(psi0)
            w21,w22 =x1*np.sin(theta0)*np.cos(psi0),x2*np.sin(theta0)*np.cos(psi0)
            print('first w1,w2 = '+str(w11)+ ' and ' +str(w21))
            print('second w1,w2 = '+str(w12)+ ' and ' +str(w22))
            return w11,w21,w12,w22

    def CalculateRhoInJCycle(self,w1,w2,w3,J_dt):
        '''Return rotation vector using J cycle calculation.
        
        Rho should have the same dimension -1 as the sampling rate.'''
        J_omega = np.array([w1,w2,w3])
        for i in range(3):
            J_alpha=np.array([np.sum(J_omega[:i+1,:],axis=0)*J_dt for i in range(3)])
            #for aa in range(N):
            #    print np.cross(alpha[aa,:],omega[aa,:])
            deltaalpha = 0.5*np.sum([np.cross(J_alpha[i,:],J_omega[i,:]) for i in range(3)],axis=0)*J_dt
#            print deltaalpha
            J_rho_i = J_alpha[-1,:] + deltaalpha
        return J_rho_i

        
'''
#%% method of omega sum start
##---------
###--------
for i in range(1,N+1):
    w[i,:]=topRK(w[i-1,:],Tau_body_temp)   
    cordvec[i,:,:3] = np.dot(
    np.transpose(CK(orien + h*np.sum(w[1:i,:],0))),
    np.eye(3))
    Tau_lab_temp = np.cross(arm*cordvec[i,:,2],F)             # change from cordvec(i-1) to i
    Tau_body_temp = np.dot(cordvec[i,:,:3],Tau_lab_temp)
    L_lab[i,:]=Tau_lab_temp*h + L_lab[i-1,:]
###--------
##---------
# method of omega sum end
'''

#plt.figure(5)#,figsize=(5,4))
#plt.plot(w_body_hasbun[:-1,0])
#print(w_body_hasbun[0,1])
#plt.plot()






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
diffindeg = [np.arccos(np.dot(a.Lzppnorm[i,:],
                              a.cordvec[i,:,2]/np.linalg.norm(a.cordvec[i,:,2]))) for i in range(0,a.N,1)]
#print(diffindeg[1000])
plt.plot(np.linspace(0.0,a.tn,a.N/1),np.degrees(diffindeg))
#plt.xlabel('xlabel', fontsize=6)
#plt.ylabel('ylabel', fontsize=6)
#plt.xlabel('time in seconds',fontsize=6)
#plt.ylabel('in degrees',fontsize=6)
plt.title('case='+str(a.classical_case)+',SampleRate='+str(a.samplerate)+'Hz,',fontsize=6)
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


