"""Solve symmetric top eular eqs ODE with eular angles and return solutions"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3

def top_der(w,t0,I,Is,ws,tau0):
    """Return the derivatives of eular angles from eular eqs"""
#w(0):phi, w(1):phi_dot, w(2):theta, w(3):theta_dot, w(4):psi
#main program produces w(5):psi_dot
#    print(shape(w))
    ders=[w[1],(Is*ws-2*I*w[1]*np.cos(w[2]))*w[3]/(I*np.sin(w[2])),
          w[3],(tau0-(Is*ws-I*w[1]*np.cos(w[2]))*w[1])*np.sin(w[2])/I,
          ws-w[1]*np.cos(w[2])];
    return ders
'''
I=1.0
Is=1.5
g=9.8
M=1.0
a=0.1                      #inertia , gravity, mass, lever arm
                         #torque value at theta=0
tmax=10
ts=0.05                             #simulation run time and tome interval
'''
def HasbunMatlabTop(I,Is,g,M,a,tau,tmax,numofsteps,
                    phi0,theta0,psi0,
                    phi0d,theta0d,psi0d):
    '''numerically solve the eular-angles Eular equations with python ODE solver odeint'''
    tr=np.linspace(0.0,tmax,numofsteps)
    tau0=M*g*a
#    print(tr)
#    ph0=np.pi
#    th0=0.6*np.pi/2+0.01
 #   ps0=0.0             #init angle values in rad
#ws=1.5+np.sqrt(4*I*tau0*np.cos(th0))/Is  #init spins ang. speed
#    ph0d=0.0 
#    th0d=0.0 #cusps example    #init angular speeds in rad/sec
#w = np.zeros((3000,5))
#ph0d=0.2; th0d=0;%waves (monotonic) example
#ph0d=0.8; th0d=0;%Looping example
    ws = psi0d
    ic1=np.array([phi0,phi0d,theta0,theta0d,psi0])        #initial conditions:
#Use MATLAB's Runge-Kutta (4,5) formula (uncomment/comment as needed)
#opt=odeset('AbsTol',1.e-8,'RelTol',1.e-5);     %user set Tolerances
#[t,w]=ode45('top_der',tr,ic1,opt,I,Is,ws,tau0);%with set tolerance
    w=odeint(top_der,ic1,tr,args=(I,Is,psi0d,tau0))  #default tolerance   
#print(shape(w))  
#w[:,5]=ws-w[:,1]*np.cos(w[:,2])                  #derivative of psi

#Next: plots of the Eulerian angles and their derivatives

#plt.plot(tr,w[:,3])


#plt.show()

#subplot(3,2,1), plot(t,w(:,1),'k'), ylabel('\phi','FontSize',14)
#title(['Top Euler Angles: \phi_0=',num2str(ph0,2),', \theta_0=',...
#        num2str(th0,2),', \psi_0=',num2str(ps0,2),' rad '],'FontSize',11)
#str1=cat(2,'I, I_s=',num2str(I,2),', ',num2str(Is,2),' kgm^2',...
#           ', M=',num2str(M,2),'kg');
#text(0.1,max(w(:,1)*(1-0.1)),str1)
#subplot(3,2,2), plot(t,w(:,2),'b'), ylabel('d\phi/dt','FontSize',14)
#title(['d\phi/dt_0=',num2str(ph0d,2),', d\theta/dt_0=',num2str(th0d,2),...
#       ', \omega_s=',num2str(ws,2),' rad/s'],'FontSize',11)
#str2=cat(2,'g=',num2str(g,2),' m/s^2, a=',num2str(a,2),' m');
#text(0.1,max(w(:,2))*(1+0.1),str2)
#%w(1):phi, w(2):phi_dot, w(3):theta, w(4):theta_dot, w(5):psi, w(6):psi_dot
#subplot(3,2,3), plot(t,w(:,3),'r'), ylabel('\theta','FontSize',14)
#subplot(3,2,4), plot(t,w(:,4),'r'), ylabel('d\theta/dt','FontSize',14)
#subplot(3,2,5), plot(t,w(:,5),'r')
#xlabel('t','FontSize',14),ylabel('\psi','FontSize',14)
#subplot(3,2,6), plot(t,w(:,6),'r')
#xlabel('t','FontSize',14),ylabel('d\psi/dt','FontSize',14)


#Angular momentum components Lx, Ly, Lz, L''z, L''y, L''z
    Lx=I*w[:,3]*np.cos(w[:,0])-(I*w[:,1]*np.cos(w[:,2])-Is*ws)*np.sin(w[:,2])*np.sin(w[:,0])
    Ly=I*w[:,3]*np.sin(w[:,0])+(I*w[:,1]*np.cos(w[:,2])-Is*ws)*np.sin(w[:,2])*np.cos(w[:,0])
    Lz=I*w[:,1]*np.square(np.sin(w[:,2]))+Is*ws*np.cos(w[:,2]);
#Lz" along x, y, z
    Lzppx=Is*ws*np.sin(w[:,2])*np.sin(w[:,0])
    Lzppy=-Is*ws*np.sin(w[:,2])*np.cos(w[:,0])
    Lzppz=Is*ws*np.cos(w[:,2])

#Effective Potential and Nutation angle Turning points 
#'''
#figure
#Vef=M*g*a*cos(w(:,3))+(Lz-Lzppz).^2./(2*I*sin(w(:,3)).^2);
#Ekp=0.5*I*w(:,4).^2; Ep=Ekp+Vef;                %kinetic, and prime energy
#% [thm1,i1]=min(Ekp);[thm2,i2]=min(Ekp(i1+1:N));%aprox Ekp zeros, near ends
#% th1=w(i1,3); th2=w(i2+i1,3);                  %theta values at which Ekp=0
#th1=w(1,3); th2=w(N,3);                         %Ekp zeros-2nd way better
#plot(w(:,3),[Ep,Vef]),xlabel('\theta','FontSize',12)
#ylabel('V_{eff}(\theta), E\prime','FontSize',12)
#str3=cat(2,'\theta_1, \theta_2 = ',num2str(th1,3),', ',...
#         num2str(th2,3),' rad',', E\prime=',num2str(Ep(1),3),' J');
#d1=abs(th2-th1); Vm=min(Vef); d2=Ep(1)-Vm;        %viewing scaling factors
#text(min(w(:,3))*(1+0.15*d1),Ep(1)*(1-0.4*d2),str1)
#text(min(w(:,3))*(1+0.15*d1),Ep(1)*(1-0.25*d2),str2)
#text(min(w(:,3))*(1+0.15*d1),Ep(1)*(1-0.1*d2),str3)%nutation angle points
#title(['Top V_{eff}(\theta): \phi_0=',num2str(ph0,2),', \theta_0=',...
#        num2str(th0,2),', \psi_0=',num2str(ps0,2),' rad, d\phi/dt_0='...
#        ,num2str(ph0d,2),', d\theta/dt_0=',num2str(th0d,2),...
#        ', \omega_s=',num2str(ws,2),' rad/s'],'FontSize',11)
#'''
#
#'''
## ================== Simulation next =================================
#figure
#va1x=min(Lx);va1y=min(Ly);va1z=min(Lz);         %get largest # for window view
#va2x=max(Lx);va2y=max(Ly);va2z=max(Lz);
#vb1x=min(Lzppx);vb1y=min(Lzppy);vb1z=min(Lzppz);
#vb2x=max(Lzppx);vb2y=max(Lzppy);vb2z=max(Lzppz);
#v1x=min(va1x,vb1x);v1y=min(va1y,vb1y);v1z=min(va1z,vb1z);
#v2x=max(va2x,vb2x);v2y=max(va2y,vb2y);v2z=max(va2z,vb2z);
#v=max([abs(v1x),v2x,abs(v1y),v2y,abs(v1z),v2z]);
#axis([-v,v,-v,v,0,v])%window size              
# for i=1:10:N
#     clf
#     axis([-v,v,-v,v,0,v])
#     hold on
#     h(1)=line([0,0],[0,0],[0,Lz(i)],'color','r');%Lz line
#     h(2)=line([0,Lx(i)],[0,Ly(i)],[0,Lz(i)],'color', 'k',...         %L
#         'LineStyle','--','linewidth', 2);
#     h(3)=line([0,0],[0,0],[0,Lzppz(i)],'color','m','Marker','*');    %Lzppz
#     h(4)=line([0,Lzppx(i)],[0,Lzppy(i)],[0,Lzppz(i)],'color', 'b',...%Lzpp
#         'LineStyle','-.','linewidth', 2);
#     box on
#     pause(.05)
# end
#h(5)=plot3(Lx,Ly,Lz,'k:');                        %plot the L trace
#h(6)=plot3(Lzppx,Lzppy,Lzppz,'b:');               %plot the Lz'' trace
#hh=legend(h,'L_z','L','L\prime\prime_z','L\prime\prime','L trace',...
#    'L\prime\prime trace',2); set(hh,'FontSize',11)
#xlabel('x','FontSize',14), ylabel('y','FontSize',14)
#zlabel('z','FontSize',14)
#title(['Top motion: \phi_0=',num2str(ph0,2),', \theta_0=',num2str(th0,2),...
#     ', \psi_0=',num2str(ps0,2),' rad, d\phi/dt_0=',num2str(ph0d,2),...
#     ', d\theta/dt_0=',num2str(th0d,2),', \omega_s=',num2str(ws,2),...
#     ' rad/s'],'FontSize',12)
#text(-v*(1-0.1),v,v*(1-.6),str1)
#text(-v*(1-0.1),v,v*(1-0.5),str2)
#'''
    return Lx,Ly,Lz,Lzppx,Lzppy,Lzppz,w
#
#fig2 = plt.figure(2)
#ax = p3.Axes3D(fig2)
#ax.view_init(elev=20, azim=135)
#plt.plot(Lx,Ly,Lz)
#plt.plot(Lzppx,Lzppy,Lzppz)
#plt.show()