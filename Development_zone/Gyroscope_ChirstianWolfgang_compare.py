import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunctionV3 as RBPlot
import RGCordTransV13 as RG

b = RG.RigidBodyObject()
#set para here
b.tn=10.2

# Christian's setup
b.g = 1.0
size = 2.0
b.M = 1.0
z0 = size/2
R = size/3
a = 2*R
bb = 2*R
c = R/5
I1 = b.M*(bb*bb+c*c)/20 + b.M*z0*z0
I2 = b.M*(a*a+c*c)/20 + b.M*z0*z0
I3 = b.M*(bb*bb+a*a)/20
b.Ix = I1
b.Iy = I2
b.Iz = I3
b.freq= 20/2/np.pi
print R
b.arm = z0


b.samplerate=500
b.orien = np.array([-np.radians(60),0,0])
b.GenerateDependentVariables()
#set case here(after dep var is created)
b.setcase(5,0,-0.5)
#b.HasbunEulerEquationODEsolve()
#b.IncludeNoiseInOmega()
b.EulerDCMiter()
#b.directDCMiter()
b.DrawOption['A_axes'] = True
b.DrawOption['A_square'] = True
b.DrawOption['A_z_axis_trace'] = True
#b.DrawOption['B_axes'] = True
#b.DrawOption['B_z_axis_trace'] = True
#b.DrawOption['C_axes'] = True
#b.DrawOption['C_z_axis_trace'] = True

#Draw the trace of the angular momentum vector
b.DrawOption['A_Angular Momentum Trace'] = True
#Draw Angular momentum vector normalized to initial value at t0
b.DrawOption['A_Angular Momentum Vec'] = True
b.DrawOption['A_Angular Velocity Vec'] = True
b.framerate = 30
b.view_azim_angle = 160

iframe = 1000
fig2 = plt.figure(2,figsize=(7,7))
ax4 = p3.Axes3D(fig2)
RBPlot.AnimationFrameSetUp(b,ax4,plt,iframe)
RBPlot.plot_body_space_cone(ax4,b.cordvec[iframe,:,2],b.w_lab_norm[iframe,:],b.L_plot[iframe])
#uncomment here to save z motion for drawing gyro logo
#np.save('outfile', b.cordvec[::20,:20,2])
'''
line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, 
         #                          [0,608],
                                  list(range(1,b.N,b.framerate)),
                                  fargs=b.linesarg,interval=100, blit=False,repeat=False)
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\AnimationFig.pgf')
#line_ani.save('EulerTop_sharp.mp4',writer = 'ffmpeg',fps='24')
'''
plt.show()