import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import sys
sys.path.insert(0, 'C:/Documents and Settings/user/My Documents/tony/Scripts/cordtrans/Physics_and_Animation_Modules_Library')
import RBPlotFunctionV4 as RBPlot
import RGCordTransV15 as RG

b = RG.RigidBodyObject()
#set para here
b.tn=1.5
b.arm = 0.2  #: location of center of mass of top from origin in meters
b.counter_weight = 0.0 #: mass doesn't spin along symmetry axis
b.counter_weight_location_from_origin = 0.0 #: location from origin
b.Iy=0.002
b.Iz= 0.0008
b.Ix=0.002
b.freq= 0.1
b.g = 20.0
b.samplerate=1000
b.orien = np.array([-np.radians(0.01),0,0])
b.GenerateDependentVariables()
#set case here(after dep var is created)
b.setcase(1,20,-10)
b.PY_ODE=1
b.HasbunEulerEquationODEsolve()
#b.IncludeNoiseInOmega()
b.EulerDCMiter()
#b.directDCMiter()
#b.DrawOption['A_axes'] = True
#b.DrawOption['A_square'] = True
b.DrawOption['B_z_axis_trace'] = True
b.DrawOption['B_axes'] = True
b.DrawOption['A_z_axis_trace'] = True
b.DrawOption['A_axes'] = True
b.DrawOption['B_Angular Velocity Vector (normalized to t0 value)'] = True
b.DrawOption['A_Angular Velocity Vec'] = True

#Draw the trace of the angular momentum vector
#b.DrawOption['A_Angular Momentum Trace'] = True
#Draw Angular momentum vector normalized to initial value at t0
#b.DrawOption['A_Angular Momentum Vec'] = True
#b.DrawOption['A_Angular Velocity Vec'] = True
b.framerate = 10
b.view_azim_angle = 160

fig2 = plt.figure(2,figsize=(7,7))
ax4 = p3.Axes3D(fig2)
RBPlot.AnimationFrameSetUp(b,ax4,plt,10)

ax4.plot(*np.transpose(b.w_lab_norm[::23,:]),linewidth=3)

#uncomment here to save z motion for drawing gyro logo
#np.save('outfile', b.cordvec[::20,:20,2])

line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, 
         #                          [0,608],
                                  list(range(1,b.N,b.framerate)),
                                  fargs=b.linesarg,interval=100, blit=False,repeat=False)

#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\AnimationFig.pgf')
#line_ani.save('EulerTop_sharp.mp4',writer = 'ffmpeg',fps='24')
plt.show()