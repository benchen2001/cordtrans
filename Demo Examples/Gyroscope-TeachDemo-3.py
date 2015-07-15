import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import sys
sys.path.insert(0, '../Physics_and_Animation_Modules_Library')
import RBPlotFunctionV3 as RBPlot
import RGCordTransV13 as RG

b = RG.RigidBodyObject()
#set para here
b.tn=1.0
b.arm = 0.2  #: location of center of mass of top from origin in meters
b.counter_weight = 0.0 #: mass doesn't spin along symmetry axis
b.counter_weight_location_from_origin = 0.0 #: location from origin
b.Iy=0.002
b.Iz= 0.0008
b.Ix=0.002
#b.freq= 5.0
#b.g = 12.0
b.samplerate=1500
b.UsePY_ODE=0
b.orien = np.array([-np.radians(55),0,0])
b.GenerateDependentVariables()
#set case here(after dep var is created)
b.setcase(3,0,-4)
b.HasbunEulerEquationODEsolve()
#b.IncludeNoiseInOmega()
b.EulerDCMiter()
#b.directDCMiter()
b.DrawOption['A_axes'] = True
b.DrawOption['A_square'] = True
b.DrawOption['A_z_axis_trace'] = True
b.DrawOption['B_axes'] = True
b.DrawOption['B_z_axis_trace'] = True
#b.DrawOption['C_axes'] = True

#Draw the trace of the angular momentum vector
#b.DrawOption['A_Angular Momentum Trace'] = True
#Draw Angular momentum vector normalized to initial value at t0
#b.DrawOption['A_Angular Momentum Vec'] = True
#b.DrawOption['A_Angular Velocity Vec'] = True
b.framerate = 5
b.view_azim_angle = 160

fig2 = plt.figure(2,figsize=(7,7))
ax4 = p3.Axes3D(fig2)
RBPlot.AnimationFrameSetUp(b,ax4,plt,10)

#uncomment here to save z motion for drawing gyro logo
#np.save('outfile', b.cordvec[::20,:20,2])

line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, 
         #                          [0,608],
                                  list(range(1,b.N,b.framerate)),
                                  fargs=b.linesarg,interval=100, blit=False,repeat=True)

#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\AnimationFig.pgf')
#line_ani.save('EulerTop_sharp.mp4',writer = 'ffmpeg',fps='24')
plt.show()