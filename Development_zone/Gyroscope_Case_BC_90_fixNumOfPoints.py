#C method test on accuracy at 55 and 90 degress with circular classical motion
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunctionV3 as RBPlot
import RGCordTransV13 as RG

b = RG.RigidBodyObject()
b.orien = np.array([-np.radians(55),0,0])
b.counter_weight = 0.05 #: mass doesn't spin along symmetry axis
b.counter_weight_location_from_origin = 0.1 #: location from origin
b.samplerate=500
b.GenerateDependentVariables()
b.setcase(1)
b.HasbunEulerEquationODEsolve()
b.eulerW2bodyW()
b.EulerDCMiter()
#b.directDCMiter()
b.DrawOption['B_axes'] = True
#b.DrawOption['A_square'] = True
b.DrawOption['B_z_axis_trace'] = True
#b.DrawOption['C_axes'] = True
#b.DrawOption['C_z_axis_trace'] = True
#b.DrawOption['A_Angular Velocity Vec'] = True
b.framerate = 20
b.view_azim_angle=170

print b.N
print np.shape(b.euler_angles_hasbun)
print np.shape(b.w_body_hasbun)
print np.shape(b.w)
print np.shape(b.cordvec[:,:,9])

print b.cordvec[-1,:,9]

fig2 = plt.figure(2,figsize=(2.5,2.5))
ax4 = p3.Axes3D(fig2)
RBPlot.AnimationFrameSetUp(b,ax4,plt)
line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, 
 range(b.N),
 fargs=b.linesarg,interval=100, blit=False,repeat=False)
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\Case_BC_d_90.pgf')
#plt.show()
