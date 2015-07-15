import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import sys
sys.path.insert(0, 'C:/Documents and Settings/user/My Documents/tony/Scripts/cordtrans/Physics_and_Animation_Modules_Library')
import RBPlotFunctionV3 as RBPlot
import RGCordTransV13 as RG

b = RG.RigidBodyObject()
b.tn=1.2
b.samplerate=2000
b.GenerateDependentVariables()
b.setcase(1)
b.HasbunEulerEquationODEsolve()
#b.eulerW2bodyW()
#b.w_body_hasbun = b.w_body_hasbun*0
b.IncludeNoiseInOmega(10)
b.EulerDCMiter()
b.directDCMiter()
b.DrawOption['C_axes'] = True
b.DrawOption['C_z_axis_trace'] = True
b.DrawOption['B_axes'] = True
b.DrawOption['B_z_axis_trace'] = True

b.framerate = 20
b.view_azim_angle=70
fig2 = plt.figure(2,figsize=(6,6))
ax4 = p3.Axes3D(fig2)
RBPlot.AnimationFrameSetUp(b,ax4,plt,10)
#trailsamplerate = 23
#ax4.plot(b.cordvec[::trailsamplerate,0,11],b.cordvec[::trailsamplerate,1,11],b.cordvec[::trailsamplerate,2,11],'g.',
#                          linewidth=1,label='C method')

line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, 
 list(range(1,b.N,b.framerate)),
 fargs=b.linesarg,interval=100, blit=False,repeat=True)
plt.show()

