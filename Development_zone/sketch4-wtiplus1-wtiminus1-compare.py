import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunctionV3 as RBPlot
import RGCordTransV13 as RG

b1 = RG.RigidBodyObject()
b1.tn=1.0
b1.freq=17.0
b1.samplerate=2000
#b.orien = np.array([-np.radians(65),0,0])
b1.GenerateDependentVariables()
b1.setcase(2)
b1.EulerDCMiter()
b1.DrawOption['A_z_axis_trace'] = True
b1.DrawOption['A_square'] = True
b1.DrawOption['A_axes'] =True

b2 = RG.RigidBodyObject()
b2.tn=1.0
b2.freq=17.0
b2.samplerate=2000
#b.orien = np.array([-np.radians(65),0,0])
b2.GenerateDependentVariables()
b2.setcase(2)
b2.EulerDCMiter_wt_i()
b2.DrawOption['A_z_axis_trace'] = True


fig2 = plt.figure(2,figsize=(4,4))
ax4 = p3.Axes3D(fig2)
RBPlot.AnimationFrameSetUp(b1,ax4,plt)
linetrace=ax4.plot(b2.cordvec[::10,0,2],b2.cordvec[::10,1,2],b2.cordvec[::10,2,2],'k:',
                  linewidth=2,label='A method with $\omega(t_i)$')
ax4.legend(loc=3,prop={'size':10})

#line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, list(range(1,b.N,60)),
#                                   fargs=b.linesarg,interval=100, blit=False,repeat=False)

plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\wti_wtiplus1.pgf')
#line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')

plt.show()