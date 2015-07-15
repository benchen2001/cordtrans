import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import sys
sys.path.insert(0, 'C:/Documents and Settings/user/My Documents/tony/Scripts/cordtrans/Physics_and_Animation_Modules_Library')

import RBPlotFunctionV5 as RBPlot
import RGCordTransV15 as RG

b = RG.RigidBodyObject()
#set para here
b.tn=1.5
b.samplerate=800
#b.orien = np.array([-np.radians(65),0,0])
b.GenerateDependentVariables()
#set case here(after dep var is created)
b.setcase(2)
#b.HasbunEulerEquationODEsolve()
b.EulerDCMiter()
#b.directDCMiter()
b.DrawOption['A_axes'] = True
b.DrawOption['A_square'] = False
b.DrawOption['A_z_axis_trace'] = True
b.DrawOption['A_cube'] = False
#Draw the trace of the angular momentum vector
b.DrawOption['A_Angular Momentum Trace'] = True
#Draw Angular momentum vector normalized to initial value at t0
b.DrawOption['A_Angular Momentum Vec'] = True
b.DrawOption['A_Angular Velocity Vec'] = True

iframe = 170
bout = b.construct_cube_4wires_method(iframe)


#%%____________________ Plot figure top animation
fig2 = plt.figure(2,figsize=(3.5,3.5))
ax4 = p3.Axes3D(fig2)
b.framerate = 20
b.view_azim_angle = 170
RBPlot.AnimationFrameSetUp(b,ax4,plt,iframe)

cubeinst1 = ax4.plot(*bout[:,0:4],c='black',linewidth=2)
cubeinst2 = ax4.plot(*bout[:,4:8],c='gray',linewidth=2)
cubeinst3 = ax4.plot(*bout[:,8:12],c='lightgray',linewidth=2)
cubeinst4 = ax4.plot(*bout[:,12:16],c='darkgray',linewidth=2)


'''
line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_multiple_segments, 
         #                          [0,608],
                                  list(range(1,b.N,b.framerate)),
                                  fargs=b.linesarg,interval=100, blit=False,repeat=False)
'''
#RBPlot.RBanimation(b,plt,fig2)
plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\AnimationFig.pgf')
#line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')

plt.show()