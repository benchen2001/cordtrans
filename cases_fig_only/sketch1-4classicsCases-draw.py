import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunctionV3 as RBPlot
import RGCordTransV13 as RG

b1 = RG.RigidBodyObject()
#set para here
b1.GenerateDependentVariables()
#set case here(after dep var is created)
b1.setcase(1)
b1.HasbunEulerEquationODEsolve()
b1.eulerW2bodyW()
b1.EulerDCMiter()
b1.DrawOption['A_z_axis_trace'] = True

b2 = RG.RigidBodyObject()
#set para here
b2.tn=1.4
b2.view_azim_angle = 160
b2.GenerateDependentVariables()
#set case here(after dep var is created)
b2.setcase(5,6,0)
b2.HasbunEulerEquationODEsolve()
b2.eulerW2bodyW()
b2.EulerDCMiter()
b2.DrawOption['A_z_axis_trace'] = True

b3 = RG.RigidBodyObject()
#set para here
b3.GenerateDependentVariables()
#set case here(after dep var is created)
b3.setcase(5,0,-3)
b3.HasbunEulerEquationODEsolve()
b3.eulerW2bodyW()
b3.EulerDCMiter()
b3.DrawOption['A_z_axis_trace'] = True

b4 = RG.RigidBodyObject()
#set para here
b4.GenerateDependentVariables()
#set case here(after dep var is created)
b4.setcase(4)
b4.HasbunEulerEquationODEsolve()
b4.eulerW2bodyW()
b4.EulerDCMiter()
b4.DrawOption['A_z_axis_trace'] = True


#%%____________________ Plot figure top animation
fig2 = plt.figure(2,figsize=(6,1.5))
gs = gridspec.GridSpec(1, 4)
gs.update(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
ax1 = fig2.add_subplot(gs[0,0],projection='3d')
#ax1 = p3.Axes3D(fig2)
#ax1.view_init(elev=20, azim=135)
#ax1.set_xlim3d([-0.8, 0.8])
#ax1.set_ylim3d([-0.8, 0.8])
#ax1.set_zlim3d([0, 0.8])
RBPlot.AnimationFrameSetUp(b1,ax1,plt)
#ax1.set_xticks([])
#ax1.set_yticks([])
#ax1.set_zticks([])

ax2 = fig2.add_subplot(gs[0,1],projection='3d')
#ax1 = p3.Axes3D(fig2)
#ax2.view_init(elev=20, azim=135)
#ax2.set_xlim3d([-0.8, 0.8])
#ax2.set_ylim3d([-0.8, 0.8])
#ax2.set_zlim3d([0, 0.8])
RBPlot.AnimationFrameSetUp(b2,ax2,plt)
#ax2.set_xticks([])
#ax2.set_yticks([])
#ax2.set_zticks([])

ax3 = fig2.add_subplot(gs[0,2],projection='3d')
##ax1 = p3.Axes3D(fig2)
#ax3.view_init(elev=20, azim=135)
#ax3.set_xlim3d([-0.8, 0.8])
#ax3.set_ylim3d([-0.8, 0.8])
#ax3.set_zlim3d([0, 0.8])
RBPlot.AnimationFrameSetUp(b3,ax3,plt)
#ax3.set_xticks([])
#ax3.set_yticks([])
#ax3.set_zticks([])

ax4 = fig2.add_subplot(gs[0,3],projection='3d')
#ax1 = p3.Axes3D(fig2)
#ax4.view_init(elev=20, azim=135)
#ax4.set_xlim3d([-0.8, 0.8])
#ax4.set_ylim3d([-0.8, 0.8])
#ax4.set_zlim3d([0, 0.8])
RBPlot.AnimationFrameSetUp(b4,ax4,plt)
#ax4.set_xticks([])
#ax4.set_yticks([])
#ax4.set_zticks([])
'''
line_ani = animation.FuncAnimation(fig2, RBPlot.update_line, list(range(1,b.N,15)),
                                   fargs=(b.baxes+b.polylines+b.baxes_hasbun,
                                          b.lineL,b.lineW,b.cordvec,b.L_plot,
                                          b.w_b_norm),interval=10, blit=False,repeat=False)
'''
#RBPlot.RBanimation(b,plt,fig2)
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\FourClassics.pgf')
plt.show()