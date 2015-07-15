import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunction as RBPlot
import RGCordTransV11 as RG

b1 = RG.RigidBodyObject()
#set para here
b1.tn = 20.0
b1.GenerateDependentVariables()
#set case here(after dep var is created)
b1.setcase(1)
b1.HasbunEulerEquationODEsolve()
b1.eulerW2bodyW()
b1.EulerDCMiter()
b1.directDCMiter()

b2 = RG.RigidBodyObject()
#set para here
b2.tn = 5
b2.GenerateDependentVariables()
#set case here(after dep var is created)
b2.setcase(2)
b2.HasbunEulerEquationODEsolve()
b2.eulerW2bodyW()
b2.EulerDCMiter()
b2.directDCMiter()

b3 = RG.RigidBodyObject()
#set para here
b3.tn = 1
b3.GenerateDependentVariables()
#set case here(after dep var is created)
b3.setcase(3)
b3.HasbunEulerEquationODEsolve()
b3.eulerW2bodyW()
b3.EulerDCMiter()
b3.directDCMiter()

b4 = RG.RigidBodyObject()
#set para here
b4.tn = 1
b4.GenerateDependentVariables()
#set case here(after dep var is created)
b4.setcase(4)
b4.HasbunEulerEquationODEsolve()
b4.eulerW2bodyW()
b4.EulerDCMiter()
b4.directDCMiter()

#%%____________________ Plot compare EOM+DCM vs Euler
fig2 = plt.figure(2,figsize=(6.4,4))
gs = gridspec.GridSpec(2, 2)
gs.update(left=0.075, bottom=0.125, right=0.95, top=0.90, wspace=0.15, hspace=0.2)
ax1 = fig2.add_subplot(gs[0,0])
ax2 = fig2.add_subplot(gs[1,0])
ax3 = fig2.add_subplot(gs[0,1])
ax4 = fig2.add_subplot(gs[1,1])

ax1.plot(np.linspace(0.0,b1.tn,b1.N/100),RBPlot.CalDiffInDeg(b1.cordvec)[::100])
ax1.plot(np.linspace(0.0,b1.tn,b1.N/100),RBPlot.CalDiffInDegBC(b1.cordvec)[::100])
ax1.locator_params(tight=True, nbins=4)
ax1.text(0.1,0.8,'(a)',fontsize=9,transform=ax1.transAxes)
ax2.plot(np.linspace(0.0,b2.tn,b2.N/100),RBPlot.CalDiffInDeg(b2.cordvec)[::100])
ax2.plot(np.linspace(0.0,b2.tn,b2.N/100),RBPlot.CalDiffInDegBC(b2.cordvec)[::100])
ax2.locator_params(tight=True, nbins=4)
ax2.text(0.1,0.8,'(b)',fontsize=9,transform=ax2.transAxes)
ax3.plot(np.linspace(0.0,b3.tn,b3.N/10),RBPlot.CalDiffInDeg(b3.cordvec)[::10])
ax3.plot(np.linspace(0.0,b3.tn,b3.N/10),RBPlot.CalDiffInDegBC(b3.cordvec)[::10])
ax3.locator_params(tight=True, nbins=4)
ax3.text(0.1,0.8,'(c)',fontsize=9,transform=ax3.transAxes)
ax4.plot(np.linspace(0.0,b4.tn,b4.N/10),RBPlot.CalDiffInDeg(b4.cordvec)[::10])
ax4.plot(np.linspace(0.0,b4.tn,b4.N/10),RBPlot.CalDiffInDegBC(b4.cordvec)[::10])
ax4.locator_params(tight=True, nbins=4)
ax4.text(0.1,0.8,'(d)',fontsize=9,transform=ax4.transAxes)

plt.suptitle(r'Angle diff in the two body z-axis',fontsize=10)
fig2.text(0.025, 0.5, 'degress',fontsize=10, va='center', rotation='vertical')
fig2.text(0.5, 0.025, 'elapsed time in second',fontsize=10, ha='center')

plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\CompareDCMEulerAngle.pgf')


'''
ax1.view_init(elev=20, azim=135)
ax1.set_xlim3d([-0.8, 0.8])
ax1.set_ylim3d([-0.8, 0.8])
ax1.set_zlim3d([0, 0.8])

ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_zticks([])

ax2 = fig2.add_subplot(gs[0,1],projection='3d')
#ax1 = p3.Axes3D(fig2)
ax2.view_init(elev=20, azim=135)
ax2.set_xlim3d([-0.8, 0.8])
ax2.set_ylim3d([-0.8, 0.8])
ax2.set_zlim3d([0, 0.8])
RBPlot.InitFirstAnimationFrame(b2,ax2,plt)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_zticks([])

ax3 = fig2.add_subplot(gs[0,2],projection='3d')
#ax1 = p3.Axes3D(fig2)
ax3.view_init(elev=20, azim=135)
ax3.set_xlim3d([-0.8, 0.8])
ax3.set_ylim3d([-0.8, 0.8])
ax3.set_zlim3d([0, 0.8])
RBPlot.InitFirstAnimationFrame(b3,ax3,plt)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_zticks([])
ax4 = fig2.add_subplot(gs[0,3],projection='3d')
#ax1 = p3.Axes3D(fig2)
ax4.view_init(elev=20, azim=135)
ax4.set_xlim3d([-0.8, 0.8])
ax4.set_ylim3d([-0.8, 0.8])
ax4.set_zlim3d([0, 0.8])
RBPlot.InitFirstAnimationFrame(b4,ax4,plt)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_zticks([])

line_ani = animation.FuncAnimation(fig2, RBPlot.update_line, list(range(1,b.N,15)),
                                   fargs=(b.baxes+b.polylines+b.baxes_hasbun,
                                          b.lineL,b.lineW,b.cordvec,b.L_plot,
                                          b.w_b_norm),interval=10, blit=False,repeat=False)

#RBPlot.RBanimation(b,plt,fig2)
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\FourClassics.pgf')
plt.show()
'''