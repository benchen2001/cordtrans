import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunction as RBPlot
import RGCordTransV11 as RG

b3 = RG.RigidBodyObject()
#set para here
b3.tn = 10.0
b3.samplerate = 2000
b3.GenerateDependentVariables()
#set case here(after dep var is created)
b3.setcase(1)
b3.HasbunEulerEquationODEsolve()
b3.eulerW2bodyW()
b3.EulerDCMiter()
b3.directDCMiter()

#%%____________________ Plot compare EOM+DCM vs Euler
fig3 = plt.figure(3,figsize=(5,4))

diffb3 = [np.arccos(np.dot(b3.cordvec[i,:,9],
                              b3.cordvec[i,:,12]/np.linalg.norm(b3.cordvec[i,:,12]))) for i in range(0,b3.N,1)]
diffindegb3 = np.degrees(diffb3)

diffb3AB = [np.arccos(np.dot(b3.cordvec[i,:,9],
                              b3.cordvec[i,:,2]/np.linalg.norm(b3.cordvec[i,:,2]))) for i in range(0,b3.N,1)]
diffindegb3AB = np.degrees(diffb3AB)



plt.plot(np.linspace(0.0,b3.tn,b3.N/1),diffindegb3[::1])
plt.plot(np.linspace(0.0,b3.tn,b3.N/1),diffindegb3AB[::1])

#plt.locator_params(tight=True, nbins=4)
#plt.text(0.1,0.8,'(c)',fontsize=9,transform=ax3.transAxes)
#plt.suptitle(r'Angle diff in the two body z-axis',fontsize=10)
#fig2.text(0.04, 0.5, 'degress',fontsize=10, va='center', rotation='vertical')
#fig2.text(0.5, 0.04, 'elapsed time in second',fontsize=10, ha='center')
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\AB_methods_cdCases_zaxis_compare.pgf')


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