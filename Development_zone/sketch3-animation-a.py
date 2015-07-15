import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunction as RBPlot
import RGCordTransV11 as RG

b = RG.RigidBodyObject()
#set para here
b.tn=1.45
b.GenerateDependentVariables()
#set case here(after dep var is created)
b.setcase(2)
b.HasbunEulerEquationODEsolve()
b.eulerW2bodyW()
b.EulerDCMiter()


#%%____________________ Plot figure top animation
fig2 = plt.figure(2,figsize=(4,2.66))
ax4 = p3.Axes3D(fig2)
ax4.view_init(elev=20, azim=125)
ax4.set_xlim3d([-0.8, 0.8])
ax4.set_ylim3d([-0.8, 0.8])
ax4.set_zlim3d([0, 0.8])
RBPlot.xthAnimationFrame(b,ax4,plt)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_zticks([])

line_ani = animation.FuncAnimation(fig2, RBPlot.update_line, list(range(1,b.N,15)),
                                   fargs=(b.baxes+b.polylines+b.baxes_hasbun,
                                          b.lineL,b.lineW,b.cordvec,b.L_plot,
                                          b.w_b_norm),interval=10, blit=False,repeat=False)

#RBPlot.RBanimation(b,plt,fig2)
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\AnimationFig.pgf')
plt.show()