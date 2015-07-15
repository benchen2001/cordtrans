import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunctionV3 as RBPlot
import RGCordTransV13 as RG

b = RG.RigidBodyObject()
#set para here
b.tn=6.0
b.freq=15.0
#b.g = 20.0
b.samplerate=2000
#b.orien = np.array([-np.radians(65),0,0])
b.GenerateDependentVariables()
#set case here(after dep var is created)
b.setcase(2,1,-6)
b.HasbunEulerEquationODEsolve()
b.eulerW2bodyW()
b.EulerDCMiter()
b.directDCMiter()
b.DrawOption['A_axes'] = True
b.DrawOption['B_axes'] = False

#Draw the trace of the angular momentum vector
#b.DrawOption['A_Angular Momentum Trace'] = True
#b.DrawOption['B_Angular Momentum Trace'] = False
#Draw Angular momentum vector normalized to initial value at t0
#b.DrawOption['A_Angular Momentum Vec'] = True
#b.DrawOption['B_Angular Momentum Vec'] = True

#b.DrawOption['A_Angular Velocity Vec'] = True

#%%____________________ Plot figure top animation
#fig3 = plt.figure(3)
#plt.plot(np.square(b.cordvec[:,0,2])+np.square(b.cordvec[:,1,2])+np.square(b.cordvec[:,2,2]))
#plt.plot(b.L_lab[:,:2])
#plt.plot(np.dot(b.cordvec[:,:2,2],b.L_lab[:,:2]))]
#x=2
#plt.plot([0,0],b.cordvec[x,:2,2],'r-')
#plt.plot([0,b.cordvec[x,0,2]],[0,b.cordvec[x,1,2]],'r-')
#plt.plot([0,b.L_lab[x,0]],[0,b.L_lab[x,1]],'b-',linewidth=5)
#print b.L_lab[x,:2]

fig2 = plt.figure(2,figsize=(7,7))
ax4 = p3.Axes3D(fig2)
RBPlot.AnimationFrameSetUp(b,ax4,plt)


#line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, list(range(1,b.N,60)),
#                                   fargs=b.linesarg,interval=100, blit=False,repeat=False)

#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\AnimationFig.pgf')
#line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')

plt.show()