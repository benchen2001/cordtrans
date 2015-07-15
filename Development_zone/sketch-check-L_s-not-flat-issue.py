import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunctionV2 as RBPlot
import RGCordTransV11 as RG

b = RG.RigidBodyObject()
#set para here
b.tn=1.0
#b.g = 30.0
b.samplerate=2000
#b.orien = np.array([-np.radians(65),0,0])
b.GenerateDependentVariables()
#set case here(after dep var is created)
b.setcase(3)
b.HasbunEulerEquationODEsolve()
b.eulerW2bodyW()
b.EulerDCMiter()
b.directDCMiter()


#%%____________________ Plot figure top animation
fig2 = plt.figure(2)#,figsize=(5,5))
#plt.plot(np.square(b.L_plot[:,0])+np.square(b.L_plot[:,1]))
plt.plot([np.linalg.norm(b.L_plot[i,:]) for i in np.arange(b.N)])
#ax4 = p3.Axes3D(fig2)np
#ax4.view_init(elev=20, azim=125)
#ax4.set_xlim3d([-0.8, 0.8])
#ax4.set_ylim3d([-0.8, 0.8])
#ax4.set_zlim3d([-0.8, 0.8])
#RBPlot.xthAnimationFrame(b,ax4,plt)
#RBPlot.getPlotAnimationLines(b)
#ax4.set_xticks([])
#ax4.set_yticks([])
#ax4.set_zticks([])
#
# line3D Object, Linedata
#linesarg = (
#            b.baxes[0],b.AlineCMxAxis,
#            b.baxes[1],b.AlineCMyAxis,
#            b.baxes[2],b.AlineCMzAxis,
#            b.polylines[0],b.SquareEdge3,
#            b.polylines[1],b.SquareEdge4,
#            b.polylines[2],b.SquareEdge5,
#            b.polylines[3],b.SquareEdge6,
#            b.baxes_hasbun[0],b.BlineCMxAxis,
#            b.baxes_hasbun[1],b.BlineCMyAxis,
#            b.baxes_hasbun[2],b.BlineCMzAxis,
#            b.lineL,b.lineLData,
#            b.lineW,b.w_bnormData,
#            b.baxes_C[0],b.ClineCMxAxis,
#            b.baxes_C[1],b.ClineCMyAxis,
#            b.baxes_C[2],b.ClineCMzAxis,
#            )

#line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, list(range(1,b.N,1)),
#                                   fargs=linesarg,interval=10, blit=False,repeat=False)

#RBPlot.RBanimation(b,plt,fig2)
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\AnimationFig.pgf')
#line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')

plt.show()