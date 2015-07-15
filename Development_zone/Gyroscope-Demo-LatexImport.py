import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunctionV3 as RBPlot
import RGCordTransV13 as RG

b = RG.RigidBodyObject()
#set para here
b.freq= 20
b.GenerateDependentVariables()
#set case here(after dep var is created)
b.setcase(1)
b.EulerDCMiter()
b.DrawOption['A_axes'] = True
b.DrawOption['A_square'] = True
b.DrawOption['A_z_axis_trace'] = True
b.DrawOption['A_Angular Momentum Trace'] = True
b.DrawOption['A_Angular Momentum Vec'] = True
b.DrawOption['A_Angular Velocity Vec'] = True
b.framerate = 20

fig2 = plt.figure(2,figsize=(7,7))
ax4 = p3.Axes3D(fig2)
RBPlot.AnimationFrameSetUp(b,ax4,plt)

line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, 
 list(range(1,b.N,b.framerate)), fargs=b.linesarg,interval=100, 
 blit=False,repeat=True)

plt.show()