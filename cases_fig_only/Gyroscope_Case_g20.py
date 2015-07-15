
import sys
sys.path.insert(0, 'C:/Documents and Settings/user/My Documents/tony/Scripts/cordtrans/')

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunctionV3 as RBPlot
import RGCordTransV13 as RG


b = RG.RigidBodyObject()
#b.freq= 0.001
#b.samplerate=1500
b.g=20
b.orien = np.array([-np.radians(55),0,0])
b.GenerateDependentVariables()
b.setcase(1)
b.EulerDCMiter()
b.DrawOption['A_axes'] = True
b.DrawOption['A_square'] = True
b.DrawOption['A_z_axis_trace'] = True
b.view_azim_angle = 150
fig2 = plt.figure(2,figsize=(2.5,2.5))
ax4 = p3.Axes3D(fig2)
RBPlot.AnimationFrameSetUp(b,ax4,plt)
b.framerate = 20
line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, list(range(1,b.N,b.framerate)),
                                  fargs=b.linesarg,interval=100, blit=False,repeat=True)
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\Case6.pgf')
#line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')
plt.show()

np.savetxt(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\otherstuff\Case_A_g20_z_datafile.txt',
           b.cordvec[::23,:,2], fmt='%.4e', delimiter=',', newline='\n', header='', footer='', comments='# ')

square_data =  np.vstack((b.cordvec[0,:,3],
                b.cordvec[0,:,4],
                b.cordvec[0,:,5],
                b.cordvec[0,:,6],
                b.cordvec[0,:,3],
                b.cordvec[0,:,5],
                b.cordvec[0,:,4],
                b.cordvec[0,:,6]))

np.savetxt(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\otherstuff\Case_A_g20_square_datafile.txt',
           square_data, fmt='%.4e', delimiter=',', newline='\n', header='', footer='', comments='# ')
#
print b.cordvec[0,:,2]