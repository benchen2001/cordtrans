import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunctionV3 as RBPlot
import RGCordTransV13 as RG

b = RG.RigidBodyObject()
b.tn=4
b.samplerate=500
b.GenerateDependentVariables()
b.setcase(1)
b.HasbunEulerEquationODEsolve()
b.w_body_hasbun = b.w_body_hasbun*0
b.IncludeNoiseInOmega(5)
b.EulerDCMiter()
b.directDCMiter()
b.DrawOption['C_axes'] = True
b.DrawOption['C_z_axis_trace'] = True
b.framerate = 200
b.view_azim_angle=70
fig2 = plt.figure(2,figsize=(2.5,2.5))
ax4 = p3.Axes3D(fig2)
RBPlot.AnimationFrameSetUp(b,ax4,plt)

b2 = RG.RigidBodyObject()
b2.tn=4
b2.samplerate=500
b2.GenerateDependentVariables()
b2.setcase(1)
b2.HasbunEulerEquationODEsolve()
#b2.eulerW2bodyW()
b2.w_body_hasbun = b.w_body_hasbun*0
b2.IncludeNoiseInOmega(1)
b2.EulerDCMiter()
b2.directDCMiter()
b2.DrawOption['C_axes'] = True
b2.DrawOption['C_z_axis_trace'] = True
b2.framerate = 200
b2.view_azim_angle=70

#line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_new, 
# list(range(1,b.N,b.framerate)),
# fargs=b.linesarg,interval=100, blit=False,repeat=True)
#plt.show()

np.savetxt(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\otherstuff\Case_C_stillnoise_amp2_datafilez.txt',
           b.cordvec[::23,:,12], fmt='%.4e', delimiter=',', newline='\n', header='', footer='', comments='# ')
np.savetxt(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\otherstuff\Case_C_stillnoise_amp2_datafiley.txt',
           b.cordvec[::23,:,11], fmt='%.4e', delimiter=',', newline='\n', header='', footer='', comments='# ')
np.savetxt(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\otherstuff\Case_C_stillnoise_amp2_datafilex.txt',
           b.cordvec[::23,:,10], fmt='%.4e', delimiter=',', newline='\n', header='', footer='', comments='# ')
'''
C_CMxaxis=(b.cordvec[0,:,10]/2+b.cordvec[0,:,12])/2
C_CMyaxis=(b.cordvec[0,:,11]/2+b.cordvec[0,:,12])/2

square_data =  np.vstack((b.cordvec[0,:,3],
                b.cordvec[0,:,4],
                b.cordvec[0,:,5],
                b.cordvec[0,:,6],
                b.cordvec[0,:,3],
                b.cordvec[0,:,5],
                b.cordvec[0,:,4],
                b.cordvec[0,:,6]))
'''
np.savetxt(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\otherstuff\Case_C_stillnoise_amp1_datafilez.txt',
           b2.cordvec[::23,:,12], fmt='%.4e', delimiter=',', newline='\n', header='', footer='', comments='# ')
np.savetxt(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\otherstuff\Case_C_stillnoise_amp1_datafiley.txt',
           b2.cordvec[::23,:,11], fmt='%.4e', delimiter=',', newline='\n', header='', footer='', comments='# ')
np.savetxt(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\otherstuff\Case_C_stillnoise_amp1_datafilex.txt',
           b2.cordvec[::23,:,10], fmt='%.4e', delimiter=',', newline='\n', header='', footer='', comments='# ')


#print b2.cordvec[::23,:,10]
#print b2.cordvec[::23,:,11]
print b.cordvec[0,:,10]
print b.cordvec[0,:,11]