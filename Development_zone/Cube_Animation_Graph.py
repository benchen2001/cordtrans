import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
import RBPlotFunctionV5 as RBPlot
import RGCordTransV15 as RG

b = RG.RigidBodyObject()
#set para here
b.tn=0.9
b.arm = 0.2  #: location of center of mass of top from origin in meters
b.counter_weight = 0.0 #: mass doesn't spin along symmetry axis
b.counter_weight_location_from_origin = 0.0 #: location from origin
b.Iy=0.002
b.Iz= 0.0008
b.Ix=0.002
#b.freq= 0.0
#b.g = 12.0
b.samplerate=500
b.orien = np.array([-np.radians(55),0,0])
b.GenerateDependentVariables()
#set case here(after dep var is created)
b.setcase(1,20,-10)
b.DrawOption['A_axes'] = True
b.DrawOption['A_square'] = False
b.DrawOption['A_z_axis_trace'] = True
b.DrawOption['A_cube'] = False

#b.HasbunEulerEquationODEsolve()
#b.IncludeNoiseInOmega()
b.EulerDCMiter()
#b.directDCMiter()
iframe = 40
testout = b.construct_cube_4wires_method(iframe)

#b.DrawOption['B_axes'] = True
#b.DrawOption['B_z_axis_trace'] = True
#b.DrawOption['C_axes'] = True

#Draw the trace of the angular momentum vector
b.DrawOption['A_Angular Momentum Trace'] = True
#Draw Angular momentum vector normalized to initial value at t0
b.DrawOption['A_Angular Momentum Vec'] = True
b.DrawOption['A_Angular Velocity Vec'] = True
b.framerate = 1
b.view_azim_angle = 45

fig2 = plt.figure(2,figsize=(6,6))
ax4 = p3.Axes3D(fig2)

#iframe_cubedata = RBPlot.plot_cube_4wires(ax4,testout,iframe,b.cordvec[iframe,:,2])
cubeinst1 = ax4.plot(*testout[:,0:4],c='black',linewidth=2)
cubeinst2 = ax4.plot(*testout[:,4:8],c='gray',linewidth=2)
cubeinst3 = ax4.plot(*testout[:,8:12],c='lightgray',linewidth=2)
cubeinst4 = ax4.plot(*testout[:,12:16],c='darkgray',linewidth=2)

RBPlot.AnimationFrameSetUp(b,ax4,plt,iframe)
RBPlot.plot_body_space_cone(ax4,b.cordvec[iframe,:,2],b.w_lab_norm[iframe,:],b.L_plot[iframe])


#uncomment here to save z motion for drawing gyro logo
#np.save('outfile', b.cordvec[::20,:20,2])
'''
line_ani = animation.FuncAnimation(fig2, RBPlot.update_line_multiple_segments, 
         #                          [0,608],
                                  list(range(1,b.N,b.framerate)),
                                  fargs=(b,
                                         cubeinst1,
                                         cubeinst2,
                                         cubeinst3,
                                         cubeinst4) + b.linesarg,
                                         interval=100, blit=False,repeat=True)
'''
#plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\AnimationFig.pgf')
#line_ani.save('EulerTop_sharp_cube.mp4',writer = 'ffmpeg',fps='24')
plt.show()