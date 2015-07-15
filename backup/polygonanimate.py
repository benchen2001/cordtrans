import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation

def update_line(x, linex,liney,linez,lineL):
    for linexi,lineyi,linezi,lineLi in zip(linex,liney,linez,lineL): #must keep the angles part for zip to work and rip the first element
        linexi.set_data([0,cordvec[x,0,0]],[0,cordvec[x,1,0]])
        linexi.set_3d_properties([0,cordvec[x,2,0]])
        lineyi.set_data([0,cordvec[x,0,1]],[0,cordvec[x,1,1]])
        lineyi.set_3d_properties([0,cordvec[x,2,1]])
        linezi.set_data([0,cordvec[x,0,2]],[0,cordvec[x,1,2]])
        linezi.set_3d_properties([0,cordvec[x,2,2]])
        lineLi.set_data([0,L_plot[x,0]],[0,L_plot[x,1]])
        lineLi.set_3d_properties([0,L_plot[x,2]])
    return linex,liney,linez

fig2 = plt.figure(2)
ax = p3.Axes3D(fig2)
ax.view_init(elev=20, azim=135) 

vtx = [cordvec[0,:,0],cordvec[0,:,1],[0,0,0]]
poly = p3.art3d.Poly3DCollection([vtx])
ax.add_collection3d(poly)

line_ani = animation.FuncAnimation(fig2, update_line, list(range(0,N,round(N/24))),
                                   fargs=(linex,liney,linez,lineL),
                                   interval=10, blit=False,repeat=None)
