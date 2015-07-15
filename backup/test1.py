# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 16:05:58 2014

@author: user
"""

def update_line(x, blines, blines2):# update baxes body axes    
    print(x)
    for ind,bline in enumerate(blines):
        for blinei, bline2i in zip(bline,blines2):
            blinei.set_data([0,cordvec[x,0,ind]],
                            [0,cordvec[x,1,ind]]) 
            blinei.set_3d_properties([0,cordvec[x,2,ind]])
    return blines

import DCMiterator
w = np.zeros([200,3])+0.8 
N = len(w[:,0])
cordvec = DCMiterator.DCMiter(w,[np.pi/2,0,0])

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3
import matplotlib.animation as animation
fig2 = plt.figure(2)
ax = p3.Axes3D(fig2)
ax.view_init(elev=20, azim=135)
ax.set_xlim3d([-2, 2])
ax.set_ylim3d([-2, 2])
ax.set_zlim3d([-2, 2])
lineL=ax.plot([0,1],
              [0,-1],
              [0,1],'k-')
baxes = [ax.plot(*[ [0,cordvec[0,j,i]] for j in range(3) ]) for i in range(3)]
# Keep "line_ani =" part, otherwise won't work, timedanimation need an return to clear the figure
line_ani = animation.FuncAnimation(fig2, update_line, list(range(0,N,1)),
                                   fargs=(baxes,lineL), 
                                   interval=100, 
                                   blit=False,
                                   repeat=None) 

  
plt.show()
