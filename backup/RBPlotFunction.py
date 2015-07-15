# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 16:05:05 2014

@author: user
"""
import numpy as np
print('RB animation and draw functions called')

polynum = [3,4,5,6,3]   # indexing vertices of a closed rectangle

def get111norm(xyzaxes):
    # Return array in size (N,3)
    #print(RBO.cordvec[:10,:3,:3])
    vec = np.sum(xyzaxes,axis=2)
    vecnorm = np.array([vec[i,:]/np.linalg.norm(vec[i,:]) for i in range(len(vec[:,0]))])    
    return vecnorm

def CalDiffInDeg(cordvec):
    DCM111norm = get111norm(cordvec[:,:,:3])
    Has111norm = get111norm(cordvec[:,:,7:10])
    print(np.shape(DCM111norm),np.shape(Has111norm))
    diffindeg = np.degrees([np.arccos(np.dot(DCM111norm[i,:],Has111norm[i,:])) for i in range(len(DCM111norm[:,0])-1)])
    print(np.shape(diffindeg))
    return diffindeg

def CalDiffInDegBC(cordvec):
    dDCM111norm = get111norm(cordvec[:,:,10:13])
    Has111norm = get111norm(cordvec[:,:,7:10])
    print(np.shape(dDCM111norm),np.shape(Has111norm))
    BCdiffindeg = np.degrees([np.arccos(np.dot(dDCM111norm[i,:],Has111norm[i,:])) for i in range(len(dDCM111norm[:,0])-1)])
    print(np.shape(BCdiffindeg))
    return BCdiffindeg


def update_line(x, blines,lineL,lineW,cordvec,L_plot,w_b_norm):
    # update top animation
    for ind,bline in enumerate(blines):
        if ind < 3:	# update baxes body axes
            for blinei,lineLi in zip(bline,lineL): 
                if ind == 2:                
                    blinei.set_data([0,(cordvec[x,0,ind]+cordvec[x,0,2])/2],
                                 [0,(cordvec[x,1,ind]+cordvec[x,1,2])/2])                
                    blinei.set_3d_properties([0,(cordvec[x,2,ind]+cordvec[x,2,2])/2])
                else:
                    blinei.set_data(
                    [0+cordvec[x,0,2]/2 , (cordvec[x,0,ind]/2+cordvec[x,0,2])/2],
                    [0+cordvec[x,1,2]/2 , (cordvec[x,1,ind]/2+cordvec[x,1,2])/2]) 
                    blinei.set_3d_properties(
                    [0+cordvec[x,2,2]/2 , (cordvec[x,2,ind]/2+cordvec[x,2,2])/2])
##            baxisi.set_data([0,cordvec[x,0,0]],[0,cordvec[x,1,0]])
##            baxisi.set_3d_properties([0,cordvec[x,2,0]])
                lineLi.set_data([0,L_plot[x,0]],[0,L_plot[x,1]])
                lineLi.set_3d_properties([0,L_plot[x,2]])
        if 3<=ind<7:           # update rectangle
            for blinei,lineLi in zip(bline,lineL):
                blinei.set_data([cordvec[x,0,polynum[ind-3]],cordvec[x,0,polynum[ind+1-3]]]
                                ,[cordvec[x,1,polynum[ind-3]],cordvec[x,1,polynum[ind+1-3]]])
                blinei.set_3d_properties([cordvec[x,2,polynum[ind-3]]
                                          ,cordvec[x,2,polynum[ind+1-3]]])	
        if ind>=7:     # update hasbun's body axes    
            for blinei,lineLi in zip(bline,lineL): 
                if ind == 9:                
                    blinei.set_data([0,(cordvec[x,0,ind]+cordvec[x,0,9])/2],
                                 [0,(cordvec[x,1,ind]+cordvec[x,1,9])/2])                
                    blinei.set_3d_properties([0,(cordvec[x,2,ind]+cordvec[x,2,9])/2])
                else:
                    blinei.set_data(
                    [0+cordvec[x,0,9]/2 , (cordvec[x,0,ind]/2+cordvec[x,0,9])/2],
                    [0+cordvec[x,1,9]/2 , (cordvec[x,1,ind]/2+cordvec[x,1,9])/2]) 
                    blinei.set_3d_properties(
                    [0+cordvec[x,2,9]/2 , (cordvec[x,2,ind]/2+cordvec[x,2,9])/2])
    #plotting angular velocity vector norm to w[0,2]
    for lineWi,lineLi in zip(lineW,lineL):
        lineWi.set_data([0,w_b_norm[x,0]],[0,w_b_norm[x,1]])
        lineWi.set_3d_properties([0,w_b_norm[x,2]])
    return blines,lineL,lineWi


def InitFirstAnimationFrame(RBO,ax,plt):
    # plot the initial body xyz axes DCM method
    cordvec,Lzppnorm,L_plot,w_b,w = RBO.cordvec,RBO.Lzppnorm,RBO.L_plot,RBO.w_b,RBO.w
    #RBO.baxes = [ax.plot(*[ [0+cordvec[0,j,2]/2,(cordvec[0,j,i]+cordvec[0,j,2])/2] for j in range(3) ]) for i in range(3)]
    #plt.setp(RBO.baxes[2],marker='o') #set z axis marker
    # unpacking list to tuple and list comprehension, see below
    # baxes[1]=ax.plot([0,cordvec[0,0,0]],[0,cordvec[0,1,0]],[0,cordvec[0,2,0]]), plot x axis line

    # plot the initial body xyz axes from Hasbun's euler angle method
    #RBO.baxes_hasbun = [ax.plot(*[ [0+cordvec[0,j,9]/2,(cordvec[0,j,i]+cordvec[0,j,9])/2] for j in range(3) ]) for i in range(7,10)]
    #plt.setp(RBO.baxes_hasbun,linestyle=':') #set z axis marker

    # Color of the rectangle
    ##
    import matplotlib.cm as mplcm
    import matplotlib.colors as colors
    cm = plt.get_cmap('Oranges')
    cNorm  = colors.Normalize(vmin=-1, vmax=3)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(4)])
    ##
    #

    # plot the initial square box
    '''
    RBO.polylines = [ ax.plot(*[ [cordvec[0,j,polynum[i]],
                                  cordvec[0,j,polynum[i+1]]] for j in range(3) ]
                            ) for i in range(4) ]
    plt.setp(RBO.polylines,lw=5)
    '''
    #change previous w(t_i-1) for w(t i-1) comparison
    linetrace=ax.plot(cordvec[:,0,2],cordvec[:,1,2],cordvec[:,2,2],'b-',
                          markersize=1,label='DCM with $\omega(t_i)$')

    #lineLtrance=ax.plot(L_plot[:,0],L_plot[:,1],L_plot[:,2],'k.',markersize=2)
    #HasbunLtrance=ax.plot(Lzppnorm[:,0],Lzppnorm[:,1],Lzppnorm[:,2],'k-',
    #                 markersize=2,label='$\omega(t_{i+1})$ and Hasbun\'s result')
#    ax.legend()
    '''
    RBO.lineL=ax.plot([0,L_plot[0,0]],
                  [0,L_plot[0,1]],
                  [0,L_plot[0,2]],'k-')
    '''
    #HasbunL=ax.plot([0,L_plot[0,0]],
    #               [0,L_plot[0,1]],
    #               [0,L_plot[0,2]],'k-')

    # plot angular velocity vector normalized to w(t0)--
    RBO.w_b_norm = w_b/w[0,2]
    #print(w[0,2],np.linalg.norm(w_b[0,:]))
    '''    
    RBO.lineW=ax.plot([0,RBO.w_b_norm[0,0]],
                  [0,RBO.w_b_norm[0,1]],
                  [0,RBO.w_b_norm[0,2]],'g-')
    '''
# ----

def xthAnimationFrame(RBO,ax,plt):
    # plot the initial body xyz axes DCM method
    cordvec,Lzppnorm,L_plot,w_b,w = RBO.cordvec,RBO.Lzppnorm,RBO.L_plot,RBO.w_b,RBO.w
    RBO.baxes = [ax.plot(*[ [0+cordvec[0,j,2]/2,(cordvec[0,j,i]+cordvec[0,j,2])/2] for j in range(3) ]) for i in range(3)]
    plt.setp(RBO.baxes[2],marker='o') #set z axis marker
    # unpacking list to tuple and list comprehension, see below
    # baxes[1]=ax.plot([0,cordvec[0,0,0]],[0,cordvec[0,1,0]],[0,cordvec[0,2,0]]), plot x axis line

    # plot the initial body xyz axes from Hasbun's euler angle method
    RBO.baxes_hasbun = [ax.plot(*[ [0+cordvec[0,j,9]/2,(cordvec[0,j,i]+cordvec[0,j,9])/2] for j in range(3) ]) for i in range(7,10)]
    plt.setp(RBO.baxes_hasbun,linestyle=':')

    # Color of the rectangle
    ##
    import matplotlib.cm as mplcm
    import matplotlib.colors as colors
    cm = plt.get_cmap('Oranges')
    cNorm  = colors.Normalize(vmin=-1, vmax=3)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(4)])
    ##
    #

    # plot the initial square box
    RBO.polylines = [ ax.plot(*[ [cordvec[0,j,polynum[i]],
                                  cordvec[0,j,polynum[i+1]]] for j in range(3) ]
                            ) for i in range(4) ]
    plt.setp(RBO.polylines,lw=5)
    
    #change previous w(t_i-1) for w(t i-1) comparison
    linetrace=ax.plot(cordvec[:,0,2],cordvec[:,1,2],cordvec[:,2,2],'b-',
                          markersize=1,label='DCM with $\omega(t_i)$')

    lineLtrance=ax.plot(L_plot[:,0],L_plot[:,1],L_plot[:,2],'g.',markersize=1)
    HasbunZtrance=ax.plot(Lzppnorm[:,0],Lzppnorm[:,1],Lzppnorm[:,2],'k.',
                     markersize=1,label='$\omega(t_{i+1})$ and Hasbun\'s result')
#    ax.legend()
    
    RBO.lineL=ax.plot([0,L_plot[0,0]],
                  [0,L_plot[0,1]],
                  [0,L_plot[0,2]],'k-')
    
    #RBO.HasbunL=ax.plot([0,L_plot[0,0]],
    #               [0,L_plot[0,1]],
     #              [0,L_plot[0,2]],'k-')

    # plot angular velocity vector normalized to w(t0)--
    RBO.w_b_norm = w_b/w[0,2]
    #print(w[0,2],np.linalg.norm(w_b[0,:]))
    RBO.lineW=ax.plot([0,RBO.w_b_norm[0,0]],
                  [0,RBO.w_b_norm[0,1]],
                  [0,RBO.w_b_norm[0,2]],'g-')
    




#def RBanimation(RBO,plt,fig2):
#    #
#    N,baxes,polylines,baxes_hasbun,lineL,lineW,cordvec,L_plot,w_b_norm = RBO.N,RBO.baxes,RBO.polylines,RBO.baxes_hasbun,RBO.lineL,RBO.lineW,RBO.cordvec,RBO.L_plot,RBO.w_b_norm
#    import matplotlib.animation as animation
#    # Keep "line_ani =" part, otherwise won't work, timedanimation need an return to clear the figure
#    line_ani = animation.FuncAnimation(fig2, update_line, list(range(1,N,1)),
#                                           fargs=(baxes+polylines+baxes_hasbun,lineL,lineW,cordvec,L_plot,w_b_norm),
#                                          interval=10, blit=False,repeat=False)
#    #line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')
#                                          
#    #plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\wti_wtiplus1.pgf')
        
