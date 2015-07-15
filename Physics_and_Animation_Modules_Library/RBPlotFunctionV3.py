# -*- coding: utf-8 -*-
"""
Plot and animation functions/environment setup for Rigid Body class

Created on Tue Nov 11 16:05:05 2014

@author: whymranderson.blogspot.tw
"""
import numpy as np
print('RB animation and draw functions called')

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d  #<- uncomment this, you need this

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

polynum = [3,4,5,6,3]   # indexing vertices of a closed rectangle

def CK(rotvec):
    '''Cayley-Klein parameters. Important: the built rotation matrix is with its couterclockwise active sense.
    This is basically the Rodriguez rotation formula in a matric form.
    '''
    amp = np.sqrt(np.dot(rotvec,rotvec))
    if amp == 0:
        ret = np.eye(3)
    else:
        axis = rotvec/amp
        phi = amp % (2*np.pi)
        a = np.cos(phi/2)
        b,c,d = axis*np.sin(phi/2)
        ret =  np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
    return ret

def get111norm(xyzaxes):
    '''Return array of body (1,1,1) unit vector as a function of time N in size (N,3)
    '''    
    vec = np.sum(xyzaxes,axis=2)
    vecnorm = np.array([vec[i,:]/np.linalg.norm(vec[i,:]) for i in range(len(vec[:,0]))])    
    return vecnorm

def CalDiffInDeg(cordvec):
    '''Calculate and return angle difference between A and B method's 111 unit vectors'''
    DCM111norm = get111norm(cordvec[:,:,:3])
    Has111norm = get111norm(cordvec[:,:,7:10])
    print(np.shape(DCM111norm),np.shape(Has111norm))
    diffindeg = np.degrees([np.arccos(np.dot(DCM111norm[i,:],Has111norm[i,:])) for i in range(len(DCM111norm[:,0])-1)])
#    print(np.shape(diffindeg))
    return diffindeg

def CalDiffInDegBC(cordvec):
    '''Calculate and return angle difference between B and C method's 111 unit vectors'''
    dDCM111norm = get111norm(cordvec[:,:,10:13])
    Has111norm = get111norm(cordvec[:,:,7:10])
    print(np.shape(dDCM111norm),np.shape(Has111norm))
    BCdiffindeg = np.degrees([np.arccos(np.dot(dDCM111norm[i,:],Has111norm[i,:])) for i in range(len(dDCM111norm[:,0])-1)])
    print(np.shape(BCdiffindeg))
    return BCdiffindeg

def getPlotAnimationLines(RBO):
    '''Construct and return the line animation data sequence in a rigid body object for 3D animation'''
    RBO.AlineCMxAxis = np.hstack([RBO.cordvec[:,:,2]/2,(RBO.cordvec[:,:,0]/2+RBO.cordvec[:,:,2])/2])
    RBO.AlineCMyAxis = np.hstack([RBO.cordvec[:,:,2]/2,(RBO.cordvec[:,:,1]/2+RBO.cordvec[:,:,2])/2])
    RBO.AlineCMzAxis = np.hstack([np.zeros((RBO.N+1,3)),RBO.cordvec[:,:,2]])    
    RBO.SquareEdge3 = np.hstack([RBO.cordvec[:,:,3],RBO.cordvec[:,:,4]])
    RBO.SquareEdge4 = np.hstack([RBO.cordvec[:,:,4],RBO.cordvec[:,:,5]])
    RBO.SquareEdge5 = np.hstack([RBO.cordvec[:,:,5],RBO.cordvec[:,:,6]])
    RBO.SquareEdge6 = np.hstack([RBO.cordvec[:,:,6],RBO.cordvec[:,:,3]])
    RBO.BlineCMxAxis = np.hstack([RBO.cordvec[:,:,9]/2,(RBO.cordvec[:,:,7]/2+RBO.cordvec[:,:,9])/2])
    RBO.BlineCMyAxis = np.hstack([RBO.cordvec[:,:,9]/2,(RBO.cordvec[:,:,8]/2+RBO.cordvec[:,:,9])/2])
    RBO.BlineCMzAxis = np.hstack([np.zeros((RBO.N+1,3)),RBO.cordvec[:,:,9]])
    RBO.lineLData = np.hstack([np.zeros((RBO.N+1,3)),RBO.L_plot])
    RBO.w_bnormData = np.hstack([np.zeros((RBO.N+1,3)),RBO.w_lab_norm])
    RBO.ClineCMxAxis = np.hstack([RBO.cordvec[:,:,12]/2,(RBO.cordvec[:,:,10]/2+RBO.cordvec[:,:,12])/2])
    RBO.ClineCMyAxis = np.hstack([RBO.cordvec[:,:,12]/2,(RBO.cordvec[:,:,11]/2+RBO.cordvec[:,:,12])/2])
    RBO.ClineCMzAxis = np.hstack([np.zeros((RBO.N+1,3)),RBO.cordvec[:,:,12]])


sphere_color = 'lightgray'#'#FFDDDD'
sphere_alpha = 0.1
frame_color = 'darkgray'
frame_alpha = 0.4
frame_width = 1
sphere_grid = 25
wirestride = 3
def plot_front(axes):
    '''plot the lucid sphere, front part'''
    u = np.linspace(-np.pi, 0, sphere_grid)
    v = np.linspace(0, np.pi, sphere_grid)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    axes.plot_surface(x, y, z, rstride=2, cstride=2,
                           color=sphere_color, linewidth=0,
                           alpha=sphere_alpha)
    # wireframe
    axes.plot_wireframe(x, y, z, rstride=wirestride, cstride=wirestride,
                             color=frame_color,
                             alpha=frame_alpha)
#    # equator
#    axes.plot(1.0 * np.cos(u), 1.0 * np.sin(u),
#                   zs=0, zdir='z', lw=frame_width,
#                   color=frame_color)
#    axes.plot(1.0 * np.cos(u), 1.0 * np.sin(u),
#                   zs=0, zdir='x', lw=frame_width,
#                   color=frame_color)

def plot_back(axes):
    '''plot the lucid sphere, back part'''
    u = np.linspace(0, np.pi, sphere_grid)
    v = np.linspace(0, np.pi, sphere_grid)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    axes.plot_surface(x, y, z, rstride=2, cstride=2,
                           color=sphere_color, linewidth=0,
                           alpha=sphere_alpha)
    # wireframe
    axes.plot_wireframe(x, y, z, rstride=wirestride, cstride=wirestride,
                             color=frame_color,
                             alpha=frame_alpha)
#    # equator
#    axes.plot(1.0 * np.cos(u), 1.0 * np.sin(u),
#                   zs=0, zdir='z', lw=frame_width,
#                   color=frame_color)
#    axes.plot(1.0 * np.cos(u), 1.0 * np.sin(u),
#                   zs=0, zdir='x', lw=frame_width,
#                   color=frame_color)


def update_line_new(x,*arg):
    '''line animation function, will be called by functionanimation, take multiple line data sequences as arguments'''
    arglen = len(arg)
    for i in np.arange(0,arglen,2):
        line3DObject = arg[i]
        line3DData = arg[i+1]   
        for linei,lineid in zip(line3DObject,line3DObject):
            linei.set_data([line3DData[x,0],line3DData[x,3]],
                           [line3DData[x,1],line3DData[x,4]])                 
            linei.set_3d_properties([line3DData[x,2],line3DData[x,5]])  


def InitFirstAnimationFrame(RBO,ax,plt):
    # plot the initial body xyz axes DCM method
    cordvec,Lzppnorm,L_plot,w_lab,w = RBO.cordvec,RBO.Lzppnorm,RBO.L_plot,RBO.w_lab,RBO.w
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
    RBO.w_lab_norm = w_lab/w[0,2]
    #print(w[0,2],np.linalg.norm(w_lab[0,:]))
    '''    
    RBO.lineW=ax.plot([0,RBO.w_lab_norm[0,0]],
                  [0,RBO.w_lab_norm[0,1]],
                  [0,RBO.w_lab_norm[0,2]],'g-')
    '''
# ----

def AnimationFrameSetUp(RBO,ax,plt,istep):
    '''plot the first frame of animation'''
    ax.set_axis_off()
    ax.view_init(elev=20, azim=RBO.view_azim_angle)
    ax.set_xlim3d([-0.57, 0.57])
    ax.set_ylim3d([-0.57, 0.57])
    ax.set_zlim3d([-0.57, 0.57])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    #ax.autoscale_view(tight=True, scalex=True, scaley=True, scalez=True)
    plot_front(ax)
    plot_back(ax)

    # plot the initial body xyz axes DCM method
    cordvec = RBO.cordvec
    if RBO.DrawOption['A_axes'] == True:
        RBO.baxes =[ax.plot(*[ [cordvec[istep,j,2]/2,(cordvec[istep,j,i]/4+cordvec[istep,j,2])/2] for j in range(3) ]) for i in range(2)] + [ax.plot(*[ [0,cordvec[istep,j,2]] for j in range(3) ])]
        plt.setp(RBO.baxes[2],marker='o') #set z axis marker
        # unpacking list to tuple and list comprehension, see below
        # baxes[1]=ax.plot([0,cordvec[0,0,0]],[0,cordvec[0,1,0]],[0,cordvec[0,2,0]]), plot x axis line
        RBO.AlineCMxAxis = np.hstack([RBO.cordvec[:,:,2]/2,(RBO.cordvec[:,:,0]/4+RBO.cordvec[:,:,2])/2])
        RBO.AlineCMyAxis = np.hstack([RBO.cordvec[:,:,2]/2,(RBO.cordvec[:,:,1]/4+RBO.cordvec[:,:,2])/2])
        RBO.AlineCMzAxis = np.hstack([np.zeros((RBO.N+1,3)),RBO.cordvec[:,:,2]])
        RBO.linesarg = RBO.linesarg + (RBO.baxes[0],RBO.AlineCMxAxis,
                                       RBO.baxes[1],RBO.AlineCMyAxis,
                                        RBO.baxes[2],RBO.AlineCMzAxis)


    # plot the initial body xyz axes from Hasbun's euler angle method
    if RBO.DrawOption['B_axes'] == True:
        RBO.baxes_hasbun = [ax.plot(*[ [0+cordvec[0,j,9]/2,(cordvec[0,j,i]+cordvec[0,j,9])/2] for j in range(3) ]) for i in range(7,10)]
        plt.setp(RBO.baxes_hasbun,linestyle=':')
        RBO.BlineCMxAxis = np.hstack([RBO.cordvec[:,:,9]/2,(RBO.cordvec[:,:,7]/2+RBO.cordvec[:,:,9])/2])
        RBO.BlineCMyAxis = np.hstack([RBO.cordvec[:,:,9]/2,(RBO.cordvec[:,:,8]/2+RBO.cordvec[:,:,9])/2])
        RBO.BlineCMzAxis = np.hstack([np.zeros((RBO.N+1,3)),RBO.cordvec[:,:,9]])
        RBO.linesarg = RBO.linesarg + (RBO.baxes_hasbun[0],RBO.BlineCMxAxis,
            RBO.baxes_hasbun[1],RBO.BlineCMyAxis,
            RBO.baxes_hasbun[2],RBO.BlineCMzAxis)

    if RBO.DrawOption['C_axes'] == True:
        RBO.baxes_C =[ax.plot(*[ [cordvec[0,j,12]/2,(cordvec[0,j,i]/2+cordvec[0,j,9])/2] for j in range(3) ]) for i in range(2)] + [ax.plot(*[ [0,cordvec[0,j,12]] for j in range(3) ])]
#        RBO.baxes_C= [ax.plot(*[ [0+cordvec[0,j,12]/2,(cordvec[0,j,i]+cordvec[0,j,9])/2] for j in range(3) ]) for i in range(10,13)]
        plt.setp(RBO.baxes_C,linestyle='--')
        RBO.ClineCMxAxis = np.hstack([RBO.cordvec[:,:,12]/2,(RBO.cordvec[:,:,10]/2+RBO.cordvec[:,:,12])/2])
        RBO.ClineCMyAxis = np.hstack([RBO.cordvec[:,:,12]/2,(RBO.cordvec[:,:,11]/2+RBO.cordvec[:,:,12])/2])
        RBO.ClineCMzAxis = np.hstack([np.zeros((RBO.N+1,3)),RBO.cordvec[:,:,12]])
        RBO.linesarg = RBO.linesarg + (RBO.baxes_C[0],RBO.ClineCMxAxis,
            RBO.baxes_C[1],RBO.ClineCMyAxis,
            RBO.baxes_C[2],RBO.ClineCMzAxis)

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
    if RBO.DrawOption['A_square'] == True:
        RBO.polylines = [ ax.plot(*[ [cordvec[istep,j,polynum[i]],
                                      cordvec[istep,j,polynum[i+1]]] for j in range(3) ]
                                      ) for i in range(4) ]
        plt.setp(RBO.polylines,lw=5)
        RBO.SquareEdge3 = np.hstack([RBO.cordvec[:,:,3],RBO.cordvec[:,:,4]])
        RBO.SquareEdge4 = np.hstack([RBO.cordvec[:,:,4],RBO.cordvec[:,:,5]])
        RBO.SquareEdge5 = np.hstack([RBO.cordvec[:,:,5],RBO.cordvec[:,:,6]])
        RBO.SquareEdge6 = np.hstack([RBO.cordvec[:,:,6],RBO.cordvec[:,:,3]])
        RBO.linesarg = RBO.linesarg+ (RBO.polylines[0],RBO.SquareEdge3,
                                        RBO.polylines[1],RBO.SquareEdge4,
                                        RBO.polylines[2],RBO.SquareEdge5,
                                        RBO.polylines[3],RBO.SquareEdge6)

    #change previous w(t_i-1) for w(t i-1) comparison
    trailsamplerate = 1    
    if RBO.DrawOption['A_z_axis_trace'] == True:
        zlineAtrace=ax.plot(cordvec[::trailsamplerate,0,2],cordvec[::trailsamplerate,1,2],cordvec[::trailsamplerate,2,2],'b-',
                          linewidth=2,label='A method with $\omega(t_{i+1})$ and B method')
    if RBO.DrawOption['A_Angular Momentum Trace'] == True:
        L_plot = RBO.L_plot
        lineLtrance=ax.plot(L_plot[::trailsamplerate,0],L_plot[::trailsamplerate,1],L_plot[::trailsamplerate,2],'k-',markersize=2)
    if RBO.DrawOption['B_z_axis_trace'] == True:
        Lzppnorm = RBO.Lzppnorm
        HasbunZtrace=ax.plot(Lzppnorm[:,0],Lzppnorm[:,1],Lzppnorm[:,2],'k.',
                     markersize=1,label='$\omega(t_{i+1})$ and Hasbun\'s result')
    if RBO.DrawOption['B_Angular Momentum Trace'] == True:
        L_Bnorm = RBO.L_Bnorm
        L_Btraceline=ax.plot(RBO.L_Bnorm[:,0],RBO.L_Bnorm[:,1],RBO.L_Bnorm[:,2],'k.',
                     markersize=1)
    if RBO.DrawOption['C_z_axis_trace'] == True:
        zlineCtrace=ax.plot(cordvec[::trailsamplerate,0,12],cordvec[::trailsamplerate,1,12],cordvec[::trailsamplerate,2,12],'g.',
                          linewidth=1,label='C method')
                     
#    ax.legend()
    if RBO.DrawOption['A_Angular Momentum Vec'] == True:
        L_plot = RBO.L_plot
        
        RBO.lineL=ax.plot([0,L_plot[istep,0]],
                          [0,L_plot[istep,1]],
                          [0,L_plot[istep,2]],'k-',marker='o')
#        print RBO.lineL
#        fanArrow3d=Arrow3D([0,L_plot[istep,0]],
#                          [0,L_plot[istep,1]],
#                          [0,L_plot[istep,2]],arrowstyle="-|>", color="b")
#        #RBO.lineL=ax.add_artist(fanArrow3d)
#        print fanArrow3d
        RBO.lineLAData = np.hstack([np.zeros((RBO.N+1,3)),RBO.L_plot])
        RBO.linesarg = RBO.linesarg + (RBO.lineL,RBO.lineLAData)        

    if RBO.DrawOption['B_Angular Momentum Vec'] == True:    
        RBO.lineL=ax.plot([0,L_Bnorm[0,0]],
                          [0,L_Bnorm[0,1]],
                          [0,L_Bnorm[0,2]],'r-',marker='o')
        #shape of L_Bnorm is (N,3), one vector less than N+1
        RBO.lineLBData = np.hstack([np.zeros((RBO.N,3)),RBO.L_Bnorm])
        RBO.linesarg = RBO.linesarg + (RBO.lineL,RBO.lineLBData)

    #RBO.HasbunL=ax.plot([0,L_plot[0,0]],
    #               [0,L_plot[0,1]],
     #              [0,L_plot[0,2]],'k-')

    # plot angular velocity vector normalized to w(t0)--
    #print(w[0,2],np.linalg.norm(w_lab[0,:]))
    if RBO.DrawOption['A_Angular Velocity Vec'] == True:
        RBO.w_lab_norm = RBO.w_lab/RBO.w[0,2]
        RBO.lineW=ax.plot([0,RBO.w_lab_norm[istep,0]],
                          [0,RBO.w_lab_norm[istep,1]],
                          [0,RBO.w_lab_norm[istep,2]],'g-',marker='o')
        RBO.w_bnormData = np.hstack([np.zeros((RBO.N+1,3)),RBO.w_lab_norm])
        RBO.linesarg = RBO.linesarg + (RBO.lineW,RBO.w_bnormData)

def rotmat_from_A_2_B(A,B):
    rotunit = np.cross(A, B)
    rotunit = rotunit/np.linalg.norm(rotunit)
    thetarot = np.arccos(np.dot(A/np.linalg.norm(A),B/np.linalg.norm(B)))
    rotmat = CK(thetarot*rotunit)
    return rotmat

def plot_body_space_cone(ax,zn,omegavec,Lvec):
    '''
    Plot the body and space cone. A static one-time plot.
    '''
    # Set up the grid in polar coordinate theta, radius
    cone_theta = np.linspace(0,2*np.pi,90)
    bcone_radius = np.abs(np.linalg.norm(zn))*np.sqrt(1-
    np.square(np.dot(zn,omegavec/np.linalg.norm(omegavec)))
    )
    bcone_len = np.abs(np.linalg.norm(zn))*np.dot(zn,omegavec/np.linalg.norm(omegavec))
    
    cone_r = np.linspace(0,bcone_radius,15)
    c1T, c1R = np.meshgrid(cone_theta, cone_r)
    # Then calculate X, Y, and Z
    c1X = c1R * np.cos(c1T)
    c1Y = c1R * np.sin(c1T)
    c1Z = np.sqrt(c1X**2 + c1Y**2)/bcone_radius*bcone_len
    cmm,cnn = np.shape(c1Z)
    mat_zn2omega = rotmat_from_A_2_B(np.array([0,0,1]),zn)
    for i in range(cmm):
        for j in range(cnn):
            c1X[i,j],c1Y[i,j],c1Z[i,j] = np.dot(mat_zn2omega,
            np.array([c1X[i,j],c1Y[i,j],c1Z[i,j]]))
    # space cone
    scone_radius = np.abs(np.linalg.norm(omegavec))*np.sqrt(1-
    np.square(np.dot(Lvec/np.linalg.norm(Lvec),omegavec/np.linalg.norm(omegavec)))
    )
    scone_len = np.abs(np.linalg.norm(omegavec))*np.dot(Lvec/np.linalg.norm(Lvec),omegavec/np.linalg.norm(omegavec))
    
    cone_r_s = np.linspace(0,scone_radius,15)
    c2T, c2R = np.meshgrid(cone_theta, cone_r_s)
    # Then calculate X, Y, and Z
    c2X = c2R * np.cos(c2T)
    c2Y = c2R * np.sin(c2T)
    c2Z = np.sqrt(c2X**2 + c2Y**2)/scone_radius*scone_len
    cmm2,cnn2 = np.shape(c2Z)
    mat_omega2L = rotmat_from_A_2_B(np.array([0,0,1]),Lvec)
    for i in range(cmm2):
        for j in range(cnn2):
            c2X[i,j],c2Y[i,j],c2Z[i,j] = np.dot(mat_omega2L,
            np.array([c2X[i,j],c2Y[i,j],c2Z[i,j]]))
            
    
    # Set the Z values outside your range to NaNs so they aren't plotted
#    c1Z[c1Z < 0] = np.nan
#    c1Z[c1Z > 2.1] = np.nan
    ax.plot_surface(c1X, c1Y, c1Z,rstride=5, cstride=5,linewidth=0,alpha=0.2,color = 'red')    
    ax.plot_surface(c2X, c2Y, c2Z,rstride=5, cstride=5,linewidth=0,alpha=0.2,color = 'red')    
#    ax.set_zlim(0,1)
    



#def RBanimation(RBO,plt,fig2):
#    #
#    N,baxes,polylines,baxes_hasbun,lineL,lineW,cordvec,L_plot,w_lab_norm = RBO.N,RBO.baxes,RBO.polylines,RBO.baxes_hasbun,RBO.lineL,RBO.lineW,RBO.cordvec,RBO.L_plot,RBO.w_lab_norm
#    import matplotlib.animation as animation
#    # Keep "line_ani =" part, otherwise won't work, timedanimation need an return to clear the figure
#    line_ani = animation.FuncAnimation(fig2, update_line, list(range(1,N,1)),
#                                           fargs=(baxes+polylines+baxes_hasbun,lineL,lineW,cordvec,L_plot,w_lab_norm),
#                                          interval=10, blit=False,repeat=False)
#    #line_ani.save('3Dtop.mp4',writer = 'ffmpeg',fps='24')
#                                          
#    #plt.savefig(r'C:\Documents and Settings\user\My Documents\tony\2014\Xelatexfolder\wti_wtiplus1.pgf')
        
