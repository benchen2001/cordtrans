import numpy as np


def uniformprecesscal(M,L,g,I3,I1,I2,w3,theta0,psi0):
    '''
    M = 1
    L= 0.04
    g = 9.8
    I3 = 0.0008
    I1 = 0.002
    I2 = 0.002
    w3 = 20*2*np.pi
    theta0 = 0.01   #check sign here
    '''
    a = I1*np.cos(theta0)
    b = -I3*w3
    c = M*g*L

    d = b**2-4*a*c # discriminant

    if d < 0:
        print('This equation has no real solution')
    elif d == 0:
        x = (-b+np.sqrt(b**2-4*a*c))/2/a
        print('This equation has one solutions: ' + str(x))
    else:
        x1 = (-b+np.sqrt(b**2-4*a*c))/2/a
        x2 = (-b-np.sqrt(b**2-4*a*c))/2/a
        print('This equation has two solutions: '+ str(x1) + 'and' + str(x2))

#    psi1,psi2=w3-x1*np.cos(theta0),w3-x2*np.cos(theta0)

    w11,w12 =x1*np.sin(theta0)*np.sin(psi0),x2*np.sin(theta0)*np.sin(psi0)
    w21,w22 =x1*np.sin(theta0)*np.cos(psi0),x2*np.sin(theta0)*np.cos(psi0)
    print('first w1,w2 = '+str(w11)+ ' and ' +str(w21))
    print('second w1,w2 = '+str(w12)+ ' and ' +str(w22))
    return w11,w21,w12,w22