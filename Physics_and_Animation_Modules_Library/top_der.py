#top_der.py: returns the derivatives for the symmetric top problem
import numpy as np
def top_der(t,w,flag,I,Is,ws,tau0):
#w(1):phi, w(2):phi_dot, w(3):theta, w(4):theta_dot, w(5):psi
#main program produces w(6):psi_dot
    print(shape(w))
    ders=[w(2),(Is*ws-2*I*w(2)*np.cos(w(3)))*w(4)/(I*np.sin(w(3))),
          w(4),(tau0-(Is*ws-I*w(2)*np.cos(w(3)))*w(2))*np.sin(w(3))/I,
          ws-w(2)*np.cos(w(3))];
    return ders
