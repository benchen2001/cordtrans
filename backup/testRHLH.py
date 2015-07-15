# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 16:24:15 2015

@author: user
"""
import numpy as np
import RGCordTransV13 as RG

rotmat = RG.CK([-pi/2,0,0])

print np.dot(rotmat,[0,1,0])