# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 12:09:54 2014

@author: user
"""

import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(5,figsize=(6.4,8))
ax = fig.add_subplot(111)
#ax.set_axis_off()
ax.tick_params(labelbottom='off', labeltop='off', labelleft='off', labelright='off')


print(fig.get_dpi())
print(fig.get_size_inches())