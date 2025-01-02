#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 17:44:24 2024

@author: chingchen
"""

import pandas as pd
import numpy as np
import numpy.ma as ma
from matplotlib import cm
import matplotlib  as mpl
from scipy.misc import derivative
import matplotlib.pyplot as plt

labelsize = 20
bwith = 3
plt.rcParams["font.family"] = "Helvetica"
path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
modelpath = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/'

colors=['#282130','#849DAB','#35838D','#CD5C5C','#97795D','#414F67','#4198B9','#2F4F4F']

melting_point_TW=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4])
melting_point_frac=np.array([0.75,0.4,0.3,0.25,0.2,0.16,0.14,0.12,0.11,0.10,0.09,0.08,0.07])
melting_point_depth=[172,243,364,375,381,320,420,370,380,350,380,320,370]

fig2,(axqq) = plt.subplots(1,1,figsize=(8,6))
axqq.scatter(melting_point_TW,100-melting_point_frac*100,color='#35838D')
axqq.set_ylim(0,100)
axqq.set_xlim(0,1.5)
axqq.set_xlabel('P$_{tide}$ (TW)',fontsize=labelsize)
axqq.set_ylabel('Fraction of tidal heating within ice shell (%)',fontsize=labelsize-2)
for aa in [axqq]:
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=True,top=True,direction='in',pad=15)
    aa.grid()
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
# fig2.savefig('/Users/chingchen/Desktop/StagYY_Works/figure/figureS4.pdf')