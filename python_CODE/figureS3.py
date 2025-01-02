#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:43:47 2024

@author: chingchen
"""

import pandas as pd
import numpy as np
import numpy.ma as ma
from matplotlib import cm
import matplotlib  as mpl
from scipy.misc import derivative
import matplotlib.pyplot as plt
from scipy.signal import find_peaks



labelsize = 20
bwith = 3

### PATH ###
path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
modelpath = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/'
colors=['#282130','#849DAB','#CD5C5C','#35838D',]

header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
               'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']


fig,(aa1,aa2,aa3) = plt.subplots(3,1,figsize=(9,21))

ax=aa1
ax3=aa2
axftop=aa3
# ---------------------------------------- figure --------------------------------
model_list = ['Europa-tidal1_eta1.0d14_P1.0TW_1.5wt%-NH3_core0.00_Hvar_2', # power
               'Europa-tidal1_eta1.0d14_P1.0TW_1.5wt%-NH3_core0.05_Hvar_2',
               'Europa-tidal1_eta1.0d14_P1.0TW_1.5wt%-NH3_core0.10_Hvar_2',
              # 'Europa-tidal1_eta5.6d13_P0.8TW_1.5wt%-NH3_core0.15_Hvar_2'
              ]
label_list=['100%','95%','90%']


min_zbot = []
max_zbot = []
amplitude_zbot=[]
for i, model in enumerate(model_list):
    i = i
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.dlid, mask = mask_cond)
    zbot_conv = ma.array(data.dlid, mask = mask_conv)
    if i < 3:
        ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    else:
        ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=1)
    #------------------------------------------------------------------------------------------
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.Tm, mask = mask_cond)
    zbot_conv = ma.array(data.Tm, mask = mask_conv)
    if i < 3:
        ax3.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    else:
        ax3.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=1)
    #------------------------------------------------------------------------------------------
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.Ftop, mask = mask_cond)
    zbot_conv = ma.array(data.Ftop, mask = mask_conv)
    if i < 3:
        axftop.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    else:
        axftop.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=1)
# ---------------------------------------- figure 2 --------------------------------
model_list = ['Europa-tidal1_eta1.0d14_P0.6TW_1.5wt%-NH3_core0.00_Hvar_2', # power
              'Europa-tidal1_eta1.0d14_P0.6TW_1.5wt%-NH3_core0.05_Hvar_2',
              'Europa-tidal1_eta1.0d14_P0.6TW_1.5wt%-NH3_core0.10_Hvar_2',
              'Europa-tidal1_eta1.0d14_P0.6TW_1.5wt%-NH3_core0.15_Hvar_2'
              ]
label_list=['100%','95%','90%','85%']
min_zbot = []
max_zbot = []
amplitude_zbot=[]
for i, model in enumerate(model_list):
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.dlid, mask = mask_cond)
    zbot_conv = ma.array(data.dlid, mask = mask_conv)
    
    if i < 3:
        ax.plot(x,zbot_conv,color=colors[i],lw=3)
    else:
        ax.plot(x,zbot_conv,color=colors[i],lw=1)
    #------------------------------------------------------------------------------------------
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.Tm, mask = mask_cond)
    zbot_conv = ma.array(data.Tm, mask = mask_conv)
    if i < 3:
        ax3.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    else:
        ax3.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=1)
    #------------------------------------------------------------------------------------------
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.Ftop, mask = mask_cond)
    zbot_conv = ma.array(data.Ftop, mask = mask_conv)
    if i < 3:
        axftop.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    else:
        axftop.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=1)
  # ------------------------------ figure setting ------------------------------
ax.set_ylim(20,0)
ax3.set_ylim(240,280)
axftop.set_ylim(20,60)
ax.legend(fontsize=labelsize)
ax.set_ylabel('stagnant thickness (km)',fontsize = labelsize)
ax3.set_ylabel('T$_m$ (K)',fontsize = labelsize)
axftop.set_xlabel('Time (Gyr)',fontsize=labelsize)
axftop.set_ylabel('surface heat flux',fontsize = labelsize)


for aa in [ax,ax3,axftop]:
    aa.minorticks_on()
    aa.tick_params(which='minor', length=5, width=2, direction='in')
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=5)
    aa.set_xlim(0,4.55)
    aa.grid()
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
# fig.savefig('/Users/chingchen/Desktop/StagYY_Works/paper_europa_ice_shell/figureS3.pdf')