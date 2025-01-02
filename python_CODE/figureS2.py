#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 18:55:19 2024

@author: chingchen
"""

import pandas as pd
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

labelsize = 20
bwith = 3

### PATH ###
path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
modelpath = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/'
colors=['#282130','#3CB371','#4682B4','#CD5C5C','#97795D','#414F67','#4198B9','#3CB371']
header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
               'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']



fig1,(aa1,aa2) = plt.subplots(2,2,figsize=(17,11))
ax1 = aa1[0]
ax2 = aa1[1]
ax3 = aa2[0]
ax4 = aa2[1]
    
label_list=['0.6 TW','0.8 TW','1.0 TW','1.2 TW']
# ---------------------------------------- figure 3 dash 1.0 TW --------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.6TW_1.5wt%_D5.0km-NH3_Hvar_2', # power
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.8TW_1.5wt%_D5.0km-NH3_Hvar_2',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3_Hvar_2',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.2TW_1.5wt%_D5.0km-NH3_Hvar_2',]
   
for i, model in enumerate(model_list):
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.Ftop, mask = data.conv).mask
    zbot_cond = ma.array(data.Ftop, mask = mask_cond)
    zbot_conv = ma.array(data.Ftop, mask = mask_conv)
    ax1.plot(x,zbot_conv,color=colors[i],lw=3,label=label_list[i])
    ax1.plot(x,zbot_cond,color='orange') 


label_list=['10$^{13}$','10$^{13.5}$','10$^{13.75}$','10$^{14}$']
# ---------------------------------------- figure 3 dash 1.0 TW --------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d13_P1.0TW_1.5wt%_D5.0km-NH3_Hvar_2',# viscosity
             'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d13_P1.0TW_1.5wt%_D5.0km-NH3_Hvar_2', 
             'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P1.0TW_1.5wt%_D5.0km-NH3_Hvar_2',
             'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3_Hvar_2',
             # 'Europa-tidal5_period0.14Gyr_emx10%_eta3.2d14_P1.0TW_1.5wt%_D5.0km-NH3',
             ]   
for i, model in enumerate(model_list):
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.Ftop, mask = data.conv).mask
    zbot_cond = ma.array(data.Ftop, mask = mask_cond)
    zbot_conv = ma.array(data.Ftop, mask = mask_conv)
    ax2.plot(x,zbot_conv,color=colors[i],lw=3,label=label_list[i])
    ax2.plot(x,zbot_cond,color='orange') 

label_list=['vol = 0.%','vol = 1.5%','vol = 3.0%']
# ---------------------------------------- figure 4 dashed 1.0 TW--------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_0.0wt%_D5.0km-NH3',# composition
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_3.0wt%_D5.0km-NH3',
              ]
for i, model in enumerate(model_list):
    data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.Ftop, mask = data.conv).mask
    zbot_cond = ma.array(data.Ftop, mask = mask_cond)
    zbot_conv = ma.array(data.Ftop, mask = mask_conv)
    ax3.plot(x,zbot_conv,color=colors[i],lw=2,label=label_list[i])
    ax3.plot(x,zbot_cond,color='orange')   
label_list = ['5%','20%','50%']
# ---------------------------------------- figure 3 dash 1.0 TW --------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx5%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              # 'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx20%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx50%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              ] 
for i, model in enumerate(model_list):
    data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.Ftop, mask = data.conv).mask
    zbot_cond = ma.array(data.Ftop, mask = mask_cond)
    zbot_conv = ma.array(data.Ftop, mask = mask_conv)
    ax4.plot(x,zbot_conv,color=colors[i],lw=3,label=label_list[i])
    ax4.plot(x,zbot_cond,color='orange') 
#  ------------------------------ figure setting ------------------------------
ax3.set_xlabel('Time (Gyr)',fontsize=labelsize)
ax4.set_xlabel('Time (Gyr)',fontsize=labelsize)
for aa in [ax1,ax2,ax3,ax4]:
    aa.legend(fontsize=labelsize)
    aa.set_ylim(0,100)
    aa.set_ylabel('Surface heat flux',fontsize = labelsize)
    aa.minorticks_on()
    aa.tick_params(which='minor', length=5, width=2, direction='in')
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
    aa.set_xlim(0,4.55)
    aa.grid()
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
# fig1.savefig('/Users/chingchen/Desktop/StagYY_Works/paper_europa_ice_shell/figureS5_v7.pdf')