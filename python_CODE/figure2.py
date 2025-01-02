#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:17:28 2024

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

### PATH ###
path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
modelpath = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/'
colors=['#282130','#849DAB','#35838D','#CD5C5C','#97795D','#414F67','#4198B9','#2F4F4F']
colors=['#B22222','#849DAB','#4B0082','#282130','#D2691E']

header_list = ['time_Gyr','Prad','Ptidal','Fcore','Pint','Hint','conv',
               'melt','P','zbot','%vol','Tbot','Tm','Fbot','Ftop','dlid','T_core']

# 
fig,(aa1,aa2) = plt.subplots(2,1,figsize=(11,18))
ax=aa1

# ---------------------------------------- figure --------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P0.6TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx5%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3'
              ]

label_list = ['5%','10%','20%','30%','50%']
label_list = ['161 km','138 km','10%','15%','20%']

min_zbot = []
max_zbot = []
amplitude_zbot=[]

i = 0
model = model_list[0]
data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
    header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
data = data.replace('D','e',regex=True).astype(float) # data convert to float
x = data.time_Gyr
mask_cond = data.conv
mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
zbot_cond = ma.array(data.Ptidal, mask = mask_cond)
zbot_conv = ma.array(data.Ptidal, mask = mask_conv)

model = 'Europa-tidal1_eta1.0d14_P1.0TW_1.5wt%-NH3'
data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
    header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
data = data.replace('D','e',regex=True).astype(float) # data convert to float
x = data.time_Gyr

ax.plot(x,data.Prad,color=colors[3],lw=3)
Radio = data.Fcore*4*np.pi**3/3/1000
aa2.plot(x,data.Fcore,color=colors[i],lw=3)

model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P0.6TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta5.6d13_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx5%_eta5.6d13_P1.0TW_1.5wt%_D5.0km-NH3'
              ]

for i, model in enumerate(model_list):
    i = i
    data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.Pint, mask = mask_cond)
    zbot_conv = ma.array(data.Pint, mask = mask_conv)
    ax.plot(x,data.Pint+data.Prad,color=colors[i],label=label_list[i],lw=3)   
    print(np.max(data.Pint))


#  ------------------------------ figure setting ------------------------------
ax.set_ylim(0,2.5)
aa2.set_ylim(0,20)
ax.set_ylabel('Total Power (TW)',fontsize = labelsize)
aa2.set_ylabel('Surface heat flux ',fontsize = labelsize)
aa2.set_xlabel('Time (Gyr)',fontsize=labelsize)

for aa in [ax,aa2]:
    aa.tick_params(which='minor', length=5, width=2, direction='in')
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=False, top=True,direction='in',pad=15)
    aa.minorticks_on()
    aa.grid()
for aa in [ax,aa2]:
    aa.set_xlim(0,4.55)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
# fig.savefig('/Users/chingchen/Desktop/StagYY_Works/paper_europa_ice_shell/figure2.pdf')