#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 25 22:27:35 2024

@author: chingchen
"""

import pandas as pd
import numpy as np
import numpy.ma as ma
from matplotlib import cm
import matplotlib  as mpl
from scipy.misc import derivative
import matplotlib.pyplot as plt
### setting ###
plt.rcParams["font.family"] = "Helvetica"
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
colors=['#282130','#849DAB','#35838D','#CD5C5C',
        '#97795D','#414F67','#4198B9','#2F4F4F']
header_list = ['Pint','etaref','pCLA','conv','P','zbot','vol',
              'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
labelsize = 30
fontsize=20
bwith = 3
fig,(ax22,ax44) = plt.subplots(2,2,figsize=(24,14))
# ----------------------------------- viscosity 1--------------------------------
model_list = [
              'Europa-tidal1_eta1.0d12_P0.0-3.0TW_0.0wt%-NH3',
              'Europa-tidal1_eta3.2d12_P0.0-3.0TW_0.0wt%-NH3',
              'Europa-tidal1_eta1.0d13_P0.0-3.0TW_0.0wt%-NH3',
              'Europa-tidal1_eta3.2d13_P0.0-3.0TW_0.0wt%-NH3',
              'Europa-tidal1_eta1.0d14_P0.0-3.0TW_0.0wt%-NH3',
              'Europa-tidal1_eta3.2d14_P0.0-3.0TW_0.0wt%-NH3',
              'Europa-tidal1_eta1.0d15_P0.0-3.0TW_0.0wt%-NH3',]
mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
vis = cm.get_cmap('cubehelix',len(model_list))
newcolors = vis(np.linspace(0, 1, len(model_list)+2))
ax=ax22[0]
ax2=ax44[0]
ax3=ax22[1]
ax4=ax44[1]
for i, model in enumerate(model_list):
   data = pd.read_csv(workpath+model+'_2_thermal-evolution.dat',
       header=None,names=header_list,delim_whitespace=True)[3:].reset_index(drop=True)
   data = data.replace('D','e',regex=True).astype(float) # data convert to float
   Pint=data.Pint/1e12
   x = Pint
   ax.plot(x[data.conv==1],data.zbot[data.conv==1],c=newcolors[i],lw=5)   
   ax.plot(x[data.conv==0],data.zbot[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
   if len(x[data.conv==0])!=0 and len(x[data.conv==1])!=0:
       xx=np.array(x[data.conv==0])[0]
       ax.vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                   ymax=data.zbot[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)  
   ax2.plot(x[data.conv==1],data.dlid[data.conv==1],c=newcolors[i],lw=5)   
   ax2.plot(x[data.conv==0],data.dlid[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
#----------------------------------- viscosity 2--------------------------------

model_list = [
              'Europa-tidal1_eta1.0d12_P0.0-3.0TW_3.0wt%-NH3',
              'Europa-tidal1_eta3.2d12_P0.0-3.0TW_3.0wt%-NH3',
              'Europa-tidal1_eta1.0d13_P0.0-3.0TW_3.0wt%-NH3',
              'Europa-tidal1_eta3.2d13_P0.0-3.0TW_3.0wt%-NH3',
              'Europa-tidal1_eta1.0d14_P0.0-3.0TW_3.0wt%-NH3',
              'Europa-tidal1_eta3.2d14_P0.0-3.0TW_3.0wt%-NH3',
              'Europa-tidal1_eta1.0d15_P0.0-3.0TW_3.0wt%-NH3',
              ]
vis = cm.get_cmap('cubehelix',len(model_list))
newcolors = vis(np.linspace(0, 1, len(model_list)+2))
for i, model in enumerate(model_list):
    data = pd.read_csv(workpath+model+'_2_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[3:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    Pint=data.Pint/1e12
    x = Pint
    ax3.plot(x[data.conv==1],data.zbot[data.conv==1],c=newcolors[i],lw=5)
    ax3.plot(x[data.conv==0],data.zbot[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
    if len(x[data.conv==0])!=0 and len(x[data.conv==1])!=0:
        xx=np.array(x[data.conv==0])[0]
        ax3.vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                    ymax=data.zbot[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)
    ax4.plot(x[data.conv==1],data.dlid[data.conv==1],c=newcolors[i],lw=5)   
    ax4.plot(x[data.conv==0],data.dlid[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
#--------------------------------- figure setting ------------------------------

for aa in [ax2,ax4]:
   aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
   aa.set_xlim(0,3)
   aa.set_ylim(60,0)
   aa.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize)
   aa.grid()
   aa.set_xlabel('P$_{tide}$',fontsize=labelsize)
   for axis in ['top','bottom','left','right']:
       aa.spines[axis].set_linewidth(bwith)
for aa in [ax,ax3]:
   aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
   aa.set_xlim(0,3)
   aa.set_ylim(161,0)
   aa.set_ylabel('Ice layer thickness (km)',fontsize = labelsize)
   aa.grid()
   for axis in ['top','bottom','left','right']:
       aa.spines[axis].set_linewidth(bwith)

# fig.savefig('/Users/chingchen/Desktop/StagYY_Works/paper_europa_ice_shell/figure3.pdf')