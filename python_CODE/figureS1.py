#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 13:36:48 2024

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
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
colors=['#282130','#849DAB','#35838D','#CD5C5C',
        '#97795D','#414F67','#4198B9','#2F4F4F']
header_list = ['Pint','etaref','pCLA','conv','P','zbot','vol',
              'Tbot','Tm','TH2O','Fbot','Ftop','dlid','T_core']
labelsize = 30
bwith = 3
fig,(ax1,ax2) = plt.subplots(2,2,figsize=(24,14))
ax = ax1[0]
ax3 = ax1[1]
axi = ax2[0]
ax4 = ax2[1]
#--------------------------------------- vol 1 ------------------------------------
model_list = [
               'Europa-tidal1_eta3.2d13_P0.0-3.0TW_0.0wt%-NH3',
                'Europa-tidal1_eta3.2d13_P0.0-3.0TW_0.5wt%-NH3',
                'Europa-tidal1_eta3.2d13_P0.0-3.0TW_1.0wt%-NH3',
              'Europa-tidal1_eta3.2d13_P0.0-3.0TW_1.5wt%-NH3',
              'Europa-tidal1_eta3.2d13_P0.0-3.0TW_2.0wt%-NH3',
              'Europa-tidal1_eta3.2d13_P0.0-3.0TW_2.5wt%-NH3',
              'Europa-tidal1_eta3.2d13_P0.0-3.0TW_3.0wt%-NH3',]
rainbow = cm.get_cmap('winter_r',len(model_list))
newcolors = rainbow(np.linspace(0, 1, len(model_list)))
label_list=['vol = 0.0%','vol = 0.5%','vol = 1.0%','vol = 1.5%','vol = 2.0%','vol = 2.5%','vol = 3.0%']
for i, model in enumerate(model_list):
   data = pd.read_csv(workpath+model+'_2_thermal-evolution.dat',
       header=None,names=header_list,delim_whitespace=True)[3:].reset_index(drop=True)
   data = data.replace('D','e',regex=True).astype(float) # data convert to float
   Pint=data.Pint/1e12
   x = Pint
   ax.plot(x[data.conv==1],data.zbot[data.conv==1],c=newcolors[i],label = label_list[i],lw=5)
   axi.plot(x[data.conv==1],data.dlid[data.conv==1],color=newcolors[i],lw=5)
   
   ax.plot(x[data.conv==0],data.zbot[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
   axi.plot(x[data.conv==0],data.dlid[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
   
   if len(x[data.conv==0])!=0 and len(x[data.conv==1])!=0:
       xx=np.array(x[data.conv==0])[0]
       
       ax.vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                  ymax=data.zbot[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)
       axi.vlines(x=xx, ymin=data.dlid[data.conv==0].iloc[0], 
                   ymax=data.dlid[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)

#--------------------------------------- vol 2 ------------------------------------
model_list = ['Europa-tidal1_eta1.0d14_P0.0-3.0TW_0.0wt%-NH3',
              'Europa-tidal1_eta1.0d14_P0.0-3.0TW_0.5wt%-NH3',
              'Europa-tidal1_eta1.0d14_P0.0-3.0TW_1.0wt%-NH3',
              'Europa-tidal1_eta1.0d14_P0.0-3.0TW_1.5wt%-NH3',
              'Europa-tidal1_eta1.0d14_P0.0-3.0TW_2.0wt%-NH3',
              'Europa-tidal1_eta1.0d14_P0.0-3.0TW_2.5wt%-NH3',
              'Europa-tidal1_eta1.0d14_P0.0-3.0TW_3.0wt%-NH3',]
# rainbow = cm.get_cmap('winter',len(model_list))
# newcolors = rainbow(np.linspace(0, 1, len(model_list)))
label_list=['vol = 0.0%','vol = 0.5%','vol = 1.0%','vol = 1.5%','vol = 2.0%','vol = 2.5%','vol = 3.0%']
for i, model in enumerate(model_list):
   data = pd.read_csv(workpath+model+'_2_thermal-evolution.dat',
       header=None,names=header_list,delim_whitespace=True)[3:].reset_index(drop=True)
   data = data.replace('D','e',regex=True).astype(float) # data convert to float
   Pint=data.Pint/1e12
   x = Pint
   ax3.plot(x[data.conv==1],data.zbot[data.conv==1],c=newcolors[i],label = label_list[i],lw=5)
   ax4.plot(x[data.conv==1],data.dlid[data.conv==1],color=newcolors[i],lw=5)
   
   ax3.plot(x[data.conv==0],data.zbot[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
   ax4.plot(x[data.conv==0],data.dlid[data.conv==0],color=newcolors[i],linestyle='dashed',lw=3)
   
   if len(x[data.conv==0])!=0 and len(x[data.conv==1])!=0:
       xx=np.array(x[data.conv==0])[0]
       
       ax3.vlines(x=x[data.conv==0].iloc[0], ymin=data.zbot[data.conv==0].iloc[0], 
                  ymax=data.zbot[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)
       ax4.vlines(x=xx, ymin=data.dlid[data.conv==0].iloc[0], 
                   ymax=data.dlid[data.conv==1].iloc[-1], colors=newcolors[i], ls='--', lw=3,)

#--------------------------------- figure setting ------------------------------
# ax.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
# ax2[0].set_ylabel('stagnant lid thickness (km)',fontsize = labelsize-5)
# ax3.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
axi.set_xlabel('P$_{tide}$',fontsize=labelsize)
ax4.set_xlabel('P$_{tide}$',fontsize=labelsize)

for aa in [ax,ax3,]:
   aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
   aa.set_xlim(0,3)
   aa.set_ylim(161,0)
   aa.grid()
   aa.set_ylabel('Ice layer thickness (km)',fontsize = labelsize)
   for axis in ['top','bottom','left','right']:
       aa.spines[axis].set_linewidth(bwith)
for aa in [axi,ax4]:
   aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
   aa.set_xlim(0,3)
   aa.set_ylim(20,0)
   aa.grid()
   aa.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize)
   for axis in ['top','bottom','left','right']:
       aa.spines[axis].set_linewidth(bwith)
       
# colorbar
# axc1  = fig.add_axes([0.95,0.516,0.03,0.363])
# norm = mpl.colors.LogNorm(vmin=3.2e12,vmax=3.2e15)
# cb1  = mpl.colorbar.ColorbarBase(axc1,cmap=vis,norm=norm,orientation='vertical')
# cb1.set_label('viscosity',fontsize=labelsize+5)
# cb1.ax.tick_params(axis='y', labelsize=labelsize-2)
# cb1.ax.yaxis.set_label_position('right')
# colorbar
axc2  = fig.add_axes([0.95,0.12,0.03,0.763])
norm = mpl.colors.Normalize(vmin=0,vmax=3)
cb1  = mpl.colorbar.ColorbarBase(axc2,cmap=rainbow,norm=norm,orientation='vertical')
cb1.set_label('vol %',fontsize=labelsize+5)
cb1.ax.tick_params(axis='y', labelsize=labelsize-2)
cb1.ax.yaxis.set_label_position('right')
fig.savefig('/Users/chingchen/Desktop/StagYY_Works/paper_europa_ice_shell/figureS2_v5.pdf')