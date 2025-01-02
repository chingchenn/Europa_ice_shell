#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 25 22:25:15 2024

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

fig3 = 0
fig3_S = 1
if fig3:
    fig,(ax3,ax4) = plt.subplots(1,2,figsize=(17,7))

    ##---------------------------------------- figure 3 1.0 TW --------------------------------
    model_list = [
                   # 'Europa-tidal1_eta1.0d14_P1.0TW_0.0wt%-NH3',
                   'Europa-tidal1_eta1.0d14_P1.0TW_1.5wt%-NH3',
                  # 'Europa-tidal1_eta1.0d14_P1.0TW_3.0wt%-NH3',
                  ]
    label_list=['vol = 0.%','vol = 1.5%','vol = 3.0%']
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        mask_cond = data.conv
        mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
        zbot_cond = ma.array(data.zbot, mask = mask_cond)
        zbot_conv = ma.array(data.zbot, mask = mask_conv)
        ax3.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
        ax3.plot(x,zbot_cond,color='darkred')  
    # ---------------------------------------- figure 3 dash --------------------------------

        mask_cond = data.conv
        mask_conv = ~ma.array(data.dlid, mask = data.conv).mask
        zbot_cond = ma.array(data.dlid, mask = mask_cond)
        zbot_conv = ma.array(data.dlid, mask = mask_conv)
        ax3.plot(x,zbot_conv,color=colors[i],lw=3,linestyle='dashed')
        ax3.plot(x,zbot_cond,color='darkred') 
        print(zbot_conv)
    # ---------------------------------------- figure 4 --------------------------------
    model_list = [
                  # 'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.2TW_0.0wt%_D5.0km-NH3',# composition
                  'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.2TW_1.5wt%_D5.0km-NH3',
                  # 'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.2TW_3.0wt%_D5.0km-NH3',
                  ]

    # ---------------------------------------- figure 4 1.0 TW --------------------------------

    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        mask_cond = data.conv
        mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
        zbot_cond = ma.array(data.zbot, mask = mask_cond)
        zbot_conv = ma.array(data.zbot, mask = mask_conv)
        ax4.plot(x,zbot_conv,color=colors[i],lw=3)
        ax4.plot(x,zbot_cond,color='darkred')   
        print(np.max(data.zbot[data.time_Gyr>2.5]),np.mean(data.zbot[data.time_Gyr>2.5]),np.min(data.zbot[data.time_Gyr>2.5]))
    
    # # ---------------------------------------- figure 4 dashed --------------------------------

        x = data.time_Gyr
        mask_cond = data.conv
        mask_conv = ~ma.array(data.dlid, mask = data.conv).mask
        zbot_cond = ma.array(data.dlid, mask = mask_cond)
        zbot_conv = ma.array(data.dlid, mask = mask_conv)
        ax4.plot(x,zbot_conv,color=colors[i],linestyle='dashed',lw=1)
        ax4.plot(x,zbot_cond,color='darkred',linestyle='dashed')   
        print(np.max(data.zbot[data.time_Gyr>2.5]),np.mean(data.zbot[data.time_Gyr>2.5]),np.min(data.zbot[data.time_Gyr>2.5]))
    #  ------------------------------ figure setting ------------------------------
    ax3.legend(fontsize=labelsize)
    ax3.set_xlabel('Time (Gyr)',fontsize=labelsize)
    ax4.set_xlabel('Time (Gyr)',fontsize=labelsize)
    for aa in [ax3,ax4]:
        aa.set_ylim(61,0)
        aa.set_ylabel('ice layer thickness (km)',fontsize = labelsize)
        aa.minorticks_on()
        aa.tick_params(which='minor', length=5, width=2, direction='in')
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        aa.set_xlim(0,4.55)
        aa.grid()
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    
    ax3.set_title('(a) viscosity = 10$^{14}$ with constant total power = 1.0 TW',fontsize=labelsize-5)
    ax4.set_title('(b) viscosity = 10$^{14}$ with varying total power = 1.0 TW ',fontsize=labelsize-5)
    # fig.savefig('/Users/chingchen/Desktop/StagYY_Works/paper_europa_ice_shell/figure3_v4.pdf')
if fig3_S:
    fig,(axa,axb) = plt.subplots(2,2,figsize=(17,11))
    ax1=axa[0]
    ax2=axa[1]
    ax3=axb[0]
    ax4=axb[1]
    ##---------------------------------------- figure 1 zbot --------------------------------
    model_list = [
                 # 'Europa-tidal1_eta1.0d14_P1.0TW_0.0wt%-NH3',
                  'Europa-tidal1_eta1.0d14_P1.0TW_1.5wt%-NH3',
                  # 'Europa-tidal1_eta1.0d14_P1.0TW_3.0wt%-NH3',
                  ]
    label_list=['vol = 0.%','vol = 1.5%','vol = 3.0%']
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        mask_cond = data.conv
        mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
        zbot_cond = ma.array(data.zbot, mask = mask_cond)
        zbot_conv = ma.array(data.zbot, mask = mask_conv)
        ax1.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
        ax1.plot(x,zbot_cond,color='darkred')  
        print(np.max(data.zbot[data.time_Gyr>2.5]),np.mean(data.zbot[data.time_Gyr>2.5]),np.min(data.zbot[data.time_Gyr>2.5]))
    # ---------------------------------------- figure 3 dlid --------------------------------
        mask_cond = data.conv
        mask_conv = ~ma.array(data.dlid, mask = data.conv).mask
        zbot_cond = ma.array(data.dlid, mask = mask_cond)
        zbot_conv = ma.array(data.dlid, mask = mask_conv)
        ax3.plot(x,zbot_conv,color=colors[i],lw=3)
        ax3.plot(x,zbot_cond,color='darkred') 
    # ---------------------------------------- figure 3 Tm --------------------------------
        mask_conv = ~ma.array(data.Tm, mask = data.conv).mask
        zbot_cond = ma.array(data.Tm, mask = mask_cond)
        zbot_conv = ma.array(data.Tm, mask = mask_conv)
        ax2.plot(x,zbot_conv,color=colors[i],lw=3)
        ax2.plot(x,zbot_cond,color='darkred')   
    # # ---------------------------------------- figure 4 Ftop --------------------------------
        mask_conv = ~ma.array(data.Ftop, mask = data.conv).mask
        zbot_cond = ma.array(data.Ftop, mask = mask_cond)
        zbot_conv = ma.array(data.Ftop, mask = mask_conv)
        ax4.plot(x,zbot_conv,color=colors[i],lw=3)
        ax4.plot(x,zbot_cond,color='darkred')  
    # ---------------------------------------- figure 2 zbot --------------------------------
    model_list = [
                  # 'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_0.0wt%_D5.0km-NH3',# composition
                  'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
                  # 'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_3.0wt%_D5.0km-NH3',
                  ]
    for i, model in enumerate(model_list):
        data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
            header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
        data = data.replace('D','e',regex=True).astype(float) # data convert to float
        x = data.time_Gyr
        mask_cond = data.conv
        mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
        zbot_cond = ma.array(data.zbot, mask = mask_cond)
        zbot_conv = ma.array(data.zbot, mask = mask_conv)
        ax1.plot(x,zbot_conv,color='gray',lw=3)
        ax1.plot(x,zbot_cond,color='darkred')   
    # # ---------------------------------------- figure 4 dild --------------------------------
        mask_conv = ~ma.array(data.dlid, mask = data.conv).mask
        zbot_cond = ma.array(data.dlid, mask = mask_cond)
        zbot_conv = ma.array(data.dlid, mask = mask_conv)
        ax3.plot(x,zbot_conv,color='gray',lw=3)
        ax3.plot(x,zbot_cond,color='darkred')   
    # ---------------------------------------- figure 4 Tm --------------------------------
        mask_conv = ~ma.array(data.Tm, mask = data.conv).mask
        zbot_cond = ma.array(data.Tm, mask = mask_cond)
        zbot_conv = ma.array(data.Tm, mask = mask_conv)
        ax2.plot(x,zbot_conv,color='gray',lw=3)
        ax2.plot(x,zbot_cond,color='darkred')   
    # # ---------------------------------------- figure 4 Ftop --------------------------------
        mask_conv = ~ma.array(data.Ftop, mask = data.conv).mask
        zbot_cond = ma.array(data.Ftop, mask = mask_cond)
        zbot_conv = ma.array(data.Ftop, mask = mask_conv)
        ax4.plot(x,zbot_conv,color='gray',lw=3)
        ax4.plot(x,zbot_cond,color='darkred')   
    #  ------------------------------ figure setting ------------------------------
    # ax1.legend(fontsize=labelsize-2)
    ax1.set_ylim(70,0)
    ax2.set_ylim(240,275)
    ax3.set_ylim(20,0)
    ax4.set_ylim(0,100)
    ax1.set_ylabel('Ice layer thickness (km)',fontsize = labelsize)
    ax2.set_ylabel('T$_m$',fontsize = labelsize)
    ax3.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize)
    ax4.set_ylabel('Surface heat flux',fontsize = labelsize)
    
    for aa in [ax1,ax2]:
        aa.minorticks_on()
        aa.tick_params(which='minor', length=5, width=2, direction='in')
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        aa.set_xlim(0,4.55)
        aa.grid()
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    for aa in [ax3,ax4]:
        aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
        aa.set_xlim(0,4.55)
        aa.grid()
        aa.set_xlabel('Time (Gyr)',fontsize=labelsize)
        for axis in ['top','bottom','left','right']:
            aa.spines[axis].set_linewidth(bwith)
    
    # ax1.set_title('viscosity = 10$^{14}$ with constant total power = 1.0 TW',fontsize=labelsize-2)
    # ax2.set_title('viscosity = 10$^{14}$ with varying total power = 1.0 TW ',fontsize=labelsize-2)
    # fig.savefig('/Users/chingchen/Desktop/StagYY_Works/paper_europa_ice_shell/figure4.pdf')