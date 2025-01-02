#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 17:08:14 2024

@author: chingchen
"""

import pandas as pd
import numpy as np
import numpy.ma as ma
import matplotlib
# from matplotlib import cm
import matplotlib  as mpl
# from scipy.misc import derivative
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

labelsize = 30
fontsize = 30
bwith = 3
tmax=150
tmin=0
smax=20
smin=0
### PATH ###
path = '/Users/chingchen/Desktop/data/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/thermal_evolution_v2/'
workpath = '/Users/chingchen/Desktop/StagYY_Works/Thermal_evolution_20240826/'
modelpath = '/Users/chingchen/Desktop/model/'
figpath = '/Users/chingchen/Desktop/figure/'
colors=['#282130','#3CB371','#4682B4','#CD5C5C','#97795D','#414F67','#4198B9','#3CB371']
header_list = ['Pint','etaref','pCLA',
               'zbot_avg','zbot_min','zbot_max']
###------------------------------ figure two meshed ------------------------------
fig2,(ax5,ax6) = plt.subplots(1,2,figsize=(26,11))
plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
# ---------------------------------------- figure 3 --------------------------------
model = 'Europa_eta1.0e13-1.0e15_P0.0-1.2TW_1.5%-NH3'
header_list = ['Pint','etaref','pCLA',
               'zbot_avg','zbot_min','zbot_max']
data_ori = pd.read_csv(workpath+model+'_Hvar_minmax_stat-ice-thickness.dat',
                   header=None,sep='\s+')[2:].dropna(axis=1)
data_ori.columns = header_list
data = data_ori.replace('D','e',regex=True).astype(float) # data convert to float

header_list = ['Pint','etaref','pCLA',
               'zlid_avg','zlid_min','zlid_max']
data2 = pd.read_csv(workpath+model+'_Hvar_minmax_stat-stagnant-lid.dat',
                   header=None,sep='\s+')[2:].dropna(axis=1)
data2.columns = header_list
data2 = data2.replace('D','e',regex=True).astype(float) # data convert to float

power_mesh = np.ones((len(data.Pint[data.Pint==0]),len(data.Pint[data.etaref==1e13]) ))
vis_mesh = np.ones((len(data.Pint[data.Pint==0]),len(data.Pint[data.etaref==1e13]) ))
ice_thickness = np.ones((len(data.Pint[data.Pint==0]),len(data.Pint[data.etaref==1e13]) ))
power_mesh2 = np.ones((len(data2.Pint[data2.Pint==0]),len(data2.Pint[data2.etaref==1e13]) ))
vis_mesh2 = np.ones((len(data2.Pint[data2.Pint==0]),len(data2.Pint[data2.etaref==1e13]) ))
stg_thickness = np.ones((len(data2.Pint[data2.Pint==0]),len(data2.Pint[data2.etaref==1e13]) ))
for uu in range(len(power_mesh[0])):
    vis_mesh[:,uu] = np.array(data.etaref[data.Pint==0])
    for mm in range(len(power_mesh)):
        power_mesh[mm,:] = np.array(data.Pint[data.etaref==1e13])
        pp = power_mesh[mm,uu]
        vv = vis_mesh[mm,uu]
        ice_thickness[mm,uu]=data.zbot_min[data.Pint==pp][data.etaref==vv]
for uu in range(len(power_mesh2[0])):
    vis_mesh2[:,uu] = np.array(data2.etaref[data2.Pint==0])
    for mm in range(len(power_mesh2)):
        power_mesh2[mm,:] = np.array(data2.Pint[data2.etaref==1e13])
        pp = power_mesh2[mm,uu]
        vv = vis_mesh2[mm,uu]        
        stg_thickness[mm,uu]=data2.zlid_min[data2.Pint==pp][data2.etaref==vv]
qqq=ax5.pcolormesh(power_mesh/1e12,vis_mesh,ice_thickness,cmap='magma_r',vmin=tmin,vmax=tmax)    
kkk=ax6.pcolormesh(power_mesh2/1e12,vis_mesh2,stg_thickness,cmap='magma_r',vmin=smin,vmax=smax)  

ax5.contour(power_mesh/1e12,vis_mesh,ice_thickness, levels=[10], colors='white',linestyle='dashed')
ax5.contour(power_mesh/1e12,vis_mesh,ice_thickness, levels=[20], colors='white')
ax5.contour(power_mesh/1e12,vis_mesh,ice_thickness, levels=[47], colors='white',linewidths=3)
ax5.contour(power_mesh/1e12,vis_mesh,ice_thickness, levels=[24.2], colors='white',linewidths=5)

ax6.contour(power_mesh/1e12,vis_mesh2,stg_thickness, levels=[6,8], colors='white',linewidths=3)
ax6.contour(power_mesh/1e12,vis_mesh2,stg_thickness, levels=[10,15], colors='white',linewidths=5)
#  ------------------------------ figure setting ------------------------------
ax5.set_xlabel('Tidal power (TW)',fontsize=labelsize)
ax6.set_xlabel('Tidal power (TW)',fontsize=labelsize)
# fig.tight_layout()
for aa in [ax5,ax6]:
    aa.set_ylim(1e13,1e15)
    aa.set_ylabel('Viscosity',fontsize = labelsize)
    # aa.yaxis.set_major_locator(ticker.FixedLocator([1e13,3.2e13,5e13,1e14,3.2e14,1e15]))
    # aa.minorticks_on()
    aa.set_yscale("log") 
    # aa.tick_params(which='minor', length=5, width=2, direction='in')
    aa.tick_params(labelsize=labelsize,width=3,length=10,)
                    # right=True, top=True,direction='in',pad=15)
    aa.set_xlim(0,1.2)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
    aa.yaxis.set_major_locator(ticker.FixedLocator([1e13,3.2e13,5e13,1e14,3.2e14,1e15]))
    # aa.set_yticks([1e13,5e13,1e14,1e15])
# colorbar figure 5
axc1  = fig2.add_axes([0.128,0.016,0.35,0.035])
# norm = mpl.colors.LogNorm(vmin=10,vmax=151)
norm = mpl.colors.Normalize(vmin=tmin,vmax=tmax)
cb1  = mpl.colorbar.ColorbarBase(axc1,cmap='magma_r',norm=norm,orientation='horizontal')
cb1.set_label('Ice shell thickness',fontsize=labelsize+5)
cb1.ax.tick_params(axis='x', labelsize=labelsize)
# cb1.ax.yaxis.set_label_position('right')
# colorbar figure 6
axc1  = fig2.add_axes([0.55,0.016,0.35,0.035])
norm = mpl.colors.Normalize(vmin=smin, vmax=smax)
cb1  = mpl.colorbar.ColorbarBase(axc1,cmap='magma_r',norm=norm,orientation='horizontal')
cb1.set_label('Stagnant lid thickness',fontsize=labelsize+5)
cb1.ax.tick_params(axis='x', labelsize=labelsize-2)
# cb1.ax.yaxis.set_label_position('bottom')

ax5.set_title('(a) min ice shell thickness',fontsize=labelsize)
ax6.set_title('(b) min stagnant lid thickness',fontsize=labelsize)
# fig2.savefig('/Users/chingchen/Desktop/StagYY_Works/paper_europa_ice_shell/figure8_v3.pdf')