#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 15:17:22 2024

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

def hat_graph(ax, xlabels, values, group_labels):
    """
    Create a hat graph.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The Axes to plot into.
    xlabels : list of str
        The category names to be displayed on the x-axis.
    values : (M, N) array-like
        The data values.
        Rows are the groups (len(group_labels) == M).
        Columns are the categories (len(xlabels) == N).
    group_labels : list of str
        The group labels displayed in the legend.
    """

    def label_bars(heights, rects):
        """Attach a text label on top of each bar."""
        for height, rect in zip(heights, rects):
            ax.annotate(f'{height}',
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 2),  # 4 points vertical offset.
                        textcoords='offset points',
                        ha='center', va='bottom',fontsize=14)

    values = np.asarray(values)
    x = np.arange(values.shape[1])
    ax.set_xticks(x, labels=xlabels)
    spacing = 0.3  # spacing between hat groups
    width = (1 - spacing) / values.shape[0]
    heights0 = values[0]
    for i, (heights, group_label) in enumerate(zip(values, group_labels)):
        style = {'fill': False} if i == 0 else {'edgecolor': 'black'}
        rects = ax.bar(x, heights - heights0,
                       width, bottom=heights0, label=group_label, **style)
        label_bars(heights, rects)
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

# 
fig,(aa1) = plt.subplots(1,1,figsize=(9,7))
ax=aa1
# ---------------------------------------- figure --------------------------------



model_list = ['Europa-tidal1_eta1.0d14_P0.6TW_1.5wt%-NH3_core0.00_Hvar_2', # power
              'Europa-tidal1_eta1.0d14_P0.6TW_1.5wt%-NH3_core0.05_Hvar_2',
              'Europa-tidal1_eta1.0d14_P0.6TW_1.5wt%-NH3_core0.10_Hvar_2',
              'Europa-tidal1_eta1.0d14_P0.6TW_1.5wt%-NH3_core0.15_Hvar_2'
              ]



model_list2 = [
              'Europa-tidal1_eta1.0d14_P1.0TW_1.5wt%-NH3_core0.00_Hvar_2', # power
               'Europa-tidal1_eta1.0d14_P1.0TW_1.5wt%-NH3_core0.05_Hvar_2',
               'Europa-tidal1_eta1.0d14_P1.0TW_1.5wt%-NH3_core0.10_Hvar_2',
              # 'Europa-tidal1_eta5.6d13_P0.8TW_1.5wt%-NH3_core0.15_Hvar_2'
              ]

label_list=['100%','95%','90%','85%']
label_list2=['100%','95%','90%']

for i, model in enumerate(model_list):
    i = i
    data = pd.read_csv(workpath+model+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.zbot, mask = mask_cond)
    zbot_conv = ma.array(data.zbot, mask = mask_conv)
    
    ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    ax.plot(x,zbot_cond,color=colors[i],linestyle='dashed')
for i, model2 in enumerate(model_list2):
    i = i    
    data = pd.read_csv(workpath+model2+'_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.zbot, mask = mask_cond)
    zbot_conv = ma.array(data.zbot, mask = mask_conv)
    
    ax.plot(x,zbot_conv,color=colors[i],lw=2,linestyle='dashed')
    if i==0:
        ax.plot(x,zbot_conv,color=colors[i],lw=3,linestyle='dashed')
    
#  ------------------------------ figure setting ------------------------------
ax.set_ylim(161,0)
ax.set_xlabel('Time (Gyr)',fontsize=labelsize)
ax.legend(fontsize=labelsize)
for aa in [ax]:
    aa.set_ylabel('Ice layer thickness (km)',fontsize = labelsize)
    aa.minorticks_on()
    aa.tick_params(which='minor', length=5, width=2, direction='in')
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
    aa.set_xlim(0,4.55)
    aa.grid()
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)

# fig.savefig('/Users/chingchen/Desktop/StagYY_Works/paper_europa_ice_shell/figure9.pdf')