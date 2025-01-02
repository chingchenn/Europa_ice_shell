#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 14:46:46 2024

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
fig,(aa1,aa2) = plt.subplots(2,2,figsize=(22,14))
ax=aa1[0]
ax2=aa1[1]
ax3=aa2[0]
ax4=aa2[1]
# ---------------------------------------- figure --------------------------------
model_list = ['Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_0.0wt%_D5.0km-NH3',# composition
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_1.5wt%_D5.0km-NH3',
              'Europa-tidal5_period0.14Gyr_emx10%_eta1.0d14_P1.0TW_3.0wt%_D5.0km-NH3',
              ]

label_list=['0%','1.5%','3.0%']

min_zbot = []
max_zbot = []
amplitude_zbot=[]
for i, model in enumerate(model_list):
    i = i
    data = pd.read_csv(workpath+model+'_Hvar_2_thermal-evolution.dat',
        header=None,names=header_list,delim_whitespace=True)[2:].reset_index(drop=True)
    data = data.replace('D','e',regex=True).astype(float) # data convert to float
    x = data.time_Gyr
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.zbot, mask = mask_cond)
    zbot_conv = ma.array(data.zbot, mask = mask_conv)
    ax.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    
    #------------------------------------------------------------------------------------------
    peaks, _ = find_peaks(zbot_conv)
    mins, _  = find_peaks(zbot_conv*-1)
    if len(zbot_conv[peaks])>20:
        amplitude_zbot.append(np.median(zbot_conv[peaks][10:19]-zbot_conv.data[mins][10:19]))
    else:
        print('---- HELP ----')      
    mask_cond = data.conv[data.time_Gyr>3.5]
    mask_conv = ~ma.array(data.zbot[data.time_Gyr>3.5], mask = data.conv[data.time_Gyr>3.5]).mask
    
    zbot_cond = ma.array(data.zbot[data.time_Gyr>3.5], mask = mask_cond)
    zbot_conv = ma.array(data.zbot[data.time_Gyr>3.5], mask = mask_conv)
    pp1 = round(np.min(zbot_conv),1)
    pp2 = round(np.max(zbot_conv),1)
    min_zbot.append(pp1)
    max_zbot.append(pp2)
    #--------------------------------------------- dlid ---------------------------------------------
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.dlid, mask = mask_cond)
    zbot_conv = ma.array(data.dlid, mask = mask_conv)
    ax3.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    
    #---------------------------------------------- Tm ---------------------------------------
    mask_cond = data.conv
    mask_conv = ~ma.array(data.zbot, mask = data.conv).mask
    zbot_cond = ma.array(data.Tm, mask = mask_cond)
    zbot_conv = ma.array(data.Tm, mask = mask_conv)
    ax4.plot(x,zbot_conv,color=colors[i],label=label_list[i],lw=3)
    
    
playerA = np.array(max_zbot)
playerB = np.array(min_zbot)
hat_graph(ax2, label_list, [playerA, playerB], ['Player A', 'Player B'])
ax2.set_ylim(161,0)
ax2.set_ylabel('Ice layer thickness (km)',fontsize = labelsize)
ax2.tick_params(labelsize=labelsize,width=3,length=10,right=False, top=True,direction='in',pad=10)


#  ------------------------------ figure setting ------------------------------
ax.legend(fontsize=labelsize)
ax.set_ylim(161,0)
ax2.set_ylim(161,0)
ax3.set_ylim(20,5)
ax4.set_ylim(240,275)

ax.set_ylabel('Ice layer thickness (km)',fontsize = labelsize)
ax2.set_xlabel('wt% of ammonia ',fontsize=labelsize)
ax3.set_ylabel('Stagnant lid thickness (km)',fontsize = labelsize)
ax4.set_ylabel('T$_m$ (K)',fontsize = labelsize)
for aa in [ax,ax3,ax4]:
    aa.minorticks_on()
    aa.set_xlim(0,4.55)
    aa.tick_params(which='minor', length=5, width=2, direction='in')
    aa.tick_params(labelsize=labelsize,width=3,length=10,right=True, top=True,direction='in',pad=15)
    aa.grid()
    aa.set_xlabel('Time (Gyr)',fontsize=labelsize)
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
for aa in [ax2]:
    aa.minorticks_on()
    for axis in ['top','bottom','left','right']:
        aa.spines[axis].set_linewidth(bwith)
    aa.grid()
    

# fig.savefig('/Users/chingchen/Desktop/StagYY_Works/paper_europa_ice_shell/figure7_v3.pdf')