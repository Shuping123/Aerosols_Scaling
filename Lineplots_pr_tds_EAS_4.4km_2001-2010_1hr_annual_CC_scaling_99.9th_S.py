#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:33:28 2019
Clausius-Claperyon equation scale
@author: shupli
"""

#from netCDF4 import Dataset
#import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.cm as cm
import xarray as xr
#from matplotlib.colors import ListedColormap
from scipy import stats
from collections import OrderedDict

import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)



dir2='/home/scripts/EAS-aerosol/subdaily/CCscale/data/'
        
        
linestyles_dict = OrderedDict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1)))])

color=['blue','green']

subreg=['NC','SE']
label =['(a)', '(b)']
label1 =['(c)', '(d)']
linestyles =['solid','dashed',linestyles_dict['densely dotted']]

timescale=['1hr','3hr','6hr','day']

run =['REF','SUx03','SUx3']
run1 =['REF','BCx03','BCx3']

###color =['green','C3','C0','m']

fig=plt.figure(figsize=(12,10))
plot_lines =[]
def CC_scaling (filename):   
  for kk, reg in enumerate (subreg):
    ax=plt.subplot (2,2, kk+1)
    for xx, runid in enumerate (run):
        scale1 = np.loadtxt(dir2+'COSMO_EAS04_1hr_pr_tds_CC_scaling_bin_'+run[xx]+'_'+reg+'_max.txt')
        scale2 = np.loadtxt(dir2+'COSMO_EAS04_3hr_pr_tds_CC_scaling_bin_'+run[xx]+'_'+reg+'_max.txt')
        scale3 = np.loadtxt(dir2+'COSMO_EAS04_6hr_pr_tds_CC_scaling_bin_'+run[xx]+'_'+reg+'_max.txt')
        scale4 = np.loadtxt(dir2+'COSMO_EAS04_day_pr_tds_CC_scaling_bin_'+run[xx]+'_'+reg+'_max.txt')
        
        line1, = plt.plot( scale1[:,0], scale1[:,2], color='green',linewidth=1.5, linestyle=linestyles[xx]) ##99.9
        line2, = plt.plot( scale2[:,0], scale2[:,2], color='C3',linewidth=1.5, linestyle=linestyles[xx])
        line3, = plt.plot( scale3[:,0], scale3[:,3], color='C0',linewidth=1.5, linestyle=linestyles[xx])
        line4, = plt.plot( scale4[:,0], scale4[:,4], color='m',linewidth=1.5, linestyle=linestyles[xx])
    
    plot_lines.append([line1, line2, line3])
    
    plt.plot(np.arange(0,41)-2,8 * (1+0.07)**(np.arange(0,41)-2), color='gray',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41)-2,3 * (1+0.07)**(np.arange(0,41)-2), color='gray',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41)-2,50 * (1+0.07)**(np.arange(0,41)-2), color='gray',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41)-2,20 * (1+0.07)**(np.arange(0,41)-2), color='gray',linestyle='dotted' ,linewidth=0.8)
    #plt.plot(np.arange(0,41),1.2* (1+0.07)**np.arange(0,41), color='black',linestyle='dotted' )
    
    plt.plot(np.arange(0,41),8 * (1+0.14)**np.arange(0,41), color='C1',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41),3 * (1+0.14)**np.arange(0,41), color='C1',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41),50 * (1+0.14)**np.arange(0,41), color='C1',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41),20 * (1+0.14)**np.arange(0,41), color='C1',linestyle='dotted' ,linewidth=0.8)
    #plt.plot(np.arange(0,41),1.2 * (1+0.14)**np.arange(0,41), color='gray',linestyle='dotted' )
    
    plt.yscale('log')
    plt.ylim(2,300)
    

    plt.xlim(-2,28)
    plt.title (label[kk], loc = 'left',fontsize=13)
    plt.title (reg, loc = 'center',fontsize=14,weight='bold')
    plt.xticks(fontsize = 13)
    ax.set_xticklabels([''])
    plt.yticks(fontsize = 13)
    
    lgnd = plt.legend(plot_lines[0], ["REF", "SU/3", "SUx3"], loc='upper left',ncol=1,fontsize=9,frameon=False)
    lgnd.legendHandles[0].set_linestyle ('solid')
    lgnd.legendHandles[1].set_linestyle ('dashed')
    lgnd.legendHandles[2].set_linestyle (linestyles_dict['densely dotted'])
    lgnd.legendHandles[0].set_color ('black')
    lgnd.legendHandles[1].set_color ('black')
    lgnd.legendHandles[2].set_color ('black')

    if kk ==0:
        plt.ylabel('Precipitation extemes (mm)',fontsize=13)

        
        lgnd1 = plt.legend([line1,line2,line3,line4 ],['1hr','3hr','6hr','day'] , loc='lower right',ncol=1,mode="expand",fontsize=10,bbox_to_anchor=(0.55, 0.05,0.2,0.46))
        plt.gca().add_artist(lgnd)
        lgnd1.legendHandles[0].set_linestyle ('solid')
        lgnd1.legendHandles[1].set_linestyle ('solid')
        lgnd1.legendHandles[2].set_linestyle ('solid')
        lgnd1.legendHandles[3].set_linestyle ('solid')
    
    plt.subplot (2,2, kk+3)
    for xx, runid in enumerate (run1):
        scale1 = np.loadtxt(dir2+'COSMO_EAS04_1hr_pr_tds_CC_scaling_bin_'+runid+'_'+reg+'_max.txt')
        scale2 = np.loadtxt(dir2+'COSMO_EAS04_3hr_pr_tds_CC_scaling_bin_'+runid+'_'+reg+'_max.txt')
        scale3 = np.loadtxt(dir2+'COSMO_EAS04_6hr_pr_tds_CC_scaling_bin_'+runid+'_'+reg+'_max.txt')
        scale4 = np.loadtxt(dir2+'COSMO_EAS04_day_pr_tds_CC_scaling_bin_'+runid+'_'+reg+'_max.txt')
        
        line1, = plt.plot( scale1[:,0], scale1[:,2], color='green',linewidth=1.5, linestyle=linestyles[xx]) ##99.9
        line2, = plt.plot( scale2[:,0], scale2[:,2], color='C3',linewidth=1.5, linestyle=linestyles[xx])
        line3, = plt.plot( scale3[:,0], scale3[:,3], color='C0',linewidth=1.5, linestyle=linestyles[xx])
        line4, = plt.plot( scale4[:,0], scale4[:,4], color='m',linewidth=1.5, linestyle=linestyles[xx])

    
    plot_lines.append([line1, line2, line3])
    
    plt.plot(np.arange(0,41)-2,8 * (1+0.07)**(np.arange(0,41)-2), color='gray',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41),3 * (1+0.07)**(np.arange(0,41)-2), color='gray',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41),50 * (1+0.07)**(np.arange(0,41)-2), color='gray',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41),20 * (1+0.07)**(np.arange(0,41)-2), color='gray',linestyle='dotted' ,linewidth=0.8)
    #plt.plot(np.arange(0,41),1.2* (1+0.07)**np.arange(0,41), color='black',linestyle='dotted' )
    
    plt.plot(np.arange(0,41)-2,8 * (1+0.14)**(np.arange(0,41)-2), color='C1',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41)-2,3 * (1+0.14)**(np.arange(0,41)-2), color='C1',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41)-2,50 * (1+0.14)**(np.arange(0,41)-2), color='C1',linestyle='dotted' ,linewidth=0.8)
    plt.plot(np.arange(0,41)-2,20 * (1+0.14)**(np.arange(0,41)-2), color='C1',linestyle='dotted' ,linewidth=0.8)
    #plt.plot(np.arange(0,41),1.2 * (1+0.14)**np.arange(0,41), color='gray',linestyle='dotted' )
    
    plt.yscale('log')
    plt.ylim(3,300)
    #if kk ==0 or kk ==2:
    #    plt.ylim(1,100)
    #else:
    #    plt.ylim(1,100)

    plt.xlim(-2,28)
    plt.xlabel('Dew point temperature ($^\circ$C)',fontsize=13)
    plt.title (label1[kk], loc = 'left',fontsize=13)
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
    
    lgnd = plt.legend(plot_lines[0], ["REF", "BC/3", "BCx3"], loc='upper left',ncol=1,fontsize=8,frameon=False)
    lgnd.legendHandles[0].set_linestyle ('solid')
    lgnd.legendHandles[1].set_linestyle ('dashed')
    lgnd.legendHandles[2].set_linestyle (linestyles_dict['densely dotted'])
    lgnd.legendHandles[0].set_color ('black')
    lgnd.legendHandles[1].set_color ('black')
    lgnd.legendHandles[2].set_color ('black')

    if kk ==0:
        plt.ylabel('Precipitation extremes (mm)',fontsize=13)

        



    fig.subplots_adjust(hspace=0.12,wspace=0.18,top =0.94, bottom=0.08, right =0.95, left=0.08)
    plt.savefig(filename, dpi =500)

CC_scaling ('Lineplot_pr_tds_4.4km_1hr_scaling_S.png')
print ('the program done')


