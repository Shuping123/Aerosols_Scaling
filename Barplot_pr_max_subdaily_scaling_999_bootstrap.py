#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:33:28 2019
Clausius-Claperyon equation scale
@author: shupli
"""


import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#import matplotlib.cm as cm
#from matplotlib.colors import ListedColormap
#from scipy import stats
import pickle
#from collections import OrderedDict
from matplotlib.ticker import FixedFormatter, FixedLocator
#import seaborn as sns
#import pandas as pd
from bootstrap_routines_3D import *

import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)


variablename='precip'
dir2='/home/EAS-aerosol/EAS-hour/'

run =['REF','SUx03','SUx3','BCx03','BCx3']


color =['silver','skyblue','C0','gold','C1']
#color =[c0,c1,c2,c3,c4]

def load(fname):
    filehandler = open(fname, 'rb')
    obj=pickle.load(filehandler)
    filehandler.close()
    return obj

subreg  =['NC', 'SE']
subreg04=xr.open_dataset('/home/EAS-aerosol/COORDINATE/Subregions_EAS004_sample_5x5.nc')
lab =['(a)', '(b)', '(c)']


median = np.zeros (shape=5) ## the median scaling rate of 
def barplot_median (filename):
    fig=plt.figure(figsize=(14,5))
    for kr, reg in enumerate (subreg):
        ax=plt.subplot( 1,2, kr+1) 
        ax = plt.gca()
        #ax.yaxis.grid(True,linestyle='dashed', linewidth=0.5, zorder=0)
        
        for xx, runid in  enumerate (run):
            slope_1hr = load (dir2+'pr_max_C-C_scaling_1hr_'+run[xx]+'_gridpoint_975_99_999.pkl')
            slope_3hr = load (dir2+'pr_max_C-C_scaling_3hr_'+run[xx]+'_gridpoint_975_99_999.pkl')
            slope_6hr = load (dir2+'pr_max_C-C_scaling_6hr_'+run[xx]+'_gridpoint_95_99_999_missing.pkl')
            slope_day = load (dir2+'pr_max_C-C_scaling_day_'+run[xx]+'_gridpoint_95_99_999.pkl')
            
            
            ### 1hr
            scale = np.ma.masked_where (slope_1hr[2,:,:] == -99999 , slope_1hr[2,:,:])
            scale = np.ma.masked_where (subreg04[reg].values !=1, scale)
            #scale = np.ma.masked_where (pval_1hr[2,:,:]  >0.05 , scale)
            scale = scale [~scale.mask]
            scale = np.ravel (scale)
            scale = scale[np.logical_not(np.isnan(scale))]
            
            ###  3hr
            scale1 = np.ma.masked_where (slope_3hr[2,:,:] == -99999 , slope_3hr[2,:,:])
            scale1 = np.ma.masked_where (subreg04[reg].values !=1, scale1)
            #scale1 = np.ma.masked_where (pval_3hr[2,:,:]  >0.05 , scale1)
            scale1 = scale1 [~scale1.mask]
            scale1 = np.ravel (scale1)
            scale1 = scale1[np.logical_not(np.isnan(scale1))]
            
            ###  6hr
            scale2 = np.ma.masked_where (slope_6hr[2,:,:] == -99999 , slope_6hr[2,:,:])
            scale2 = np.ma.masked_where (subreg04[reg].values !=1, scale2)
            #scale2 = np.ma.masked_where (pval_6hr[2,:,:]  >0.05 , scale2)
            scale2 = scale2 [~scale2.mask]
            scale2 = np.ravel (scale2)
            scale2 = scale2[np.logical_not(np.isnan(scale2))]
            
            #### day
            
            scale3 = np.ma.masked_where (slope_day[2,:,:] == -99999 , slope_day[2,:,:])
            scale3 = np.ma.masked_where (subreg04[reg].values !=1, scale3)
            #scale3 = np.ma.masked_where (pval_day[2,:,:]  >0.05 , scale3)
            scale3 = scale3 [~scale3.mask]
            scale3 = np.ravel (scale3)
            scale3 = scale3[np.logical_not(np.isnan(scale3))]
            
            
           
            scale  = scale *100
            scale1 = scale1 *100
            scale2 = scale2 *100
            scale3 = scale3 *100
            
            #print (np.shape(scale), np.shape(scale1), np.shape(scale2), np.shape(scale3))
            
            lower, upper, best = bootci(scale, stat=np.median, nboot=4000, replacement=True, alpha=0.1)
            plt.bar (xx*0.2,best, color=color[xx],edgecolor='black', width =0.2, linewidth=0.3,alpha=0.7)
            plt.plot ([xx*0.2,xx*0.2], [lower, upper],color='black',linewidth=0.6,marker='_',markersize=4)
            #plt.text  (xx*0.22-0.13, upper+0.1, str('%.2f' % round(best,2)),color= 'blue', fontsize=6)
            
            lower, upper, best = bootci(scale1, stat=np.median, nboot=4000, replacement=True, alpha=0.1)
            plt.bar (xx*0.2+1.8,best, color=color[xx],edgecolor='black', width =0.2,linewidth=0.3,alpha=0.7)
            plt.plot ([xx*0.2+1.8,xx*0.2+1.8], [lower, upper],color='black',linewidth=0.6,marker='_',markersize=4)
            
            
            lower, upper, best = bootci(scale2, stat=np.median, nboot=4000, replacement=True, alpha=0.1)
            plt.bar (xx*0.2+3.6,best, color=color[xx],edgecolor= 'black', width =0.2, linewidth=0.3,alpha=0.7)
            plt.plot ([xx*0.2+3.6,xx*0.2+3.6], [lower, upper],color='black',linewidth=0.6,marker='_',markersize=4)
            
            lower, upper, best = bootci(scale3, stat=np.median, nboot=4000, replacement=True, alpha=0.1)
            dot, = plt.bar (xx*0.2+5.4,best, color=color[xx],edgecolor='black', width =0.2, linewidth=0.3,alpha=0.7)
            plt.plot ([xx*0.2+5.4,xx*0.2+5.4], [lower, upper],color='black',linewidth=0.6,marker='_',markersize=4)
            #plt.text  (xx*0.22+5.27, upper+0.1, str('%.2f' % round(best,2)),color= 'blue', fontsize=6)
            
 
        plt.xlim(-0.4,6.6)
        plt.ylim(6,13)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.axvline (x=6.6, linestyle ='dashed',color='gray',linewidth = 0.6)
        plt.axhline(y=7.0, linestyle ='dashed',color='m',linewidth =1.0)
        plt.axhline(y=13.0, linestyle ='dashed',color='gray',linewidth =0.5)
        plt.axhline(y=8.0, linestyle ='dashed',color='gray',linewidth =0.5)
        plt.axhline(y=9.0, linestyle ='dashed',color='gray',linewidth =0.5)
        plt.axhline(y=10.0, linestyle ='dashed',color='gray',linewidth =0.5)
        plt.axhline(y=11.0, linestyle ='dashed',color='gray',linewidth =0.5)
        plt.axhline(y=12.0, linestyle ='dashed',color='gray',linewidth =0.5)
        
        plt.xlabel ('Precipitation timescale', fontsize=12)
        plt.xticks( fontsize=12)
        if kr ==0:
                #plt.ylim(5,11)
                plt.ylabel ('Scaling rate (%/K)', fontsize=12)
                plt.yticks( fontsize=12)
                lgnd = ax.legend([dot,dot,dot,dot,dot],['REF','SU/3', 'SUx3','BC/3','BCx3'],
                                 loc='upper right',fontsize=8,frameon=True,numpoints=1)
                lgnd.legendHandles[0].set_color('silver')
                lgnd.legendHandles[1].set_color('skyblue')
                lgnd.legendHandles[2].set_color('C0')
                lgnd.legendHandles[3].set_color('gold')
                lgnd.legendHandles[4].set_color('C1')
                lgnd.legendHandles[4].set_edgecolor('black')
                lgnd.legendHandles[3].set_edgecolor('black')
                lgnd.legendHandles[2].set_edgecolor('black')
                lgnd.legendHandles[1].set_edgecolor('black')
                lgnd.legendHandles[0].set_edgecolor('black')
                
        else:
            ax.set_yticklabels([''])

        ax.xaxis.set_major_locator(FixedLocator([0.4,2.2,4.0,5.8]))
        ax.xaxis.set_major_formatter(FixedFormatter(['1hr','3hr','6hr','day']))
        plt.title (' '+lab[kr]+' '+reg, loc='left',fontsize=12)
            
    
    fig.subplots_adjust(hspace=0.15,wspace=0.13,top =0.92, bottom=0.13, right =0.96, left=0.07)
    plt.savefig(filename, dpi=500)
        
        
barplot_median ('Barplot_pr_max_subdaily_scaling_NC_SE_bootstrap.png') 


print ('the program done')
