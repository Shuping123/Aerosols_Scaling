#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:33:28 2019

@author: shupli
"""

#from netCDF4 import Dataset
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.legend_handler import HandlerLine2D
import numpy as np
#from search_station_in_model import lookupNearest, lookup_radius
#import dypy.small_tools as ds
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)


dir1='/home/scripts/EAS-aerosol/1hr/MODEL/'

##### validation for daily precipitation
nd = 87648
nsta = 106 


data1 = np.loadtxt (dir1+'Stations_info_selected_from_166stations.txt')


### read observational hourly precipitation
obssta = np.zeros (shape=(nd,nsta))
for ista in np.arange (0,nsta):
    ### read observations for 106 stations
    data = np.loadtxt ('/home/OBS_1hr/SURF_CLI_CHN_PRE_FTM_H24-QC-'+str(int(data1[ista,0]))+'.TXT')
    
    obssta [:,ista] = data[:,7:31].reshape(nd)  ### daily precipitation no NAN values
    if np.min(obssta[:,ista]) == np.max(obssta[:,ista]):
        print (np.max(obssta[:,ista]),data[0,ista])
    
    
        
modsta04 = np.loadtxt (dir1 +'COSMO_EAS004_1hr_precipitation_106station_2001-2010.txt')
modsta11 = np.loadtxt (dir1 +'COSMO_EAS11_1hr_precipitation_106station_2001-2010.txt')

### shift Beijing Time to UTC

### UTC Time: 2001.01.02.(00-01)UTC - 2010.12.30.(23-24)UTC (hours: nd-48)
### BJ Time: 2001.01.02.(08-09)BJS - 2010.12.31.(07-08)BJS

modsta04_UTC = modsta04[24:nd-24,:]
modsta11_UTC = modsta11[24:nd-24,:]
obssta_UTC = obssta[36:nd-12,:]

#modsta04_3h = np.zeros (shape=(int((nd-48)/3),nsta))
#modsta04_6h = np.zeros (shape=(int((nd-48)/6),nsta))

#obs_mask =  np.ma.masked_where(obssta_UTC== 32766,obssta_UTC)

obssta_3h = np.zeros (shape=(int((nd-48)/3),nsta))
obssta_6h = np.zeros (shape=(int((nd-48)/6),nsta))

### 3-hourly precipitation
modsta04_3h = (np.concatenate([sum(modsta04_UTC[i:i+3,:]) for i in range(0, len(modsta04_UTC[:,0]), 3)])).reshape(int((nd-48)/3),106)
modsta04_6h = (np.concatenate([sum(modsta04_UTC[i:i+6,:]) for i in range(0, len(modsta04_UTC[:,0]), 6)])).reshape(int((nd-48)/6),106)

modsta11_3h = (np.concatenate([sum(modsta11_UTC[i:i+3,:]) for i in range(0, len(modsta11_UTC[:,0]), 3)])).reshape(int((nd-48)/3),106)
modsta11_6h = (np.concatenate([sum(modsta11_UTC[i:i+6,:]) for i in range(0, len(modsta11_UTC[:,0]), 6)])).reshape(int((nd-48)/6),106)


for ista in np.arange (0,nsta):
    print (ista)
    ### 3/6-hourly precipitation model 
        
    #modsta04_3h[:,ista] = [sum(modsta04_UTC[i:i+3,ista]) for i in range(0, len(modsta04_UTC[:,0]), 3)]
    #modsta04_6h[:,ista] = [sum(modsta04_UTC[i:i+6,ista]) for i in range(0, len(modsta04_UTC[:,0]), 6)]
    
    ### 3/6-hourly observational data
    for i in np.arange (0,int((nd-48)/3)):
        if np.min(obssta_UTC[i*3:(i+1)*3,ista]) == 32766 or np.max(obssta_UTC[i*3:(i+1)*3,ista]) == 32766:
            obssta_3h[i,ista] = 32766
        else:
            obssta_3h[i,ista] = np.sum(obssta_UTC[i*3:(i+1)*3,ista])
    for i in np.arange (0,int((nd-48)/6)):
        if np.min(obssta_UTC[i*6:(i+1)*6,ista]) == 32766 or np.max(obssta_UTC[i*6:(i+1)*6,ista]) == 32766:
            obssta_6h[i,ista] = 32766
        else:
            obssta_6h[i,ista] = np.sum(obssta_UTC[i*6:(i+1)*6,ista])
            

pctl =[99.0,99.9,99.99,99.999]
pctl1 =[95.0,99.0,99.9,99.99]

obs_percent = np.zeros (shape = (nsta,len(pctl)))  ### 4: 99%, 99.9%, 99.99%, 99.999%
mod04_percent= np.zeros (shape =(nsta,len(pctl)))
mod11_percent= np.zeros (shape =(nsta,len(pctl)))

obs_3hpctl = np.zeros (shape = (nsta,len(pctl)))
obs_6hpctl = np.zeros (shape = (nsta,len(pctl)))

mod04_3hpctl = np.zeros (shape = (nsta,len(pctl)))
mod04_6hpctl = np.zeros (shape = (nsta,len(pctl)))
mod11_3hpctl = np.zeros (shape = (nsta,len(pctl)))
mod11_6hpctl = np.zeros (shape = (nsta,len(pctl)))


#### mask

for ista in np.arange (0,nsta):
    mod04 = np.ma.masked_where(obssta_UTC[:,ista]==32766,modsta04_UTC[:,ista])
    mod04_3h = np.ma.masked_where(obssta_3h[:,ista]==32766,modsta04_3h[:,ista])
    mod04_6h = np.ma.masked_where(obssta_6h[:,ista]==32766,modsta04_6h[:,ista])
    
    mod11 = np.ma.masked_where(obssta_UTC[:,ista]==32766,modsta11_UTC[:,ista])
    mod11_3h = np.ma.masked_where(obssta_3h[:,ista]==32766,modsta11_3h[:,ista])
    mod11_6h = np.ma.masked_where(obssta_6h[:,ista]==32766,modsta11_6h[:,ista])
    
    
    obs = np.ma.masked_where(obssta_UTC[:,ista]== 32766,obssta_UTC[:,ista])
    obs_3h = np.ma.masked_where(obssta_3h[:,ista]==32766,obssta_3h[:,ista])
    obs_6h = np.ma.masked_where(obssta_6h[:,ista]==32766,obssta_6h[:,ista])
    
    
    mod04=mod04[~mod04.mask] 
    mod04_3h = mod04_3h[~mod04_3h.mask]
    mod04_6h = mod04_6h[~mod04_6h.mask]
    
    mod11=mod11[~mod11.mask] 
    mod11_3h = mod11_3h[~mod11_3h.mask]
    mod11_6h = mod11_6h[~mod11_6h.mask]
    
    obs=obs[~obs.mask]
    obs_3h = obs_3h[~obs_3h.mask]
    obs_6h = obs_6h[~obs_6h.mask]
    obs = obs *0.1
    obs_3h = obs_3h *0.1
    obs_6h = obs_6h *0.1
    
    for kk in np.arange (0, len(pctl)):
        obs_percent [ista,kk] = np.nanpercentile (obs,q=pctl[kk])
        obs_3hpctl [ista,kk] = np.nanpercentile (obs_3h,q=pctl[kk])
        obs_6hpctl [ista,kk] = np.nanpercentile (obs_6h,q=pctl1[kk])
    
        mod04_percent [ista,kk] = np.nanpercentile (mod04,q=pctl[kk])
        mod04_3hpctl[ista,kk] = np.nanpercentile (mod04_3h,q=pctl[kk])
        mod04_6hpctl[ista,kk] = np.nanpercentile (mod04_6h,q=pctl1[kk])
        
        mod11_percent [ista,kk] = np.nanpercentile (mod11,q=pctl[kk])
        mod11_3hpctl[ista,kk] = np.nanpercentile (mod11_3h,q=pctl[kk])
        mod11_6hpctl[ista,kk] = np.nanpercentile (mod11_6h,q=pctl1[kk])
   

obs_pctl = np.zeros (shape =(nsta, len(pctl),3)) ## 3 -> 1hr, 3hr, and 6hr
mod11_pctl = np.zeros (shape =(nsta, len(pctl),3))
mod04_pctl = np.zeros (shape =(nsta, len(pctl),3))


obs_pctl[:,:,0] = obs_percent ; obs_pctl[:,:,1] = obs_3hpctl ; obs_pctl[:,:,2] = obs_6hpctl
mod11_pctl[:,:,0] = mod11_percent ; mod11_pctl[:,:,1] = mod11_3hpctl; mod11_pctl[:,:,2] = mod11_6hpctl
mod04_pctl[:,:,0] = mod04_percent ; mod04_pctl[:,:,1] = mod04_3hpctl; mod04_pctl[:,:,2] = mod04_6hpctl

obs_daypctl = np.loadtxt ('/home/scripts/EAS-aerosol/1hr/MODEL/OBS_day_precipitation_106station_2001-2010_percentiles.txt')
mod11_daypctl = np.loadtxt ('/home/scripts/EAS-aerosol/1hr/MODEL/COSMO_EAS11_day_precipitation_106station_2001-2010_percentiles.txt')
mod04_daypctl = np.loadtxt ('/home/scripts/EAS-aerosol/1hr/MODEL/COSMO_EAS04_day_precipitation_106station_2001-2010_percentiles.txt')

    
print (obs_percent[:,2])
print (mod04_percent[:,2])
#color =['green','red','blue','m']
color =['green','C3','C0','m']
tab = ['Hourly', '3-Hourly','6-Hourly', 'Daily']

label = ['(a)','(b)','(c)','(d)']
label1 = ['(e)','(f)','(g)','(h)']


def scatter_percentile (filename):
    fig, axs = plt.subplots(figsize=(10,10))
    
    for ii in np.arange (0,3):
        ax = plt.subplot (2,2,1+ii)
        if ii != 2:
            plt.scatter (obs_pctl[:,1,ii],mod11_pctl[:,1,ii] , color='green',marker ='o',s=28)
            plt.scatter (obs_pctl[:,1,ii],mod04_pctl[:,1,ii] , color='blue',marker ='D',s=20)
            r2 = np.corrcoef(obs_pctl[:,1,ii],mod11_pctl[:,1,ii])[0,1]**2
            bias = np.nanmean(mod11_pctl[:,1,ii]-obs_pctl[:,1,ii])
            r2_04 = np.corrcoef(obs_pctl[:,1,ii],mod04_pctl[:,1,ii])[0,1]**2
            bias_04 = np.nanmean (mod04_pctl[:,1,ii] -obs_pctl[:,1,ii])
        else:
            plt.scatter (obs_pctl[:,2,ii],mod11_pctl[:,2,ii] , color='green',marker ='o',s=28)
            plt.scatter (obs_pctl[:,2,ii],mod04_pctl[:,2,ii] , color='blue',marker ='D',s=20)
            r2 = np.corrcoef(obs_pctl[:,2,ii],mod11_pctl[:,2,ii])[0,1]**2
            bias = np.nanmean(mod11_pctl[:,2,ii]-obs_pctl[:,2,ii])
            r2_04 = np.corrcoef(obs_pctl[:,2,ii],mod04_pctl[:,2,ii])[0,1]**2
            bias_04 = np.nanmean (mod04_pctl[:,2,ii] -obs_pctl[:,2,ii])
            
        ax.text (0.06,0.80,'BIAS: '+str('%.2f' % round(bias,2))+' mm',color= 'green',fontsize=12,transform=ax.transAxes)
        ax.text (0.06,0.74,'R$^2$: '+str('%.2f' % round(r2,2)),color= 'green',fontsize=12,transform=ax.transAxes)

        #ax1.text (0.65,0.18,'CTRL04',color= 'blue',fontsize=12,transform=ax1.transAxes)
        ax.text (0.65,0.12,'BIAS: '+str('%.2f' % round(bias_04,2))+' mm',color= 'blue',fontsize=12,transform=ax.transAxes)
        ax.text (0.65,0.06,'R$^2$: '+str('%.2f' % round(r2_04,2)),color= 'blue',fontsize=12,transform=ax.transAxes)
        
        if ii ==0:
            plt.xlim (0.5,80)
            plt.ylim (0.5,80)
            ax.legend (['REF_0.11','REF_0.04'], loc='upper left',fontsize=12)
        elif ii ==1:
            plt.xlim (1,100)
            plt.ylim (1,100)
        else:
            plt.xlim (5,200)
            plt.ylim (5,200)
        
       
        ax.tick_params(labelsize=12, width = 1.2)
        lims = [ np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
        plt.plot(lims, lims, linestyle='--', linewidth=0.8, color='black' )
        plt.xscale('log')
        plt.yscale('log')
        
        
        plt.ylabel('Modeled precipitation (mm)' ,fontsize=12)
        plt.xlabel('Observed precipitation (mm)' ,fontsize=12)
            
        #plt.title (tab[ii], loc ='center' ,weight ='bold')
        plt.title (label[ii]+' '+tab[ii], loc='left',fontsize=12)
        
        
    ax = plt.subplot (2,2,4)
    for kk in np.arange (0,len(pctl)):
        plt.scatter (obs_daypctl[:,3],mod11_daypctl[:,3] , color='green',marker ='o',s=28)
        plt.scatter (obs_daypctl[:,3],mod04_daypctl[:,3] , color='blue',marker ='D',s=20)
        r2 = np.corrcoef(obs_daypctl[:,3],mod11_daypctl[:,3])[0,1]**2
        bias = np.nanmean(mod11_daypctl[:,3]-obs_daypctl[:,3])
        
        r2_04 = np.corrcoef(obs_daypctl[:,3],mod04_daypctl[:,3])[0,1]**2
        bias_04 = np.nanmean(mod04_daypctl[:,3]-obs_daypctl[:,3])
        
        ax.text (0.06,0.80,'BIAS: '+str('%.2f' % round(bias,2))+' mm',color= 'green',fontsize=12,transform=ax.transAxes)
        ax.text (0.06,0.74,'R$^2$: '+str('%.2f' % round(r2,2)),color= 'green',fontsize=12,transform=ax.transAxes)

        #ax1.text (0.65,0.18,'CTRL04',color= 'blue',fontsize=12,transform=ax1.transAxes)
        ax.text (0.65,0.12,'BIAS: '+str('%.2f' % round(bias_04,2))+' mm',color= 'blue',fontsize=12,transform=ax.transAxes)
        ax.text (0.65,0.06,'R$^2$: '+str('%.2f' % round(r2_04,2)),color= 'blue',fontsize=12,transform=ax.transAxes)
     
    
    #plt.legend(['90.0','95.0','99.0','99.9'],fontsize=7,frameon=False)    
    plt.xlim (5,400)
    plt.ylim (5,400)
    ax.tick_params(labelsize=12, width = 1.2)

    lims = [ np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    plt.plot(lims, lims, linestyle='--', linewidth=0.8, color='black' )
    plt.xscale('log')
    plt.yscale('log')
    plt.title ('(d) Daily', loc ='left' ,fontsize=12)
    plt.ylabel('Modeled precipitation (mm)' ,fontsize=12)
    plt.xlabel('Observed precipitation (mm)' ,fontsize=12)
    
    

    fig.subplots_adjust(wspace=0.22, hspace =0.22, top =0.95, bottom=0.07, right =0.97, left=0.08)
    plt.savefig(filename, dpi=500)
 
scatter_percentile ('Scatter_pr_daily_hourly_CMAOBS_percentiles_validation.png')    
        


print ('the program done')
