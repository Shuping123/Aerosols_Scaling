#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:33:28 2019

@author: shupli
"""

#import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
from matplotlib.colors import ListedColormap
#from scipy import stats
from matplotlib.ticker import FixedFormatter, FixedLocator
#import xesmf as xe
from matplotlib.ticker import (MultipleLocator)

dir2='/home/EAS-aerosol/EAS-hour/diurnal_cycle/'

dir1='/home/scripts/EAS-aerosol/1hr/MODEL/'

nsta= 106
nt =24 
nd = 87648

data1 = np.loadtxt ('/home/EAS-aerosol/Validation/Stations_info_selected_from_166stations.txt')
#eas11 = np.loadtxt ('/home/EAS-aerosol/Validation/EAS11_indices_for_106station.txt')
#eas04 = np.loadtxt ('/home/EAS-aerosol/Validation/EAS004_indices_for_106station.txt')

modsta04 = np.loadtxt (dir1 +'COSMO_EAS004_1hr_precipitation_106station_2001-2010.txt')
modsta11 = np.loadtxt (dir1 +'COSMO_EAS11_1hr_precipitation_106station_2001-2010.txt')

pmax11 = np.loadtxt(dir1+'COSMO_EAS11_prhmax_precipitation_106station_2001-2010.txt')
pmax04 = np.loadtxt(dir1+'COSMO_EAS004_prhmax_precipitation_106station_2001-2010.txt')

        
obssta = np.zeros (shape=(nd,nsta))       

for ista in np.arange (0, nsta): 
    ### read observations for 106 stations
    data = np.loadtxt ('/home/EAS-aerosol/Validation/OBS_1hr/SURF_CLI_CHN_PRE_FTM_H24-QC-'+str(int(data1[ista,0]))+'.TXT')
    
    obssta [:,ista] = data[:,7:31].reshape(nd)  ### daily precipitation no NAN values
    if np.min(obssta[:,ista]) == np.max(obssta[:,ista]):
        print (np.max(obssta[:,ista]),data[0,ista])


### shift Beijing Time to UTC
### UTC Time: 2001.06.01.(00-01)UTC - 2010.8.31.(23-24)UTC 
### BJ Time: 2001.06.01.(08-09)BJS - 2010.9.1.(07-08)BJS
nsday = 22080 ; nwday = 21624
obssta_UTC = np.zeros (shape=(nsday,nsta))
obsstawin_UTC = np.zeros (shape=(nwday,nsta))

mod11 = np.zeros (shape=(nsday,nsta))
mod04 = np.zeros (shape=(nsday,nsta))

nn = 0
nn1 =0 ; kk1 =0
for i in np.arange (0,10):
    print (i)
    #if i ==0:
    #    obssta_UTC[i*92*24:(i+1)*92*24,:] = obssta[151*24+12:(151+92)*24+12,:]
    #    nn =365*24
    if i==3 or i ==7:
        obssta_UTC[i*92*24:(i+1)*92*24,:] = obssta[nn+152*24+12:nn+(152+92)*24+12,:]
        mod11[i*92*24:(i+1)*92*24,:] = modsta11[nn+152*24:nn+(152+92)*24,:]
        mod04[i*92*24:(i+1)*92*24,:] = modsta04[nn+152*24:nn+(152+92)*24,:]
        nn = nn+366*24
        
    else:
        obssta_UTC[i*92*24:(i+1)*92*24,:] = obssta[nn+151*24+12:nn+(151+92)*24+12,:]
        mod11[i*92*24:(i+1)*92*24,:] = modsta11[nn+151*24:nn+(151+92)*24,:]
        mod04[i*92*24:(i+1)*92*24,:] = modsta04[nn+151*24:nn+(151+92)*24,:]
        nn = nn+365*24
        
    if i ==3 or i ==7:
        obsstawin_UTC[kk1:kk1+60*24,:] = obssta[nn1+12:nn1+60*24+12,:]
        obsstawin_UTC[kk1+60*24:kk1+91*24,:] = obssta[nn1+335*24+12:nn1+366*24+12,:]
        kk1 = kk1+ 91*24
        nn1 = nn1 +366*24
        print (kk1+60*24,kk1+91*24)
    else:
        obsstawin_UTC[kk1:kk1+59*24,:] = obssta[nn1+12:nn1+59*24+12,:]
        if i !=9 :
            obsstawin_UTC[kk1+59*24:kk1+90*24,:] = obssta[nn1+334*24+12:nn1+365*24+12,:]
            print (kk1+59*24,kk1+90*24)
        else:
            obsstawin_UTC[kk1+59*24:kk1+89*24,:] = obssta[nn1+334*24+12:nn1+364*24+12,:]
            print (kk1+59*24,kk1+89*24)
        
        kk1 = kk1 + 90*24
        nn1 = nn1 + 365*24
#### calculate daily maximum hour precipitation in summer
obs1 = np.where  (obssta_UTC== 32766, -999, obssta_UTC)
        
### mask missing values in the model
mod11= np.ma.masked_where(obssta_UTC== 32766,mod11) 
mod04= np.ma.masked_where(obssta_UTC== 32766,mod04)       
    
obs = np.ma.masked_where(obssta_UTC== 32766,obssta_UTC)
obswin = np.ma.masked_where(obsstawin_UTC== 32766,obsstawin_UTC)


obs_dcycle = np.zeros (shape =(nt,nsta)) ### mean precipitation, diurnal cycle
mod11_dcycle = np.zeros (shape =(nt,nsta))
mod04_dcycle = np.zeros (shape =(nt,nsta))

obs_whour = np.zeros (shape =(nt,nsta)) ### wet hour frequency,  diurnal cycle
mod11_dhour = np.zeros (shape =(nt,nsta))
mod04_dhour = np.zeros (shape =(nt,nsta))


r2_eas11 = np.zeros (shape =(nsta,2)) 
r2_eas04 = np.zeros (shape =(nsta,2))

obs_prmax = np.zeros (shape=(920,nsta))
prmax04_sum = np.zeros (shape=(920,nsta))
prmax11_sum = np.zeros (shape=(920,nsta))

nx =0 ; nk =0
for i in np.arange (0,10):
    if i ==3 or i ==7:
        prmax11_sum[i*92:(i+1)*92,:] = pmax11[nx+152:nx+152+92,:]
        prmax04_sum[i*92:(i+1)*92,:] = pmax04[nx+152:nx+152+92,:]
        nx = nx + 366
    else:
        prmax11_sum[i*92:(i+1)*92,:] = pmax11[nx+151:nx+151+92,:]
        prmax04_sum[i*92:(i+1)*92,:] = pmax04[nx+151:nx+151+92,:]
        nx = nx + 365

      
for i in np.arange(0,920):
    obs_prmax[i,:] = np.max (obs1[i*24:(i+1)*24,:],axis=0)
    print (obs_prmax)

#obs_prmax = obs_prmax*0.1
#print (np.max(obs_prmax))
# Calculate cumulative frequency
def get_cumfreq(val):
    weights = np.ones_like(val)/float(len(val))
    max_int = int(np.max(val)) + 1
    num_bins = 10 * max_int
    freq, base = np.histogram(val, bins=num_bins, range=(0, max_int),
                              weights=weights)
    freq = freq[::-1]
    cumfreq = np.cumsum(freq)
    cumfreq = cumfreq[::-1]
    base = base[1:-1]
    cumfreq = cumfreq[1:]
    return base, cumfreq

prmax11_sum = np.where (obs_prmax>=0, prmax11_sum,obs_prmax)
prmax04_sum = np.where (obs_prmax>=0, prmax04_sum,obs_prmax)

obs_prmax = np.ma.masked_where(obs_prmax == -999, obs_prmax)
prmax11_sum = np.ma.masked_where(prmax11_sum == -999, prmax11_sum)
prmax04_sum = np.ma.masked_where(prmax04_sum == -999, prmax04_sum)

obs_prmax = np.ravel (obs_prmax)
prmax11_sum = np.ravel (prmax11_sum)
prmax04_sum = np.ravel (prmax04_sum)

obs_prmax = obs_prmax[~obs_prmax.mask]
prmax11_sum = prmax11_sum[~prmax11_sum.mask]
prmax04_sum = prmax04_sum[~prmax04_sum.mask]


for i in np.arange (0,nt):
    obs_dcycle[i,:] = np.nanmean (obs[i:nsday:24,:],axis=0)
    mod11_dcycle[i,:] = np.nanmean (mod11[i:nsday:24,:],axis=0)
    mod04_dcycle[i,:] = np.nanmean (mod04[i:nsday:24,:],axis=0)
    
    for ist in np.arange (0,nsta):
        temp = obs[i:nsday:24,ist]*0.1
        print (np.shape(temp))
        temp = temp[~temp.mask]
        obs_whour[i,ist] = np.count_nonzero (temp>0.1)/(len(temp)*1.0)
        mod11_dhour[i,ist]= np.count_nonzero (mod11[i:nsday:24,ist]>0.1)/(len(temp)*1.0)
        mod04_dhour[i,ist]= np.count_nonzero (mod04[i:nsday:24,ist]>0.1)/(len(temp)*1.0)
        
        #obs_whour[i,ist]= np.count_nonzero ((obs[i:nsday:24,ist]*0.1)>0.1)
        #mod11_whour[i,ist]= np.count_nonzero (mod11[i:nsday:24,ist]>0.1)
        #mod04_whour[i,ist]= np.count_nonzero (mod04[i:nsday:24,ist]>0.1)


obs_dcycle = obs_dcycle * 0.1 ###mm/h
for ista in np.arange (0,nsta):
    r2_eas11[ista,0] = np.corrcoef (mod11_dcycle[:,ista], obs_dcycle[:,ista])[0,1]
    r2_eas04[ista,0] = np.corrcoef (mod04_dcycle[:,ista], obs_dcycle[:,ista])[0,1]
    
    ### wet hour frequency
    r2_eas11[ista,1] = np.corrcoef (mod11_dhour[:,ista], obs_whour[:,ista])[0,1]
    r2_eas04[ista,1] = np.corrcoef (mod04_dhour[:,ista], obs_whour[:,ista])[0,1]
    




tt = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]


mod_dcycle = np.zeros (shape =(len(tt),nsta,2)) ### EAS11 and EAS04 in summer
mod_whour = np.zeros (shape =(len(tt),nsta,2))

obscycle = np.zeros (shape =(len(tt),nsta))
obswhour = np.zeros (shape =(len(tt),nsta))

mod_dcycle[0:24,:,0] = mod11_dcycle ; mod_dcycle[24,:,0] = mod11_dcycle[0,:]
mod_dcycle[0:24,:,1] = mod04_dcycle ; mod_dcycle[24,:,1] = mod04_dcycle[0,:]

mod_whour[0:24,:,0] = mod11_dhour ; mod_whour[24,:,0] = mod11_dhour[0,:]
mod_whour[0:24,:,1] = mod04_dhour ; mod_whour[24,:,1] = mod04_dhour[0,:]

obscycle[0:24,:] = obs_dcycle ; obscycle[24,:] = obs_dcycle[0,:]
obswhour[0:24,:] = obs_whour ; obswhour[24,:] = obs_whour[0,:]

"""
label=[' (d)',' (e)']
#fig = plt.figure(figsize=(12,8))
def lineplots(filename):
    
    fig = plt.figure(figsize=(12,8),constrained_layout=False)
    gs0 = fig.add_gridspec(nrows =2, ncols=3,left=0.08, right=0.94, wspace=0.28, hspace =0.28, top = 0.94,bottom=0.08)
    gs00 = gs0[0,:].subgridspec(1, 2)
   
    ax=fig.add_subplot(gs00[0, 0])
    l1, =plt.plot(tt,np.nanmean(mod_dcycle[:,:,0],axis=1),color='C0',linewidth=2.0)
    l2, =plt.plot(tt,np.nanmean(mod_dcycle[:,:,1],axis=1),color='C3',linewidth=2.0)
    l3, =plt.plot(tt,np.nanmean(obscycle,axis=1),color='gray',linewidth=3)
       
    
    plt.grid(linestyle='dashed', linewidth= 0.5)           
    plt.xlim(tt[0],tt[10])
    plt.ylim(0.1,0.24)
    #plt.ylim (ymax=0.24)
    plt.xlabel ('Time (UTC)',fontsize=12)
    plt.ylabel ('Precipitation (mm/h)',fontsize=12)
    plt.title (' (a) Mean', loc ='left')
  
    plt.xticks ([0,4,8,12,16,20,24])
    
    #ax.xaxis.set_minor_locator(MultipleLocator(1))
    #ax.yaxis.set_major_locator(MultipleLocator(0.005))
    plt.legend ([l3,l1,l2],['OBS','REF_0.11','REF_0.04'], loc ='lower right',frameon=False,facecolor="white",fontsize=10)
    ax.xaxis.set_major_locator(FixedLocator([0,4,8,12,16,20,24]))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FixedFormatter(['00','04','08','12','16','20','24']))
    
    
    
    ax=fig.add_subplot(gs00[0, 1])
    l1, =plt.plot(tt,np.nanmean(mod_whour[:,:,0],axis=1),color='C0',linewidth=2.0)
    l2, =plt.plot(tt,np.nanmean(mod_whour[:,:,1],axis=1),color='C3',linewidth=2.0)
    l3, =plt.plot(tt,np.nanmean(obswhour,axis=1),color='gray',linewidth=3)
       
    
    plt.grid(linestyle='dashed', linewidth= 0.5)           
    plt.xlim(tt[0],tt[10])
    plt.ylim(0.06,0.2)
    #plt.ylim (ymax=0.24)
    plt.xlabel ('Time (UTC)',fontsize=12)
    plt.ylabel ('Frequency',fontsize=12)
    plt.title (' (b) Wet hour frequency', loc ='left')
  
    plt.xticks ([0,4,8,12,16,20,24])
    plt.legend ([l3,l1,l2],['OBS','REF_0.11','REF_0.04'], loc ='upper right',frameon=False,facecolor="white",fontsize=10)
    ax.xaxis.set_major_locator(FixedLocator([0,4,8,12,16,20,24]))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FixedFormatter(['00','04','08','12','16','20','24']))
    
    #ax.xaxis.set_minor_locator(MultipleLocator(1))
    #ax.yaxis.set_major_locator(MultipleLocator(0.005))
   
    
    
    gs01 = gs0[1,:].subgridspec(1, 3, wspace=0.26)
    for kk in np.arange(0,2):
        ax=fig.add_subplot(gs01[0, kk+1])
        plt.scatter (r2_eas11[:,kk]**2,r2_eas04[:,kk]**2 , marker ='o',s=32, color='C0',edgecolors='C0',linewidth=0.8)
        #plt.scatter (r2_eas11,r2_eas04 , color='C3',marker ='o',s=3)
        plt.ylim(-0.1,1)
        plt.xlim(-0.1,1)
        lims = [ np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
        plt.plot(lims, lims, linestyle='--', linewidth=0.8, color='black' )
        if kk  ==  0:
            plt.title (label[kk]+' R$^2$ of mean', loc ='left')
        else:
            plt.title (label[kk]+' R$^2$ of frequency', loc ='left')
        plt.ylabel ('R$^2$ (REF_0.04 vs. OBS)',fontsize=11)
        plt.xlabel ('R$^2$ (REF_0.11 vs. OBS)',fontsize=11)
    
    
    
    ax=fig.add_subplot(gs01[0, 0])
    base1,cumfreq1 = get_cumfreq(obs_prmax*0.1)
    base2,cumfreq2 = get_cumfreq(prmax11_sum)
    base3,cumfreq3 = get_cumfreq(prmax04_sum)
    
    l3, = plt.plot(base1,cumfreq1,color='gray',linewidth=3.0)
    l1, = plt.plot(base2,cumfreq2,color='C0',linewidth=2.0)
    l2, = plt.plot(base3,cumfreq3,color='C3',linewidth=2.0)
    
    plt.ylim(ymin=10e-4,ymax = 1)
    ax.set_yscale('log')
    plt.xlim(xmin=0,xmax = 60)
    plt.title (' (c) pHmax', loc ='left')
    plt.ylabel ('Cumulative frequency',fontsize=12)
    plt.xlabel ('Precipitation (mm/h)',fontsize=12)
    plt.legend ([l3,l1,l2],['OBS','REF_0.11','REF_0.04'], loc ='upper right',frameon=False,facecolor="white",fontsize=10)
    
    plt.subplots_adjust(right=0.96, bottom=0.15, left =0.08, top=0.9)        
    plt.savefig(filename,dpi=500)  

        
lineplots('Lineplots_pr_whour_phmax_EAS_2001-2010_diurnal_cycle_validataion.png')
"""

print ('the program done')
