#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:33:28 2019

@author: shupli
"""


import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from mpl_toolkits.basemap import Basemap
from scipy import stats
import matplotlib.gridspec as gridspec
#import xesmf as xe
from mask_sea_land import  mask_sea, normalize180
#from bootstrap_routines_3D import *
from bootstrap_routines import *
import mask_sea_land
from matplotlib.ticker import (MultipleLocator)
from matplotlib.ticker import FixedFormatter, FixedLocator
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)

season=['DJF','MAM','JJA','SON']

var ='prhmax'

dir2='/home/EAS-aerosol/EAS-aerosol/'

dir1='/home/EAS-aerosol/EAS-hour/'

land_sea = xr.open_dataset('/home/EAS-aerosol/COORDINATE/land_sea_eas004.nc')

data = ['_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_mon_2001-2010.nc']

data1 =['_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux03_mon_2001-2010.nc',
        '_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux3_mon_2001-2010.nc',
        '_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx03_mon_2001-2010.nc',
        '_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx3_mon_2001-2010.nc']


datafre = ['pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_1hr_2001-2010_JJA_mm_']

datafre1 =['pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux03_1hr_2001-2010_JJA_mm_',
        'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux3_1hr_2001-2010_JJA_mm_',
        'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx03_1hr_2001-2010_JJA_mm_',
        'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx3_1hr_2001-2010_JJA_mm_']

rotpole = land_sea['rotated_pole']
lon = land_sea['lon'].data
lat = land_sea['lat'].data

lon_0 = normalize180(rotpole.grid_north_pole_longitude-180.)
rotlon = rotpole.grid_north_pole_longitude
rotlat = rotpole.grid_north_pole_latitude

aero = xr.open_dataset(dir2+'extpar_4.4km_EAS004_720x660_soil_landuse.nc', decode_times=False)
bc03 = aero['AER_BC12']*(-2./3.)
bc3 = aero['AER_BC12']*2.

su03 = aero['AER_SO412']*(-2./3.)
su3  = aero['AER_SO412']*2.

#sumean =su.sel(time=su['time.season']=='JJA').mean('time')
bc03 = bc03.where (land_sea['HSURF'].values==1)
bc03_mean =bc03[5:8,:,:].mean ('time')
bc03_mean = np.ravel (bc03_mean)
bc03_mean = bc03_mean[np.logical_not(np.isnan(bc03_mean))]



bc3 = bc3.where (land_sea['HSURF'].values==1)
bc3_mean =bc3[5:8,:,:].mean ('time')
bc3_mean = np.ravel (bc3_mean)
bc3_mean = bc3_mean[np.logical_not(np.isnan(bc3_mean))]


su03 = su03.where (land_sea['HSURF'].values==1)
su03_mean =su03[5:8,:,:].mean ('time')
su03_mean = np.ravel (su03_mean)
su03_mean = su03_mean[np.logical_not(np.isnan(su03_mean))]


su3 = su3.where (land_sea['HSURF'].values==1)
su3_mean =su3[5:8,:,:].mean ('time')
su3_mean = np.ravel (su3_mean)
su3_mean = su3_mean[np.logical_not(np.isnan(su3_mean))]


print (np.min(bc03_mean), np.max(bc03_mean))
print (np.min(bc3_mean), np.max(bc3_mean))

print (np.min(su03_mean), np.max(su03_mean))
print (np.min(su3_mean), np.max(su3_mean))


levels_diff =[-280,-30,-25,-20,-15,-10,-5,-2,2,5,10,15,20,25,30,290]
norm_diff = matplotlib.colors.BoundaryNorm(levels_diff,len(levels_diff))

run = ['SU/3' ,'SUx3','BC/3','BCx3']

def bootstrap_test (kk, mod, obs):
     #### student's t-test

    obs1= obs.sel(time=obs['time.season']==season[kk])
    mod1= mod.sel(time=mod['time.season']==season[kk])
    if kk ==0:
        obssum1 = obs1.groupby(obs1.time.dt.year,restore_coord_dims=True).mean("time")
        modsum1 = mod1.groupby(mod1.time.dt.year,restore_coord_dims=True).mean("time")
        obssum1.values[0,:,:] = (obs1.values[0,:,:]+obs1.values[1,:,:]+obs1.values[29,:,:])/3.
        modsum1.values[0,:,:] = (mod1.values[0,:,:]+mod1.values[1,:,:]+mod1.values[29,:,:])/3.
        for i in np.arange(0,9):
            obssum1.values[i+1,:,:] = (obs1.values[i*3+2,:,:]+obs1.values[i*3+3,:,:]+obs1.values[i*3+4,:,:])/3.                   
            modsum1.values[i+1,:,:] = (mod1.values[i*3+2,:,:]+mod1.values[i*3+3,:,:]+mod1.values[i*3+4,:,:])/3.     
                 
    else:
        obssum1=obs1.groupby(obs1.time.dt.year,restore_coord_dims=True).mean("time")
        modsum1=mod1.groupby(mod1.time.dt.year,restore_coord_dims=True).mean("time")
    
    lower_ci, upper_ci =bootci_diff(modsum1.data, obssum1.data, nboot=30,replacement=True,alpha=0.1)
        
   
    return lower_ci, upper_ci


bins_su03 =[-0.4,-0.2,-0.15,-0.05,0]
#bins_su03 =[-0.3,-0.2,-0.1,-0.05,0]
#bins_su3 =[0, 0.25, 0.4, 0.6, 0.9]
bins_su3 =[0, 0.2, 0.4, 0.6, 0.9]

bins_bc03 =[-0.03,-0.02,-0.01,-0.005,0]
bins_bc3 = [0,0.02,0.04,0.06,0.08]
#bins_su03 =[-0.3,-0.25,-0.2,-0.15,-0.1,0]
#bins_su3 =[0,0.2,0.4,0.6,0.8,1.0]

#bins_bc03 =[-0.03,-0.025,-0.02,-0.015,-0.01,0]
#bins_bc3 =[0,0.02,0.03, 0.04,0.06,0.08]
label=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)']

def contour_pr(filename):
    fig=plt.figure(figsize=(12,9),constrained_layout=False)
    #gs1 = gridspec.GridSpec(2, 4)
    #gs2 = gridspec.GridSpec(1, 4)
    #gs0 = fig.add_gridspec(nrows =3, ncols=3, wspace=0.25, hspace =0.15)
    #gs00 = gs0[-1,:].subgridspec(1, 4)
    
    cont04 = xr.open_dataset(dir2+var+data[0])[var]
    cont04 = cont04 * 86400.
    contsum04_sum = cont04.sel(time=cont04['time.season']=='JJA').mean('time')
    contsum04_sum = contsum04_sum.where (land_sea['HSURF'].values[0,:,:] ==1)
    
    prbinmean = np.zeros (shape=(len(bins_su3)-1,4))
    prbin25  = np.zeros (shape=(len(bins_su3)-1,4))
    prbin75 = np.zeros (shape=(len(bins_su3)-1,4))
    for xx, model in enumerate(data1):
        
        mod = xr.open_dataset(dir2 + var + model)[var]
        mod = mod *86400.
        
        modsum = mod.sel(time=mod['time.season']=='JJA').mean('time')
        modsum = modsum.where (land_sea['HSURF'].values[0,:,:] ==1)
  
        lower_ci, upper_ci = bootstrap_test(2, mod, cont04)
        diff = (modsum - contsum04_sum)/contsum04_sum *100
        
        
        ax=plt.subplot( 3, 4, xx+1)
        #ax=fig.add_subplot(gs00[0, xx])
        
        m= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=rotlon,o_lat_p=rotlat,llcrnrlat = lat[0,0], 
                    urcrnrlat = lat[-1,-1],llcrnrlon = lon[0,0], urcrnrlon = lon[-1,-1],resolution='l',area_thresh=10000.)
        m.drawcoastlines(linewidth=0.3)
        x,y = m(lon,lat)
        im1 =m.contourf(x,y,diff,levels_diff,norm=norm_diff,cmap=mask_sea_land.prdiffcolor1)
        lon1 = lon[::12,::12]
        lat1 = lat[::12,::12]
        m.scatter( lon1[(0.<=lower_ci[::12,::12]) | (upper_ci[::12,::12] <=0.)], lat1[(0.<=lower_ci[::12,::12]) | (upper_ci[::12,::12] <=0.)],s=0.3, marker='o',c='black',latlon=True)
       
        m.drawmeridians([100,110,120], labels=[0,0,0,0], color='grey', linewidth=1.2,dashes =[1,1],fontsize=12)
       
        if xx == 0:    
            m.drawparallels([20,30,40],labels=[1,0,0,0],color='grey', linewidth=1.2,dashes =[1,1],fontsize=12)
            plt.ylabel('pHmax (%)',fontsize=12,labelpad =40, weight ='bold') 
        else:
            m.drawparallels([20,30,40],labels=[0,0,0,0],color='grey', linewidth=1.2,dashes =[1,1],fontsize=12)
        
        if xx ==0 or xx==1:
            ax.text(0.28, 1.13, run[xx]+' - REF', transform=ax.transAxes, fontsize=12, weight='bold')
        if xx==2 or xx==3:
            ax.text(0.28, 1.13, run[xx]+' - REF', transform=ax.transAxes, fontsize=12, weight='bold')
            
        mask_sea (ax,m)
        aveland = np.ma.array((diff.mean()))
        plt.title ('Mean: '+str(aveland.round(2))+' %',fontsize=12,loc='right',pad=2)
        plt.title (label[xx],loc='left',fontsize=12,pad=2)


        diff = np.ravel (diff)
        diff = diff[np.logical_not(np.isnan(diff))]
        
        
        if xx==0:
            
            pw_temp = np.ravel (su03_mean)
            bins = bins_su03
        elif xx==1:
            pw_temp = np.ravel (su3_mean)
            bins = bins_su3
        elif xx==2:
           
            pw_temp = np.ravel (bc03_mean)
            bins = bins_bc03
        else:
            
            pw_temp = np.ravel (bc3_mean)
            bins = bins_bc3
        
        for ik in np.arange (0, len(bins)-1):
            prbin = np.ma.masked_where ( (pw_temp<bins[ik]) | (pw_temp>=bins[ik+1]), diff)
            prbin = prbin[~prbin.mask]
            
            prbin = prbin[np.logical_not(np.isnan(prbin))]
            
            print (len(prbin))
            
            prbinmean[ik,xx] = np.mean(prbin)
            prbin25[ik,xx] = np.percentile(prbin,25)
            prbin75[ik,xx] = np.percentile(prbin,75)
            
        
    cont04 = xr.open_dataset(dir1+'frequency/'+datafre[0]+'fre_monsum.nc')['pr']
    contsum04 = cont04.sel(time=cont04['time.season']=='JJA').sum('time')
    contsum04 = contsum04.where (land_sea['HSURF'].values[0,:,:] ==1)
    for xx, model in enumerate(datafre1):

            mod = xr.open_dataset(dir1 + 'frequency/' + model+'fre_monsum.nc')['pr']

            modsum = mod.sel(time=mod['time.season']=='JJA').sum('time')
            modsum = modsum.where (land_sea['HSURF'].values[0,:,:] ==1)


            lower_ci, upper_ci = bootstrap_test(2, mod, cont04)
            diff = (modsum - contsum04)/contsum04 *100

            ax=plt.subplot( 3, 4, xx+5)
            m= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=rotlon,o_lat_p=rotlat,llcrnrlat = lat[0,0],
                       urcrnrlat = lat[-1,-1],llcrnrlon = lon[0,0], urcrnrlon = lon[-1,-1],resolution='l',area_thresh=10000.)
            m.drawcoastlines(linewidth=0.3)
            x,y = m(lon,lat)
            im1 =m.contourf(x,y,diff,levels_diff,norm=norm_diff,cmap=mask_sea_land.prdiffcolor1)
            lon1 = lon[::12,::12]
            lat1 = lat[::12,::12]
            m.scatter( lon1[(0.<=lower_ci[::12,::12]) | (upper_ci[::12,::12] <=0.)], lat1[(0.<=lower_ci[::12,::12]) | (upper_ci[::12,::12] <=0.)],s=0.3, marker='o',c='black',latlon=True)

            m.drawmeridians([100,110,120], labels=[0,0,0,1], color='grey', linewidth=1.2,dashes =[1,1],fontsize=12)

            if xx == 0:
                m.drawparallels([20,30,40],labels=[1,0,0,0],color='grey', linewidth=1.2,dashes =[1,1],fontsize=12)
                plt.ylabel('Frequency (%)',fontsize=12,labelpad =40, weight ='bold')
            else:
                m.drawparallels([20,30,40],labels=[0,0,0,0],color='grey', linewidth=1.2,dashes =[1,1],fontsize=12)

            mask_sea (ax,m)
            aveland = np.ma.array((diff.mean()))
            plt.title ('Mean: '+str(aveland.round(2))+' %',fontsize=12,loc='right',pad=2)
            plt.title (label[xx+4],loc='left',fontsize=12,pad=2)
            
            #ax=fig.add_subplot(gs00[0, xx])
            #ax=plt.subplot( 3, 4, xx+9)
            ax= fig.add_axes([0.08+xx*0.218, 0.1, 0.198, 0.22])
            diff = np.ravel (diff)
            diff = diff[np.logical_not(np.isnan(diff))]
        
            if xx==0:
                pw_temp = np.ravel (su03_mean)
                bins = bins_su03
            elif xx==1:
                pw_temp = np.ravel (su3_mean)
                bins = bins_su3
            elif xx==2:
                pw_temp = np.ravel (bc03_mean)
                bins = bins_bc03
            else:
                pw_temp = np.ravel (bc3_mean)
                bins = bins_bc3
        
            for ik in np.arange (0, len(bins)-1):
                prbin = np.ma.masked_where ( (pw_temp<bins[ik]) | (pw_temp>=bins[ik+1]), diff)
                prbin = prbin[~prbin.mask]
            
                prbin = prbin[np.logical_not(np.isnan(prbin))]
                print (len(prbin))
                            
                l1, =plt.plot (ik*1.2, prbinmean[ik,xx],marker='o',markersize=10, color = 'green',linestyle = 'None')
                plt.plot ([ik*1.2,ik*1.2], [prbin25[ik,xx],prbin75[ik,xx]], color = 'gray',linestyle = 'solid')
                l2, =plt.plot (ik*1.2+0.3, np.mean(prbin),marker='o',markersize=10, color = 'blue',linestyle = 'None')
                plt.plot ([ik*1.2+0.3,ik*1.2+0.3], [np.percentile(prbin,25),np.percentile(prbin,75)], color = 'gray',linestyle = 'solid')
               
        
            plt.ylim(-20,20)
            plt.xlim(-0.2,4.2)
            plt.axhline (y=0,linestyle='dashed',color='m',linewidth=0.8) 
            plt.axvline (x=4.2,linestyle='dashed',color='gray',linewidth=0.5) 
           # plt.title (label[xx+8],loc='left',fontsize=12,pad=2)
            ax.tick_params(width = 1.2,labelsize=11)
            ax.text (0.05,16, label[xx+8],fontsize=12)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.yaxis.set_major_locator(MultipleLocator(10))
            ax.yaxis.set_minor_locator(MultipleLocator(5))
            ax = plt.gca()
            ax.yaxis.grid(True,linestyle='dashed', linewidth=0.5)
            ax.xaxis.set_major_locator(FixedLocator([0.15,1.35,2.55,3.75]))
            if xx !=0:
                plt.yticks([-20,-10,0,10,20],[' ',' ',' ',' ',' '], fontsize=12)
            else:
                plt.ylabel('Relative change (%)',fontsize=12,labelpad =9, weight ='bold')
                lgnd=ax.legend ([l1,l2],['pHmax','Frequency'], loc ='upper right',frameon=False,facecolor="white",fontsize=9)
                
                lgnd.legendHandles[0]._legmarker.set_markersize(6)
                lgnd.legendHandles[1]._legmarker.set_markersize(6)
                
            if xx==0:
                ax.xaxis.set_major_formatter(FixedFormatter(['[-0.3,-0.2]','[-0.2,-0.1]','[-0.1,-0.05]','[-0.05,0]']))
                plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
                plt.xlabel('Reduced AOD',fontsize=12)
            elif xx==1:
                ax.xaxis.set_major_formatter(FixedFormatter(['[0,0.2]','[0.2,0.4]','[0.4,0.6]','[0.6,0.9]']))
                plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
                plt.xlabel('Increased AOD',fontsize=12)
            
            elif xx==2:
                ax.xaxis.set_major_formatter(FixedFormatter(['[-0.3,-0.2]','[-0.2,-0.1]','[-0.1,-0.05]','[-0.05,0]']))
                plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
                plt.xlabel('Reduced AOD ($10^{-1}$)',fontsize=12)
            
            else:
                ax.xaxis.set_major_formatter(FixedFormatter(['[0,0.02]','[0.2,0.4]','[0.4,0.6]','[0.6,0.8]']))
                plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
                plt.xlabel('Increased AOD ($10^{-1}$)',fontsize=12)
                
                
            
            
        
   
        
    #cb_ax.text(1.02, 0.2,'SSR (W/m$^2$)',transform = cb_ax.transAxes,fontsize=12)
    cb_ax=fig.add_axes([0.94, 0.43, 0.013, 0.42])
    cb=plt.colorbar(im1, cax=cb_ax, pad=0.05,shrink = 0.6, orientation = 'vertical',ticks =[-30,-25,-20,-15,-10,-5,-2,2,5,10,15,20,25,30])
    cb.ax.tick_params(labelsize=11)
    #cb.set_label('Relative difference in daily maximum hourly precipitation (%)',fontsize=10)
    
    
    
    fig.subplots_adjust(hspace=0.06,wspace=0.1,top =0.92, bottom=0.06,left=0.08, right=0.93)
    plt.savefig(filename, dpi=500)
              
contour_pr('Contour_plots_prhmax_wethour_difference_aerosol.png')

print ('the program done')
