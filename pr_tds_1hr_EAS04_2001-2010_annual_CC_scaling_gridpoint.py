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
from Regression import  lag_linregress_3D
#from matplotlib.colors import ListedColormap
from scipy import stats
from collections import OrderedDict
import pickle

import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)

nlon = 144
nlat = 132
nh = 87648

dir1='/home/EAS-aerosol/EAS-hour/'
dir2='/home/EAS-aerosol/EAS-aerosol/'


data =['pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_1hr_2001-2010_5x5max.nc',
       'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux03_1hr_2001-2010_5x5max.nc',
       'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux3_1hr_2001-2010_5x5max.nc',
       'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx03_1hr_2001-2010_5x5max.nc',
       'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx3_1hr_2001-2010_5x5max.nc'
       ]

temp =[
       'tds_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_day_2001-2010_5x5mean.nc',
       'tds_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux03_day_2001-2010_5x5mean.nc'
       'tds_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux3_day_2001-2010_5x5mean.nc',
       'tds_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx03_day_2001-2010_5x5mean.nc',
       'tds_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx3_day_2001-2010_5x5mean.nc'
       ]

temp1 =['tas_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_day_2001-2010_5x5mean.nc',
       'tas_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux03_day_2001-2010_5x5mean.nc',
       'tas_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux3_day_2001-2010_5x5mean.nc',
       'tas_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx03_day_2001-2010_5x5mean.nc',
       'tas_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx3_day_2001-2010_5x5mean.nc'
       ]




def histedges_equalN(x, nbin):
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

def cc_scale (pre, ts, pth, bins ):
    ## pre: precipitation (one dimension)
    ## ts: air temperature (dew point temperature) (one dimension)
    ## pth: percentile
    ## one can change temperature range

    bin_means, bin_edges, binnumber = stats.binned_statistic(ts, pre, statistic=lambda y: np.percentile(y, pth), bins=bins)
    bin_width = (bin_edges[1] - bin_edges[0])
    binmid= (bin_edges[:-1]+ bin_edges[1:])/2.
    
    ## binmid: the meddile value in each bin_edges
    ## bin_means: statistics in each bin
    
    binmid.shape = len (binmid), 1
    bin_means.shape = len (bin_means),1
    
    return binmid, bin_means

#subreg=['NC','CC','SE']
#linestyles =['dashed',linestyles_dict['densely dotted']]

run =['SUx03']


def save(data,file_name):
    filehandler = open(file_name, 'wb')
    pickle.dump(data,filehandler)
    filehandler.close()
    

slope_t = np.zeros (shape =(3, nlat,nlon)) ###99.0,99.9,99.99 
pval_t = np.zeros (shape =(3, nlat,nlon))

for ks, exp in enumerate (run):
    print (exp)
    pr = xr.open_dataset (dir2+data[ks])['pr']
    
    tds= xr.open_dataset (dir2+temp[ks])['tds'].values
    tds = tds -273.15
    tds = np.repeat (tds, 24, axis=0)
    
        #### select temperature above zero
    tas = xr.open_dataset (dir2+temp1[ks])['tas'].values
    tas = tas-273.15
    tas = np.repeat (tas, 24, axis=0)
    
    
    ### to calculate C-C scaling at each grid point
    for ilon in np.arange (0,nlon):
        print (ilon)
        for ilat in np.arange (0,nlat):
            pre_reg = np.ma.masked_where ( (tas[:,ilat,ilon]<=0) , pr[:,ilat,ilon])
            tds_reg = np.ma.masked_where ( (tas[:,ilat,ilon]<=0) , tds[:,ilat,ilon])
            
            pre_reg=pre_reg[~pre_reg.mask]
            tds_reg=tds_reg[~tds_reg.mask]
            
            pre_reg = np.ravel (pre_reg)
            tds_reg = np.ravel (tds_reg)
            #print (np.shape (pre_reg))
            
            if (len (pre_reg) >=5000):
            
                bins= histedges_equalN(tds_reg, nbin=10)
            
                
                binmid0, bin_pr975 = cc_scale (pre_reg, tds_reg, 97.5, bins)
                binmid1, bin_pr99 = cc_scale (pre_reg, tds_reg, 99.0, bins)
                binmid2, bin_pr999 = cc_scale (pre_reg, tds_reg, 99.9, bins)
                
                y0 = np.log10(bin_pr975) 
                y1 = np.log10(bin_pr99)
                y2 = np.log10(bin_pr999)
                
            
                cov,cor,slope,intercept,pval,stderr = lag_linregress_3D(len(binmid0),binmid0,y0)
                slope_t[0, ilat, ilon] =  10**slope -1
                pval_t[0, ilat, ilon] =  pval
                
                cov,cor,slope,intercept,pval,stderr = lag_linregress_3D(len(binmid1),binmid1,y1)
                slope_t[1, ilat, ilon] =  10**slope -1
                pval_t[1, ilat, ilon] =  pval
                
                cov,cor,slope,intercept,pval,stderr = lag_linregress_3D(len(binmid2),binmid2,y2)
                slope_t[2, ilat, ilon] =  10**slope -1
                pval_t[2, ilat, ilon] =  pval
            

            else:
                slope_t[:, ilat, ilon] = [-99999, -99999, -99999]
                pval_t[:, ilat, ilon] =  [99999, 99999, 99999]
                
            
     
            
    
    save (slope_t,dir1+'pr_max_C-C_scaling_1hr_'+run[ks]+'_gridpoint_975_99_999.pkl')
    save (pval_t,dir1+'pr_max_C-C_scaling_1hr_'+run[ks]+'_pval_gridpoint_975_99_999.pkl')
    
 
print ('the program done')


