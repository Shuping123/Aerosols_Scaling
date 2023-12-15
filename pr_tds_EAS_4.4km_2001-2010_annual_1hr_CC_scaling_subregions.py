#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:33:28 2019
Clausius-Claperyon equation scale
@author: shupli
"""

#from netCDF4 import Dataset
#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.cm as cm
import xarray as xr
#from matplotlib.colors import ListedColormap
from scipy import stats
#from collections import OrderedDict

import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)



dir2='/home/EAS-aerosol/EAS-aerosol/'
dir1='/home/scripts/EAS-aerosol/subdaily/CCscale/data/'



data =['pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_1hr_2001-2010_5x5max.nc',
       'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux03_1hr_2001-2010_5x5max.nc',
       'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux3_1hr_2001-2010_5x5max.nc',
       'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx03_1hr_2001-2010_5x5max.nc',
       'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx3_1hr_2001-2010_5x5max.nc']

temp =['tds_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_day_2001-2010_5x5mean.nc',
       'tds_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux03_day_2001-2010_5x5mean.nc',
       'tds_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux3_day_2001-2010_5x5mean.nc',
       'tds_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx03_day_2001-2010_5x5mean.nc',
       'tds_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx3_day_2001-2010_5x5mean.nc']

temp1 =['tas_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_day_2001-2010_5x5mean.nc',
       'tas_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux03_day_2001-2010_5x5mean.nc',
       'tas_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_sux3_day_2001-2010_5x5mean.nc',
       'tas_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx03_day_2001-2010_5x5mean.nc',
       'tas_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_bcx3_day_2001-2010_5x5mean.nc']



subreg04=xr.open_dataset('/home/EAS-aerosol/COORDINATE/Subregions_EAS004_sample_5x5.nc')


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


run =['REF','SUx03','SUx3','BCx03','BCx3']

subreg=['NC','SE']
def CC_scaling_subreg (reg):
 for xx, model in enumerate (data):
    pr = xr.open_dataset (dir2+data[xx])['pr']
    ta= xr.open_dataset (dir2+temp[xx])['tds']
    ta = ta -273.15
    
    #### select temperature above zero
    tas = xr.open_dataset (dir2+temp1[xx])['tas']
    tas = tas-273.15
    tas = np.repeat (tas, 24, axis=0)
    
    modt1 = xr.open_dataset (dir2+data[0])['pr']
    modt1.values = tas
    
    
    ta = np.repeat (ta, 24, axis=0)
    modt = xr.open_dataset (dir2+data[0])['pr']
    modt.values = ta
    
    
    pr = pr.where ((subreg04[reg].values==1) & (modt1.values>0))
    modt = modt.where((subreg04[reg].values==1) & (modt1.values>0))
        
    pr= np.ravel(pr)
    pr = pr[np.logical_not(np.isnan(pr))]
    
    modt = np.ravel(modt)
    modt = modt[np.logical_not(np.isnan(modt))]
    
    bins= histedges_equalN(modt, 20)
   
    
    ### 99th percentile
    binmid_cont, bin_means_cont= cc_scale (pr, modt, 99.0, bins)
    ### 95th percentile
    binmid_cont9, bin_means_cont9 = cc_scale (pr, modt, 99.9, bins)
    ### 90th percentile
    binmid_cont99, bin_means_cont99 = cc_scale (pr, modt, 99.99, bins)
    
    binmid_cont999, bin_means_cont999 = cc_scale (pr, modt, 99.999, bins)
    
    np.savetxt(dir1+'COSMO_EAS04_1hr_pr_tds_CC_scaling_bin_'+run[xx]+'_'+reg+'_max.txt', 
                   np.hstack((binmid_cont,bin_means_cont,bin_means_cont9,bin_means_cont99,bin_means_cont999)), delimiter=',', fmt=' '.join(['%8.4f']*5) )
CC_scaling_subreg ('NC')
CC_scaling_subreg ('SE')
    
print ('the program done')


