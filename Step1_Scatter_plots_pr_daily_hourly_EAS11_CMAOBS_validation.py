#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:33:28 2019

Select hourly precipitation for 106 stations in REF_0.04, REF_0.11, and CMORPH
@author: shupli
"""


import xarray as xr

import numpy as np
#from search_station_in_model import lookupNearest, lookup_radius

import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)



dir2='/home/EAS-aerosol/EAS-11_aerosol/'

nd = 87648
nsta = 106


###### in this section, select gridpoint in CMORPH and model are close to 106 observed stations
#data = ['pr_EAS-11_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_day_2001-2010.nc',
#        'pr_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_day_2001-2010.nc']


eas11 = np.loadtxt ('/home/scripts/EAS-aerosol/1hr/MODEL/EAS11_indices_for_106station.txt')

print (eas11[:,0])



modsta11 =np.zeros (shape =(nd,nsta))
tt =0
for i in np.arange (0,10):
    mod11 =  xr.open_dataset (dir2+'REF_1hr/pr_EAS-11_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_1hr_'+str(2001+i)+'01010030-'+str(2001+i)+'12312330.nc')['pr']
    mod11 = mod11 *3600.
    
    print (i)
    for ista  in np.arange (0,nsta):    ## 106 observed stations
       
        modsta11[tt:tt+len(mod11[:,0,0]),ista] = mod11[:,int(eas11[ista,2]),int(eas11[ista,3])].data
    tt = tt + len(mod11[:,0,0])
     

np.savetxt ('/home/scripts/EAS-aerosol/1hr/MODEL/COSMO_EAS11_1hr_precipitation_106station_2001-2010.txt',
            modsta11, delimiter=',', fmt=' '.join(['%7.3f']*106))
    




print ('the program done')
