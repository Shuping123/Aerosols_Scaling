#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:33:28 2019

@author: shupli
"""

from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from mpl_toolkits.basemap import Basemap
import xarray as xr
import mask_sea_land
#import xesmf as xe
#from matplotlib.ticker import ( NullFormatter,ScalarFormatter, AutoMinorLocator)



dir1 = '/home/EAS-aerosol/EAS-aerosol/'

nlat = 660
nlon = 720
nmon = 12


##### AerCom-1 aerosol data (contains five aerosol components)
var=['BC12','ORG12','SO412','SS12','DUST12']

#extpar= xr.open_dataset(dir1+'extpar_4.4km_EAS004_720x660_soil_landuse.nc')['AER_BC12']
land_sea = xr.open_dataset ('/home/EAS-aerosol/COORDINATE/land_sea_eas004.nc')['HSURF']

aer = xr.open_dataset(dir1+'extpar_4.4km_EAS004_720x660_soil_landuse.nc', decode_times=False)

aod_com_ann = np.zeros(shape=(5,nmon,nlat,nlon))  # 5 types of aerosol
for i in np.arange(0,5):
    print (i,var[i])
    aod = aer['AER_'+var[i]]
    aod_com_ann[i,:,:,:] = aod.values

hurs = xr.open_dataset(dir1+'hurs_EAS-004_ECMWF-ERA5_evaluation_r1i1p1_CLMcom-ETH-COSMO-crCLIM-v1-1_v1_mon_2001-2010.nc')['hurs']
#hursmean = hurs.sel(time=hurs['time.season']=='JJA').mean('time')
hursmean  = hurs.mean(axis=0)

levels2 = [0,0.01,0.02,0.03,0.04,0.05,0.10,0.15,0.2,0.25,0.3,0.35,1.2]
norm2=matplotlib.colors.BoundaryNorm(levels2,len(levels2))

levels1 = [0,10,20,30,40,45,50,60,70,75,80,85,90,100]
norm1 =  matplotlib.colors.BoundaryNorm(levels1,len(levels1))

jet = cm.get_cmap('jet', 256)
jetBig = jet (np.linspace(0.4, 1.0, 256))
jetBig1 = jet (np.linspace(0.0, 0.7, 256))
white = np.array([255/256, 255/256, 255/256, 1])
new =np.zeros(shape=(14,4))
for i in np.arange(0,13):
    new[i+1,:]=jetBig[60+i*15,:]
new[0,:]=white
newcmp1=ListedColormap(new)


data = Dataset(dir1+'extpar_4.4km_EAS004_720x660_soil_landuse.nc')
lat = data.variables['lat'][:,:]
lon = data.variables['lon'][:,:]
rotpole = data.variables['rotated_pole']
def normalize180(lon):
    """Normalize lon to range [180, 180)"""
    lower = -180.; upper = 180.
    if lon > upper or lon == lower:
        lon = lower + abs(lon + upper) % (abs(lower) + abs(upper))
    if lon < lower or lon == upper:
        lon = upper - abs(lon - lower) % (abs(lower) + abs(upper))
    return lower if lon == upper else lon

lon_0 = normalize180(rotpole.grid_north_pole_longitude-180.)
o_lon_p = rotpole.grid_north_pole_longitude
o_lat_p = rotpole.grid_north_pole_latitude

run = ['BC','OC','SU','SS','DU']
tt= np.arange(0,12)+1
color=['black','green','red','m','blue']
#labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']  
label =['(a) ', '(b) ','(c) ','(d) ','(e) ']  

def subreg (ax, m):
    lon_sc1 = np.linspace(105,105)
    lat_sc1 = np.linspace(18,28)
    lon_sc2 = np.linspace(105,122)
    lat_sc2 = np.linspace(28,28)
    lon_sc3 = np.linspace(122,122)
    lat_sc3 = np.linspace(18,28)
    lon_sc4 = np.linspace(105,122)
    lat_sc4 = np.linspace(18,18)

    x1,y1 = m(lon_sc1,lat_sc1)
    m.plot (x1, y1,  linewidth=1.0 ,color='black')
    x1,y1 = m(lon_sc2,lat_sc2)
    m.plot (x1, y1,  linewidth=1.0, color='black')
    x1,y1 = m(lon_sc3,lat_sc3)
    m.plot (x1, y1,  linewidth=1.0, color='black')
    x1,y1 = m(lon_sc4,lat_sc4)
    m.plot (x1, y1,  linewidth=1.0, color='black')
    ax.annotate('SE', xy=(0.37, 0.36), xycoords='axes fraction',fontsize=10,color='blue',weight='bold')

    ## north China
    lon_sc1 = np.linspace(113,113)
    lat_sc1 = np.linspace(30,41)
    lon_sc2 = np.linspace(113,122)
    lat_sc2 = np.linspace(41,41)
    lon_sc3 = np.linspace(122,122)
    lat_sc3 = np.linspace(30,41)
    lon_sc4 = np.linspace(113,122)
    lat_sc4 = np.linspace(30,30)

    x1,y1 = m(lon_sc1,lat_sc1)
    m.plot (x1, y1,  linewidth=1.0 ,color='black')
    x1,y1 = m(lon_sc2,lat_sc2)
    m.plot (x1, y1,  linewidth=1.0, color='black')
    x1,y1 = m(lon_sc3,lat_sc3)
    m.plot (x1, y1,  linewidth=1.0, color='black')
    x1,y1 = m(lon_sc4,lat_sc4)
    m.plot (x1, y1,  linewidth=1.0, color='black')
    ax.annotate('NC', xy=(0.62, 0.823), xycoords='axes fraction',fontsize=10,color='blue',weight='bold')
 
    

def contour_aod (filename):
    fig = plt.figure(figsize=(8,4))
    aod = np.mean(aod_com_ann,axis=1)
    
    ### Sulfate aerosols
    ax=plt.subplot(1,2,1)
    m= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,
               llcrnrlat = lat[0,0], urcrnrlat = lat[-1,-1],
               llcrnrlon = lon[0,0], urcrnrlon = lon[-1,-1],resolution='l',area_thresh=10000.)
    m.drawcoastlines(linewidth=0.3)
    x,y = m(lon,lat)
    im=m.contourf(x,y,aod[2,:,:],levels2,norm=norm2,cmap=newcmp1)
    m.drawmeridians([100,110,120], labels=[0,0,0,1], color='grey', linewidth=1.2,dashes =[1,1],fontsize=10)
    m.drawparallels([20,30,40],labels=[1,0,0,0],color='grey', linewidth=1.2,dashes= [1,1],fontsize=10)
    subreg (ax, m)

    
    plt.title('(a) SU ',fontsize=12, loc ='left')
    
    #cb_ax=fig.add_axes([0.08, 0.11, 0.38, 0.025])
    #plt.colorbar(im, cax=cb_ax,pad=0.05,orientation = 'horizontal',shrink = 0.6, ticks=[0,0.2,0.25,0.3,0.35,0.4])
    #cb.ax.tick_params(labelsize=12)
    
    ### black carbon aerosol
    ax= plt.subplot(1,2,2)
    m= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,
               llcrnrlat = lat[0,0], urcrnrlat = lat[-1,-1],
               llcrnrlon = lon[0,0], urcrnrlon = lon[-1,-1],resolution='l',area_thresh=10000.)
    m.drawcoastlines(linewidth=0.3)
    x,y = m(lon,lat)
    im=m.contourf(x,y,aod[0,:,:],levels2,norm=norm2,cmap=newcmp1)
    m.drawmeridians([100,110,120], labels=[0,0,0,1], color='grey', linewidth=1.2,dashes =[1,1],fontsize=10)
    m.drawparallels([20,30,40],labels=[0,0,0,0],color='grey', linewidth=1.2,dashes= [1,1],fontsize=10)
    subreg (ax, m)
    plt.title('(b) BC ',fontsize=12, loc ='left')
    
    cb_ax=fig.add_axes([0.15, 0.1, 0.7, 0.025])
    plt.colorbar(im, cax=cb_ax,pad=0.05,orientation = 'horizontal',shrink = 0.6, ticks=[0,0.01,0.02,0.03,0.04,0.05,0.10,0.15,0.2,0.25,0.3,0.35])
    #cb.ax.tick_params(labelsize=12)
    
    
    ### black carbon aerosol
    #ax= plt.subplot(1,3,3)
    #m= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,
    #           llcrnrlat = lat[0,0], urcrnrlat = lat[-1,-1],
    #           llcrnrlon = lon[0,0], urcrnrlon = lon[-1,-1],resolution='l',area_thresh=10000.)
    #m.drawcoastlines(linewidth=0.3)
    #x,y = m(lon,lat)
    #im=m.contourf(x,y,hursmean,levels1,norm=norm1,cmap= mask_sea_land.prcolor)
    #m.drawmeridians([100,110,120], labels=[0,0,0,1], color='grey', linewidth=1.2,dashes =[1,1],fontsize=10)
    #m.drawparallels([20,30,40],labels=[0,0,0,0],color='grey', linewidth=1.2,dashes= [1,1],fontsize=10)
    #subreg (ax, m)
    #plt.title('(C) Near-surface RH ',fontsize=12, loc ='left')
    
    

    
    fig.subplots_adjust(hspace=0.08,wspace=0.2, top =0.95, bottom =0.12, left =0.08, right =0.94)
    plt.savefig(filename,dpi =500)
      
        
contour_aod ('Contour_plots_AeroCom_components_EAS004_2000_BC_SU.png')

print ('the program done')
    
