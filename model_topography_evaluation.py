from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import Validation as Val
from matplotlib.ticker import FixedFormatter, FixedLocator
from matplotlib.ticker import (MultipleLocator)

### reference: https://www.dkrz.de/up/de-services/de-analysis/de-vis/vis-sw/de-pyngl-pynio
### http://www.pyngl.ucar.edu/index.shtml
path = '/home/EAS-aerosol/COORDINATE/'
file1 = 'HSURF_878x590_EAS11.nc'
file2 = 'coordinates_eas004.nc'

data_50 = Dataset(path+file1)
lon= data_50.variables['rlon']
lat= data_50.variables['rlat']
lon_50 = data_50.variables['lon'][:,:]
lat_50 = data_50.variables['lat'][:,:]
topo_50= data_50.variables['HSURF'][:,:]
rotpole = data_50.variables['rotated_pole']


#### lon_12km
data=Dataset(path+'HSURF_EAS11.nc')
lon12 = data.variables['lon']
lat12 = data.variables['lat']

data = Dataset(path+file2)

lon_004 = data.variables['lon'][:,:]
lat_004 = data.variables['lat'][:,:]
topo_004= data.variables['HSURF'][0,:,:]


def normalize180(lon):
    """Normalize lon to range [180, 180)"""
    lower = -180.; upper = 180.
    if lon > upper or lon == lower:
        lon = lower + abs(lon + upper) % (abs(lower) + abs(upper))
    if lon < lower or lon == upper:
        lon = upper - abs(lon - lower) % (abs(lower) + abs(upper))
    return lower if lon == upper else lon

lon_0 = normalize180(rotpole.grid_north_pole_longitude-180.)
rotlon = rotpole.grid_north_pole_longitude
rotlat = rotpole.grid_north_pole_latitude

levels=[-250,0,250,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,10000]
norm = matplotlib.colors.BoundaryNorm(levels,len(levels))


blue=cm.get_cmap('terrain',256)
color=blue(np.linspace (0, 1, 256))
white=np.array([255/256, 255/256, 255/256, 1])

new =np.zeros (shape=(14,4))
for i in np.arange(0,4):
    new[i+1,:]=color[60+i*18,:]
    new[i+5,:]=color[130+i*20,:]
    

new[0,:] =np.array([151/256, 182/256, 225/256, 1])
new[9,:]=color[200,:]  
new[10,:]=color[215,:]
new[11,:]=color[230,:] 
new[12,:]=color[246,:]
new[13,:]=color[250,:]

newcmp = ListedColormap(new)


def draw_screen_poly( lat1, lon1, m):
    x1, y1 = m( lon1, lat1 )
    
    xy = zip(x1,y1)
    poly = Polygon( list(xy), facecolor='none',edgecolor='black',linewidth=0.6 )
    plt.gca().add_patch(poly)

### station longitude and latitude, and its number
data1 = np.loadtxt ('/home/EAS-aerosol/Validation/Stations_info_selected_from_166stations.txt')

fig=plt.figure(figsize=(12,7),constrained_layout=False)
ax2= fig.add_axes([0.02, 0.556, 0.44,0.41])
m1= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=rotlon,o_lat_p=rotlat,llcrnrlat = lat_004[0,0], 
           urcrnrlat = lat_004[-1,-1],llcrnrlon = lon_004[0,0], urcrnrlon = lon_004[-1,-1],resolution='l',area_thresh=10000.)
m1.drawcoastlines(linewidth=0.2)
x,y = m1(lon_004,lat_004)
im =m1.contourf(x,y,topo_004,levels,norm=norm,cmap=newcmp)
m1.drawmeridians([100,110,120], labels=[0,0,0,1], color='black', linewidth=0.5,dashes =[3,2], fontsize=12)
m1.drawparallels([20,30,40],labels=[1,0,0,0],color='black', linewidth=0.5,dashes =[3,2],fontsize=12)
### station location
X, Y = m1(data1[:,2]*0.01,data1[:,1]*0.01)
ax2.scatter (X,Y, marker='o', s=8,facecolor='none' ,edgecolor='black')
ax2.axis("off")
plt.title('(a) Analysis domain',fontsize=12, loc ='left',pad=0.1)

### subregions
lon_sc1 = np.linspace(105,105)
lat_sc1 = np.linspace(18,28)
lon_sc2 = np.linspace(105,122)
lat_sc2 = np.linspace(28,28)
lon_sc3 = np.linspace(122,122)
lat_sc3 = np.linspace(18,28)
lon_sc4 = np.linspace(105,122)
lat_sc4 = np.linspace(18,18)

x1,y1 = m1(lon_sc1,lat_sc1)
m1.plot (x1, y1,  linewidth=1.0 ,color='black')
x1,y1 = m1(lon_sc2,lat_sc2)
m1.plot (x1, y1,  linewidth=1.0, color='black')
x1,y1 = m1(lon_sc3,lat_sc3)
m1.plot (x1, y1,  linewidth=1.0, color='black')
x1,y1 = m1(lon_sc4,lat_sc4)
m1.plot (x1, y1,  linewidth=1.0, color='black')
ax2.annotate('SE', xy=(0.37, 0.37), xycoords='axes fraction',fontsize=10,color='blue',weight='bold')

## north China
lon_sc1 = np.linspace(113,113)
lat_sc1 = np.linspace(30,41)
lon_sc2 = np.linspace(113,122)
lat_sc2 = np.linspace(41,41)
lon_sc3 = np.linspace(122,122)
lat_sc3 = np.linspace(30,41)
lon_sc4 = np.linspace(113,122)
lat_sc4 = np.linspace(30,30)

x1,y1 = m1(lon_sc1,lat_sc1)
m1.plot (x1, y1,  linewidth=1.0 ,color='black')
x1,y1 = m1(lon_sc2,lat_sc2)
m1.plot (x1, y1,  linewidth=1.0, color='black')
x1,y1 = m1(lon_sc3,lat_sc3)
m1.plot (x1, y1,  linewidth=1.0, color='black')
x1,y1 = m1(lon_sc4,lat_sc4)
m1.plot (x1, y1,  linewidth=1.0, color='black')
ax2.annotate('NC', xy=(0.62, 0.82), xycoords='axes fraction',fontsize=10,color='blue',weight='bold')




tt = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]

ax=plt.subplot( 2, 2, 3)

l1, =plt.plot(tt,np.nanmean(Val.mod_dcycle[:,:,0],axis=1),color='C0',linewidth=2.0)
l2, =plt.plot(tt,np.nanmean(Val.mod_dcycle[:,:,1],axis=1),color='C3',linewidth=2.0)
l3, =plt.plot(tt,np.nanmean(Val.obscycle,axis=1),color='gray',linewidth=3)
       
plt.grid(linestyle='dashed', linewidth= 0.5)           
plt.xlim(tt[0],tt[10])
plt.ylim(0.1,0.24)
plt.xlabel ('Time (UTC)',fontsize=12)
plt.ylabel ('Precipitation (mm/h)',fontsize=12)
plt.title ('(c) Mean', loc ='left',fontsize=12)
plt.xticks ([0,4,8,12,16,20,24])
#plt.legend ([l3,l1,l2],['OBS','REF_0.11','REF_0.04'], loc ='lower right',frameon=False,facecolor="white",fontsize=14)
ax.tick_params(labelsize=12,width=1.2)
r2= np.corrcoef (np.nanmean(Val.mod_dcycle[:,:,0],axis=1),np.nanmean(Val.obscycle,axis=1))[0,1]**2
r2_04= np.corrcoef (np.nanmean(Val.mod_dcycle[:,:,1],axis=1),np.nanmean(Val.obscycle,axis=1))[0,1]**2

ax.text (12,0.225, 'R$^2$: '+str('%.2f' % round(r2,2)),color= 'C0',fontsize=12)
ax.text (12,0.212, 'R$^2$: '+str('%.2f' % round(r2_04,2)),color= 'C3',fontsize=12)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.2)
    
#plt.legend ([l3,l1,l2],['OBS','REF_0.11','REF_0.04'], loc ='lower right',frameon=False,facecolor="white",fontsize=10)
ax.xaxis.set_major_locator(FixedLocator([0,4,8,12,16,20,24]))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(0.02))
ax.xaxis.set_major_formatter(FixedFormatter(['00','04','08','12','16','20','24']))

ax1=plt.subplot( 2, 2, 4)
plt.plot(tt,np.nanmean(Val.mod_whour[:,:,0],axis=1),color='C0',linewidth=2.0)
plt.plot(tt,np.nanmean(Val.mod_whour[:,:,1],axis=1),color='C3',linewidth=2.0)
plt.plot(tt,np.nanmean(Val.obswhour,axis=1),color='gray',linewidth=3)
       
plt.grid(linestyle='dashed', linewidth= 0.5)           
plt.xlim(tt[0],tt[10])
plt.ylim(0.06,0.2)
plt.xlabel ('Time (UTC)',fontsize=14)
plt.ylabel ('Frequency',fontsize=14)
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(1.2)
r2= np.corrcoef (np.nanmean(Val.mod_whour[:,:,0],axis=1),np.nanmean(Val.obswhour,axis=1))[0,1]**2
r2_04= np.corrcoef (np.nanmean(Val.mod_whour[:,:,1],axis=1),np.nanmean(Val.obswhour,axis=1))[0,1]**2

ax1.text (12,0.185, 'R$^2$: '+str('%.2f' % round(r2,2)),color= 'C0',fontsize=12)
ax1.text (12,0.172, 'R$^2$: '+str('%.2f' % round(r2_04,2)),color= 'C3',fontsize=12)
plt.title ('(d) Wet hour frequency', loc ='left',fontsize=12)
  
plt.xticks ([0,4,8,12,16,20,24])
#plt.legend ([l3,l1,l2],['OBS','REF_0.11','REF_0.04'], loc ='upper right',frameon=False,facecolor="white",fontsize=14)
ax1.xaxis.set_major_locator(FixedLocator([0,4,8,12,16,20,24]))
ax1.xaxis.set_minor_locator(MultipleLocator(1))
ax1.xaxis.set_major_formatter(FixedFormatter(['00','04','08','12','16','20','24']))
ax1.tick_params(labelsize=12, width =1.2)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.2)
#ax1.set_facecolor(np.array([245/256, 250/256, 254/256, 1]))


#ax=fig.add_axes([0.59, 0.12, 0.38,0.38])
ax=plt.subplot( 2, 2, 2)
base1,cumfreq1 = Val.get_cumfreq(Val.obs_prmax*0.1)
base2,cumfreq2 = Val.get_cumfreq(Val.prmax11_sum)
base3,cumfreq3 = Val.get_cumfreq(Val.prmax04_sum)
l3, = plt.plot(base1,cumfreq1,color='gray',linewidth=3.0)
l1, = plt.plot(base2,cumfreq2,color='C0',linewidth=2.0)
l2, = plt.plot(base3,cumfreq3,color='C3',linewidth=2.0)
    
plt.ylim(ymin=10e-4,ymax = 1)
ax.set_yscale('log')
plt.xlim(xmin=0,xmax = 60)
ax.xaxis.set_minor_locator(MultipleLocator(5))
plt.title (' (b) pHmax', loc ='left',fontsize=12)
plt.ylabel ('Cumulative frequency',fontsize=12)
plt.xlabel ('Precipitation (mm/h)',fontsize=12)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.2)
ax.tick_params(labelsize=12, width = 1.2)
#ax.set_facecolor(np.array([245/256, 250/256, 254/256, 1]))
plt.legend ([l3,l1,l2],['OBS','REF_0.11','REF_0.04'], loc ='upper right',frameon=False,facecolor="white",fontsize=12)


"""
ax=fig.add_axes([0.08, 0.06, 0.3,0.3])
plt.scatter (Val.r2_eas11[:,0]**2,Val.r2_eas04[:,0]**2 , marker ='o',s=32, color='C0',edgecolors='C0',linewidth=0.8)
plt.scatter (Val.r2_eas11[:,1]**2,Val.r2_eas04[:,1]**2 , marker ='o',s=32, color='C3',edgecolors='C3',linewidth=0.8)
plt.ylim(-0.1,1)
plt.xlim(-0.1,1)
lims = [ np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
plt.plot(lims, lims, linestyle='--', linewidth=0.8, color='black' )

plt.title ('(e)'+' R$^2$', loc ='left')

plt.ylabel ('R$^2$ (REF_0.04 vs. OBS)',fontsize=14)
plt.xlabel ('R$^2$ (REF_0.11 vs. OBS)',fontsize=14)
"""








cb_ax=fig.add_axes([0.39, 0.56, 0.015, 0.4])
cb=plt.colorbar(im, cax=cb_ax,pad=0.5,orientation = 'vertical',shrink = 0.6,ticks=[0,250,500,1000,1500,2000,2500,3000,3500,4000,4500,5000])
cb.ax.tick_params(labelsize=10)
#cb.ax.set_xlabel('Terrain (m)',fontsize=14)


#cb_ax=fig.add_axes([0.2, 0.65, 0.6, 0.015])
#cb=plt.colorbar(im, cax=cb_ax,pad=0.5,orientation = 'horizontal',shrink = 0.6,ticks=[0,250,500,1000,1500,2000,2500,3000,3500,4000,4500,5000])
#cb.ax.tick_params(labelsize=14)
#cb.ax.set_xlabel('Terrain (m)',fontsize=14)

fig.subplots_adjust(bottom=0.08, top  =0.96, right =0.97, left =0.08,wspace=0.23, hspace=0.28)
plt.savefig('Model_Topography_subregions1.png',dpi =600)



print ('the program done')



