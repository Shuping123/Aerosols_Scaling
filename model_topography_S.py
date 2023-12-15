from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
from matplotlib.patches import ConnectionPatch

### reference: https://www.dkrz.de/up/de-services/de-analysis/de-vis/vis-sw/de-pyngl-pynio
### http://www.pyngl.ucar.edu/index.shtml
path = '/home/LiShuPing/ETH_work/EAS-aerosol/COORDINATE/'
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
#print (lon12[-1,0],lat12[-1,0]) ## top left
#print (lon12[-1,-1],lat12[-1,-1]) ## top right
#print (lon12[0,-1],lat12[0,-1]) ## bottom right
#print (lon12[0,0],lat12[0,0])

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
    
#new[0,:]=color[25,:]
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
    poly = Polygon( list(xy), facecolor='none',edgecolor='red',linewidth=0.6 )
    plt.gca().add_patch(poly)

### station longitude and latitude, and its number
data1 = np.loadtxt ('/home/LiShuPing/ETH_work/EAS-aerosol/Validation/Stations_info_selected_from_166stations.txt')


fig=plt.figure(figsize=(8,6),constrained_layout=False)
ax1 = fig.add_axes([0.1,0.2,0.85,0.65])
m= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=rotlon,o_lat_p=rotlat,llcrnrlat = lat_50[0,0], 
           urcrnrlat = lat_50[-1,-1],llcrnrlon = lon_50[0,0], urcrnrlon = lon_50[-1,-1],resolution='l',area_thresh=10000.)
m.drawcoastlines(linewidth=0.2)
x,y = m(lon_50,lat_50)
im =m.contourf(x,y,topo_50,levels,norm=norm,cmap=newcmp)
m.drawmeridians([50,70,90,110,130,150,170], labels=[0,0,0,1], color='black', linewidth=0.5,dashes =[3,2], fontsize=10)
m.drawparallels([0,20,40,60],labels=[1,0,0,0],color='black', linewidth=0.5,dashes =[3,2],fontsize=10)
ax1.axis("off")


x1,y1 = m(lon12[-1,0],lat12[-1,0])  # left lower point
x2,y2 = m(lon12[-1,-1],lat12[-1,-1])  # left upper point
x3,y3 = m(lon12[0,-1],lat12[0,-1])  # right upper point
x4,y4 = m(lon12[0,0],lat12[0,0])  # right lower point
poly = Polygon([(x1,y1),(x2,y2),(x3, y3),(x4,y4)],facecolor='none',edgecolor='blue',linewidth=0.8)
plt.gca().add_patch(poly)

#m.text(145, 8, '12 km', transform=ax1.transAxes, fontsize=4)
plt.annotate('12 km', xy=(0.85, 0.06), xycoords='axes fraction',fontsize=11,color='blue')


#props = dict( facecolor='white', alpha=0.9,pad=0.4,linewidth=0)
#ax1.text(0.52, 0.706, '4.4 km', transform=ax1.transAxes, fontsize=10,
#        verticalalignment='top', bbox=props)
x1,y1 = m(95.24,15.21)  # left lower point
x2,y2 = m(88.90,42.06)  # left upper point
x3,y3 = m(128.45,43.99)  # right upper point
x4,y4 = m(125.53,16.70)  # right lower point
poly = Polygon([(x1,y1),(x2,y2),(x3, y3),(x4,y4)],facecolor='none',edgecolor='black',linewidth=0.8)
plt.gca().add_patch(poly)
plt.annotate('4.4 km', xy=(0.48, 0.24), xycoords='axes fraction',fontsize=11,color='black')

#x1,y1 = m(95.66,15.81)  # left lower point
#x2,y2 = m(89.72,41.67)  # left upper point
#x3,y3 = m(127.69,43.52)  # right upper point
#x4,y4 = m(125.05,17.26)  # right lower point
#poly = Polygon([(x1,y1),(x2,y2),(x3, y3),(x4,y4)],facecolor='none',edgecolor='gray',linestyle='dashed',linewidth=1.2)
#plt.gca().add_patch(poly)


plt.title(' Model domains',fontsize=12, loc ='center',pad=0.1)







cb_ax=fig.add_axes([0.18, 0.11, 0.68, 0.03])
cb=plt.colorbar(im, cax=cb_ax,pad=0.5,orientation = 'horizontal',shrink = 0.6,ticks=[0,250,500,1000,1500,2000,2500,3000,3500,4000,4500,5000])
cb.ax.tick_params(labelsize=10)
cb.ax.set_xlabel('Terrain (m)',fontsize=10)
fig.subplots_adjust(bottom=0.15, top  =0.96, right =0.94, left =0.05)
plt.savefig('Model_Topography_S.png',dpi =600)



print ('the program done')



