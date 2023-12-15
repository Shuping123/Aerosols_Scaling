import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.patches import Path, PathPatch


blue=cm.get_cmap('Blues',256)
color=blue(np.linspace (0, 1, 256))
red=cm.get_cmap('nipy_spectral',256)
color1=red(np.linspace (0, 1, 256))

white=np.array([255/256, 255/256, 255/256, 1])
new =np.zeros (shape=(17,4))
for i in np.arange(1,7):
    new[i,:]=color[10+i*28]
for i in np.arange(0,10):
    new[i+7,:]=color1[115+i*12,:]
new[0,:]=white
prcolor = ListedColormap(new)
    
### temperature color
blue=cm.get_cmap('hot',256)
color=blue(np.linspace (0, 1, 256))
red=cm.get_cmap('nipy_spectral',256)
color1=red(np.linspace (0, 1, 256))

white=np.array([255/256, 255/256, 255/256, 1])
new =np.zeros (shape=(17,4))
new[0,:]=color1[30,:]
new[8,:]=color[240,:]
for i in np.arange(0,8):
    new[i+9,:]=color[200-i*21]
for i in np.arange(1,4):
    new[i,:]=color1[30+i*10,:]
    new[i+3,:] =color1[90+i*25,:]
new[7,:]=white
tascolor = ListedColormap(new)

### total precipitation difference (color)
blue=cm.get_cmap('BrBG',256)
color=blue(np.linspace (0, 1, 256))
white=np.array([255/256, 255/256, 255/256, 1])
new =np.zeros (shape=(14,4))
for i in np.arange(0,7):
    new[i,:]=color[30+i*12,:]
for i in np.arange(0,7):
    new[i+7,:]=color[135+i*15]

prdiffcolor = ListedColormap(new)

### color for total precipitation with white color
new = np.zeros (shape=(15,4))
for i in np.arange(0,7):
    new[i,:]=color[30+i*12,:]
for i in np.arange(0,7):
    new[i+8,:]=color[135+i*15]
new[7,:]=white

prdiffcolor1 = ListedColormap(new)

### 2-m temperature difference (color)
blue=cm.get_cmap('Blues',150)
color=blue(np.linspace (0, 1, 150))
red=cm.get_cmap('OrRd',150)
color1=red(np.linspace (0, 1, 150))

white=np.array([255/256, 255/256, 255/256, 1])
new =np.zeros (shape=(14,4))
for i in np.arange(0,7):
    new[i,:]=color[140-i*21,:]
for i in np.arange(0,7):
    new[i+7,:]=color1[10+i*20]
#new[7,:]=white
tasdiffcolor = ListedColormap(new)

#### colors for temperature difference with white color 
new = np.zeros (shape=(15,4))
for i in np.arange(0,7):
    new[i,:]=color[140-i*21,:]
for i in np.arange(0,7):
    new[i+8,:]=color1[10+i*20]
new[7,:]=white
tasdiffcolor1 = ListedColormap(new)


def mask_sea (ax,m):
    ##getting the limits of the map:
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    map_edges = np.array([[x0,y0],[x1,y0],[x1,y1],[x0,y1]])
    polys = [p.boundary for p in m.landpolygons]
    ##combining with map edges
    polys = [map_edges]+polys[:]
    ##creating a PathPatch
    codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]] for p in polys]
    polys_lin = [v for p in polys for v in p]
    codes_lin = [c for cs in codes for c in cs]
    path = Path(polys_lin, codes_lin)
    patch = PathPatch(path,facecolor='lightgray', lw=0)
    ##masking the data:
    ax.add_patch(patch)

def normalize180(lon):
    #Normalize lon to range [180, 180)
    lower = -180.; upper = 180.
    if lon > upper or lon == lower:
        lon = lower + abs(lon + upper) % (abs(lower) + abs(upper))
    if lon < lower or lon == upper:
        lon = upper - abs(lon - lower) % (abs(lower) + abs(upper))
    return lower if lon == upper else lon

