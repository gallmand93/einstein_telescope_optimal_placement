import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap


plt.close("all")


r = np.loadtxt("HLE localisation results 20 deg n=30 ntrials=1e5.txt")
#l = r[:,2]

#l1 =l


n = 30  #Size of array
N = n**2     #No of values
longmin = -11
longmax = 45
deltalong = (longmax - longmin)/float(n)
latmin = 35
latmax = 72
deltalat = (latmax - latmin)/float(n)



X = np.zeros([n,n])
Y = np.zeros([n,n])
Z = np.zeros([n,n])
for i in range(n):
    for j in range(n):
        X[i,j] = r[i+n*j,0]
        Y[i,j] = r[i+n*j,1]
        Z[i,j] = r[i+n*j,2]




plt.figure()

map = Basemap(projection="cyl", lat_0=(latmax+latmin)/2, lon_0=(longmax+longmin)/2,
    resolution = 'l', area_thresh = 1000.0,
    llcrnrlon=longmin, llcrnrlat=latmin,
    urcrnrlon=longmax-deltalong, urcrnrlat=latmax-deltalat)


map.drawcoastlines()
map.drawcountries()
map.fillcontinents(color = "none", lake_color = "none")


map.drawmeridians(np.arange(0, 360, 5))
map.drawparallels(np.arange(-90, 90, 5))
#map.drawrivers()
minZ = np.min(Z)
maxZ = np.max(Z)
diffZ = maxZ - minZ
levels = np.linspace(minZ, maxZ, 100)
cs = map.contourf(X,Y,Z, levels = levels, cmap='rainbow', alpha = 1)

cbar = map.colorbar(cs,location='bottom',pad="5%", ticks=[minZ, minZ+0.2*diffZ, minZ+0.4*diffZ, minZ+0.6*diffZ, minZ+0.8*diffZ, maxZ])
plt.title("Contour plot showing the localisation power\nof the Einstein Telescope by location")
#map.fillcontinents(color = "none")
map.drawmapboundary(fill_color = "none")

plt.savefig("20degcontour900 n=30 ntrials=1e5.png")

plt.show()
