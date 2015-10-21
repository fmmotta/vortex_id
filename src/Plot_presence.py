import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
from scipy.spatial import Voronoi, voronoi_plot_2d

RawData =np.loadtxt("data/presentFOAMsw.txt")
X,Y,Z	=np.transpose(RawData)

xmax = np.max(X); xmin = np.min(X); Dx = xmax-xmin; xmed=(xmax+xmin)/2;
ymax = np.max(Y); ymin = np.min(Y); Dy = ymax-ymin; ymed=(ymax+ymin)/2;
zmax = np.max(Z); zmin = np.min(Z); Dz = zmax-zmin; zmed=(zmax+zmin)/2;
per = 0.20

print str(xmin)+" "+str(xmax)
print str(ymin)+" "+str(ymax)
print str(zmin)+" "+str(zmax)

#Data 	= [ [x,y,z] for [x,y,z] in RawData if x<xmax and x> xmin and y<ymax and y>ymin]
#X,Y,Z	= np.transpose(Data)
#Z	= np.log(Z)
Nsample= 150

xi = np.linspace(X.min(),X.max(),Nsample)
yi = np.linspace(Y.min(),Y.max(),Nsample)
zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')

fig = plt.figure()


plt.contourf(xi, yi, zi ,Nsample)                             
#plt.pcolormesh(xi,yi,zi)
#plt.imshow(zi, extent = (xmin, xmax, ymin, ymax))
#plt.clim(0,zmax)
cb=plt.colorbar() #aspect='equal',ticks=[0.,11.56,23.11,34.67,46.228])
plt.clim(zmin,zmax)    
#plt.surface().set_clim([zmin,zmax])    
#cb.set_clim(0, zmax)

#for [x,y,z] in Data:
#	s=plt.scatter(x,y,color='black',s=0.01)

#plt.tight_layout()
plt.axes().set_aspect('equal')
plt.savefig('data/presentFOAMsw.pdf',bbox_inches='tight')#,dpi=600)#,bbox_inches='tight')

#plt.show()
