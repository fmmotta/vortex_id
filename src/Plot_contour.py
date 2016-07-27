import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.patches as mpatches

#RawData =np.loadtxt("data/initFOAMsw.txt")
RawData =np.loadtxt("data/initUSplit-3.txt")
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
Xsample = 200#256
Ysample = 200#96

xi = np.linspace(X.min(),X.max(),Xsample)
yi = np.linspace(Y.min(),Y.max(),Ysample)
zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')

fig = plt.figure()#(num=None, figsize=(32, 24), dpi=300, facecolor='w', edgecolor='k')

#plt.contourf(xi, yi, zi ,Nsample)                             
#plt.imshow(zi, extent=(0,2.*math.pi,0,1), origin='lower',cmap='nipy_spectral')

#plt.imshow(zi, extent=(0,2*math.pi,0,1), origin='lower',cmap='nipy_spectral')

plt.imshow(zi, extent=(-5,4.95,-5,4.95), origin='lower',cmap='nipy_spectral')

#plt.pcolormesh(xi,yi,zi)
#plt.imshow(zi, extent = (xmin, xmax, ymin, ymax))
#plt.clim(0,zmax)
cb=plt.colorbar(aspect='equal')#,ticks=[0.,11.56,23.11,34.67,46.228])
cb.ax.tick_params(labelsize=8) 
plt.clim(zmin,zmax)
ax = fig.add_subplot(1,1,1)

#fig.patches.append(mpatches.Circle([0.445, 0.45], 0.03, transform=fig.transFigure,facecolor="none",linestyle="dashed",edgecolor="r",alpha=1))

#plt.surface().set_clim([zmin,zmax])    
#cb.set_clim(0, zmax)

#for [x,y,z] in Data:
#	s=plt.scatter(x,y,color='black',s=0.01)

#plt.tight_layout()
#plt.axes().set_aspect('equal')

RawData =np.loadtxt("data/vortices.dat")
W,R,X,Y =np.transpose(RawData)

plt.scatter(X,Y,s=5, lw = 0,color='green')

plt.savefig('data/initFOAMsw.pdf',bbox_inches='tight',transparent=True)#,dpi=600)#,bbox_inches='tight')

#plt.show()
