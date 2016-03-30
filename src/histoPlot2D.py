import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.patches as mpatches

n = 0

RawData =np.loadtxt("totalVortices.dat")
# X,Y,Ux,Uy,SigUx,SigUy,W,SigW,Gamma,Beta=np.transpose(RawData)
G,rc,a,b,M00,M01,M10,M11,L1,L2,G00,G01,G10,G11 = np.transpose(RawData)

xmax = np.max(X); xmin = np.min(X); Dx = xmax-xmin; xmed=(xmax+xmin)/2;
ymax = np.max(Y); ymin = np.min(Y); Dy = ymax-ymin; ymed=(ymax+ymin)/2;
#zmax = np.max(Ux); zmin = np.min(Ux); Dz = zmax-zmin; zmed=(zmax+zmin)/2;
#per = 0.20

#print str(xmin)+" "+str(xmax)
#print str(ymin)+" "+str(ymax)
#print str(zmin)+" "+str(zmax)

#Xsample = 200
#Ysample = 200

#xi = np.linspace(X.min(),X.max(),Xsample)
#yi = np.linspace(Y.min(),Y.max(),Ysample)
#zi = griddata((X, Y), Ux, (xi[None,:], yi[:,None]), method='cubic')

#fig = plt.figure()

#plt.imshow(zi, extent=(0,2*math.pi,0,1), origin='lower',cmap='nipy_spectral')

#cb=plt.colorbar(aspect='equal')
#cb.ax.tick_params(labelsize=8) 
#plt.clim(zmin,zmax)
#ax = fig.add_subplot(1,1,1)

plt.hist2d(b,G,bins=100)

plt.savefig("GammaHeightJhisto.pdf",bbox_inches='tight',transparent=True)
