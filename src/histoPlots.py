import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from scipy.interpolate import griddata
from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.patches as mpatches

tag = 'neguSub'
tag2 = 'posneguSub'

smLabel0 = 'Input Values'
smLabel1 = 'Output Values without shear'
smLabel2 = 'Output Values with shear'

RawData =np.loadtxt("negNoShear/histoIna-negNoShear.txt")
valueLeft,valueRight,count = np.transpose(RawData)

f, ax = plt.subplots(1)

plt.xlabel('X position',fontsize=14,color='black')
plt.ylabel('Counts',fontsize=14,color='black')

ax.plot((valueLeft+valueRight)/2,count,'r--',label=smLabel0)

RawData2 =np.loadtxt("negNoShear/histoOua-negNoShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'bo',label=smLabel1)

RawData2 =np.loadtxt("negStrongShear/histoOua-negStrongShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'g^',label=smLabel2)

plt.legend(loc=4)

ax.set_xlim(xmin=-10,xmax=10)
ax.set_ylim(ymin=0)
plt.savefig("histoa-"+tag+".pdf",bbox_inches='tight',transparent=True)

###########################################################################

fig = plt.figure()
RawData =np.loadtxt("negNoShear/histoInb-negNoShear.txt")
valueLeft,valueRight,count = np.transpose(RawData)

f, ax = plt.subplots(1)

plt.xlabel('Y position',fontsize=14,color='black')
plt.ylabel('Counts',fontsize=14,color='black')

ax.plot((valueLeft+valueRight)/2,count,'r--',label=smLabel0)

RawData2 =np.loadtxt("negNoShear/histoOub-negNoShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'bo',label=smLabel1)

RawData2 =np.loadtxt("negStrongShear/histoOub-negStrongShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'g^',label=smLabel2)

plt.legend(loc=4)

ax.set_xlim(xmin=-10,xmax=10)
ax.set_ylim(ymin=0)
plt.savefig("histob-"+tag+".pdf",bbox_inches='tight',transparent=True)

###########################################################################

fig = plt.figure()
RawData =np.loadtxt("negNoShear/histoInG-negNoShear.txt")
valueLeft,valueRight,count = np.transpose(RawData)

f, ax = plt.subplots(1)

plt.xlabel('Circulation',fontsize=14,color='black')
plt.ylabel('Counts',fontsize=14,color='black')

ax.plot((valueLeft+valueRight)/2,count,'r--',label=smLabel0)

RawData2 =np.loadtxt("negNoShear/histoOuG-negNoShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'bo',label=smLabel1)

RawData2 =np.loadtxt("negStrongShear/histoOuG-negStrongShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'g^',label=smLabel2)

plt.legend(loc=2)

ax.set_xlim(xmin=-40,xmax=10)
ax.set_ylim(ymin=0)
plt.savefig("histoG-"+tag+".pdf",bbox_inches='tight',transparent=True)

###########################################################################

fig = plt.figure()
RawData =np.loadtxt("negNoShear/histoInRc-negNoShear.txt")
valueLeft,valueRight,count = np.transpose(RawData)

f, ax = plt.subplots(1)

plt.xlabel('Radius',fontsize=14,color='black')
plt.ylabel('Counts',fontsize=14,color='black')

ax.plot((valueLeft+valueRight)/2,count,'r--',label=smLabel0)

RawData2 =np.loadtxt("negNoShear/histoOuRc-negNoShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'bo',label=smLabel1)

RawData2 =np.loadtxt("negStrongShear/histoOuRc-negStrongShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'g^',label=smLabel2)

plt.legend(loc=1)

ax.set_xlim(xmin=0,xmax=2.5)
ax.set_ylim(ymin=0)
plt.savefig("histoRc-"+tag+".pdf",bbox_inches='tight',transparent=True)

###########################################################################

RawData =np.loadtxt("posnegNoShear/histoIna-posnegNoShear.txt")
valueLeft,valueRight,count = np.transpose(RawData)

f, ax = plt.subplots(1)

plt.xlabel('X position',fontsize=14,color='black')
plt.ylabel('Counts',fontsize=14,color='black')

ax.plot((valueLeft+valueRight)/2,count,'r--',label=smLabel0)

RawData2 =np.loadtxt("posnegNoShear/histoOua-posnegNoShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'bo',label=smLabel1)

RawData2 =np.loadtxt("posnegStrongShear/histoOua-posnegStrongShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'g^',label=smLabel2)

plt.legend(loc=4)
ax.set_xlim(xmin=-10,xmax=10)
ax.set_ylim(ymin=0)
plt.savefig("histoa-"+tag2+".pdf",bbox_inches='tight',transparent=True)

###########################################################################

fig = plt.figure()
RawData =np.loadtxt("posnegNoShear/histoInb-posnegNoShear.txt")
valueLeft,valueRight,count = np.transpose(RawData)

f, ax = plt.subplots(1)

plt.xlabel('Y position',fontsize=14,color='black')
plt.ylabel('Counts',fontsize=14,color='black')

ax.plot((valueLeft+valueRight)/2,count,'r--',label=smLabel0)

RawData2 =np.loadtxt("posnegNoShear/histoOub-posnegNoShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'bo',label=smLabel1)

RawData2 =np.loadtxt("posnegStrongShear/histoOub-posnegStrongShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'g^',label=smLabel2)

plt.legend(loc=4)
ax.set_xlim(xmin=-10,xmax=10)
ax.set_ylim(ymin=0)
plt.savefig("histob-"+tag2+".pdf",bbox_inches='tight',transparent=True)

###########################################################################

fig = plt.figure()
RawData =np.loadtxt("posnegNoShear/histoInG-posnegNoShear.txt")
valueLeft,valueRight,count = np.transpose(RawData)

f, ax = plt.subplots(1)

plt.xlabel('Circulation',fontsize=14,color='black')
plt.ylabel('Counts',fontsize=14,color='black')

ax.plot((valueLeft+valueRight)/2,count,'r--',label=smLabel0)

RawData2 =np.loadtxt("posnegNoShear/histoOuG-posnegNoShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'bo',label=smLabel1)

RawData2 =np.loadtxt("posnegStrongShear/histoOuG-posnegStrongShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'g^',label=smLabel2)

plt.legend(loc=2)
ax.set_xlim(xmin=-40,xmax=25)
ax.set_ylim(ymin=0,ymax=35000)
plt.savefig("histoG-"+tag2+".pdf",bbox_inches='tight',transparent=True)

###########################################################################

fig = plt.figure()
RawData =np.loadtxt("posnegNoShear/histoInRc-posnegNoShear.txt")
valueLeft,valueRight,count = np.transpose(RawData)

f, ax = plt.subplots(1)

plt.xlabel('Radius',fontsize=14,color='black')
plt.ylabel('Counts',fontsize=14,color='black')

ax.plot((valueLeft+valueRight)/2,count,'r--',label=smLabel0)

RawData2 =np.loadtxt("posnegNoShear/histoOuRc-posnegNoShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'bo',label=smLabel1)

RawData2 =np.loadtxt("posnegStrongShear/histoOuRc-posnegStrongShear.txt")
valueLeft2,valueRight2,count2 = np.transpose(RawData2)

ax.plot((valueLeft2+valueRight2)/2,count2,'g^',label=smLabel2)

plt.legend(loc=1)
ax.set_xlim(xmin=0,xmax=2.5)
ax.set_ylim(ymin=0)
plt.savefig("histoRc-"+tag2+".pdf",bbox_inches='tight',transparent=True)

###########################################################################