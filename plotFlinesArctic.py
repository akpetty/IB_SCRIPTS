############################################################## 
# Date: 20/01/16
# Name: plot_ridges_bulk.py
# Author: Alek Petty
# Description: Script to plot wind and ice type and IB flight lines
# Input requirements: ERA-I wind data, ice type, and IB flines
# Extra info: check the wind/ice_type/IB flightline functions for more info on where to put the data.

import matplotlib
matplotlib.use("AGG")

# basemap import
from mpl_toolkits.basemap import Basemap, shiftgrid
# Numpy import
import numpy as np
from pylab import *
from scipy.io import netcdf
import numpy.ma as ma
import string
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import mark_inset
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
from scipy import stats
import IB_functions as ro
from matplotlib import rc
from netCDF4 import Dataset
from glob import glob


rcParams['axes.labelsize'] =8
rcParams['xtick.labelsize']=8
rcParams['ytick.labelsize']=8
rcParams['legend.fontsize']=8
rcParams['font.size']=8
rcParams['lines.linewidth'] = .5
rcParams['patch.linewidth'] = .5

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

#mpl.rc('text', usetex=True)
m=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)

rawdatapath='../../../DATA/'
figpath='./Figures/'

startYear=2009
endYear=2016
numYears=endYear-startYear+1

xpts_all=[]
ypts_all=[]
for year in xrange(startYear, endYear+1, 1):
	print year
	lonsT, latsT = ro.calc_icebridge_flights(year,rawdatapath, 'GR')
	xptsT, yptsT=m(lonsT, latsT)
	xpts_all.append(xptsT)
	ypts_all.append(yptsT)

ice_type=[]
for year in xrange(startYear, 2014+1):
	ice_typeT, xpts_type, ypts_type = ro.get_mean_ice_type(m, rawdatapath, year, res=1, extrapolate=0)
	ice_type.append(ice_typeT)
ice_type_meanT=np.mean(ice_type, axis=0)
ice_type_mean=np.copy(ice_type_meanT)
ice_type_mean[where(ice_type_meanT>0.9)]=0.75
ice_type_mean[where((ice_type_meanT<0.9) & (ice_type_meanT>0.6))]=0.5
ice_type_mean[where((ice_type_meanT<0.6) & (ice_type_meanT>0.4))]=0.25
#mask open water
ice_type_mean=ma.masked_where(ice_type_mean<0.25, ice_type_mean)

etopo=Dataset(rawdatapath+'/OTHER/etopo5.nc', 'r')
topo = etopo.variables['topo'][::1, ::1]
topo_lon = etopo.variables['topo_lon'][::1]
topo_lat = etopo.variables['topo_lat'][::1]
xptsTopo, yptsTopo= m(*np.meshgrid(topo_lon, topo_lat))


colors=['#00BFFF','k', 'm', 'b', 'y', 'r', 'g', 'c']

aspect = m.ymax/m.xmax
textwidth=4.5
fig = figure(figsize=(textwidth,(textwidth*aspect)))
ax1 = gca()

#im5 = m.pcolormesh(xptsTopo, yptsTopo, topo, vmin=-5000, vmax=5000, cmap=cm.Greys_r,shading='gouraud', zorder=1)
m.contour(xptsTopo, yptsTopo, topo,levels=[0], colors='k', zorder=2)
im0 = m.pcolormesh(xpts_type , ypts_type, ice_type_mean, 
	edgecolors='white', vmin=0.25, vmax=0.75, cmap=cm.Greys,shading='gouraud', zorder=2)
plts=[]
for x in xrange(numYears):
	vars()['im'+str(x+1)] = m.plot(xpts_all[x] , ypts_all[x], 
		color = colors[x], linewidth=0.5, zorder=4)
	plts=plts+vars()['im'+str(x+1)]

varnames=[str(x) for x in xrange(startYear, endYear+1)]

leg = ax1.legend(plts, varnames, 
	loc=1, ncol=1,columnspacing=0.8, labelspacing=0.5,handletextpad=0.1, borderaxespad=0.05,bbox_to_anchor=(1.15, 0.8), frameon=False)
llines = leg.get_lines()
setp(llines, linewidth=2.0)
leg.set_zorder(20)

m.drawparallels(np.arange(90,-90,-10), 
	labels=[False,False,True,False], fontsize=8,linewidth = 0.25, zorder=5)
m.drawmeridians(np.arange(-180.,180.,30.), 
	linewidth = 0.25, zorder=10)
#m.fillcontinents(color='white',lake_color='white', zorder=2)
#m.drawcoastlines(linewidth=.5, zorder=5)
#m.bluemarble(zorder=1)


cax = fig.add_axes([0.9, 0.9, 0.07, 0.025])
cbar = colorbar(im0,cax=cax, orientation='horizontal', use_gridspec=True)
cbar.set_ticks([0.25, 0.75])
cbar.set_ticklabels(['FY', 'MY'])
cbar.set_clim(0., 1.)

subplots_adjust( bottom=0.04, top=0.92, left = 0.01, right=0.87, hspace=0.05)
savefig(figpath+'FlinesArctic.png', dpi=300)
close(fig)





