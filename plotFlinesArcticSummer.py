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

def calc_icebridge_flights_traj(year, rawdatapath, loc):
    init_path = rawdatapath+'/ICEBRIDGE/POSAV/SEA_ICE/'+loc+'/'+str(year)+'_'+loc+'_NASA/'

    files = glob(init_path+'*.thin')
    lons=[]
    lats=[]

    for x in xrange(size(files)):
        data = loadtxt(files[x], skiprows=1)
        lats_t = data[::10, 3]
        lons_t = data[::10, 4]
        for x in xrange(size(lats_t)):
            lats.append(lats_t[x])
            lons.append(lons_t[x])

    return lons, lats

#mpl.rc('text', usetex=True)
m=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=100, urcrnrlat=88)

rawdatapath='../../../DATA/'
figpath='./Figures/'


xpts_all=[]
ypts_all=[]

#year='2016SUMMER'
#loc='GR'
#init_path = rawdatapath+'/ICEBRIDGE/POSAV/SEA_ICE/'+loc+'/'+str(year)+'_'+loc+'_NASA/'

lonsT, latsT = calc_icebridge_flights_traj('2016SUMMER', rawdatapath, 'GR')
xptsT, yptsT=m(lonsT, latsT)
xpts_all.append(xptsT)
ypts_all.append(yptsT)

etopo=Dataset(rawdatapath+'/OTHER/etopo5.nc', 'r')
topo = etopo.variables['topo'][::1, ::1]
topo_lon = etopo.variables['topo_lon'][::1]
topo_lat = etopo.variables['topo_lat'][::1]
xptsTopo, yptsTopo= m(*np.meshgrid(topo_lon, topo_lat))


colors=['#00BFFF','k', 'm', 'b', 'y', 'r', 'g', 'c']

aspect = m.ymax/m.xmax
textwidth=2.5
fig = figure(figsize=(textwidth,(textwidth*aspect)))
ax1 = gca()

#im5 = m.pcolormesh(xptsTopo, yptsTopo, topo, vmin=-5000, vmax=5000, cmap=cm.Greys_r,shading='gouraud', zorder=1)
m.contour(xptsTopo, yptsTopo, topo,levels=[0], colors='k', zorder=2)

im2 = m.plot(xpts_all[0] , ypts_all[0], 
	color = colors[-1], linewidth=0.5, zorder=4)
plts=im2

varnames=['2016\n(Summer)']

leg = ax1.legend(plts, varnames, 
	loc=1, ncol=1,columnspacing=0.8, labelspacing=0.5,handletextpad=0.1, borderaxespad=0.05,bbox_to_anchor=(0.95, 0.93), frameon=False)
llines = leg.get_lines()
setp(llines, linewidth=2.0)
leg.set_zorder(20)

m.drawparallels(np.arange(90,-90,-10), 
	labels=[False,False,False,False], fontsize=8,linewidth = 0.25, zorder=5)
m.drawmeridians(np.arange(-180.,180.,30.), 
	linewidth = 0.25, zorder=10)

subplots_adjust( bottom=0.04, top=0.98, left = 0.01, right=0.98, hspace=0.05)
savefig(figpath+'FlinesArcticSummer.png', dpi=300)
close(fig)





