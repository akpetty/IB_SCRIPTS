import matplotlib
matplotlib.use("AGG")
import sys
sys.path.append('/Users/apetty/GoogleDrive/GPYTHON/')
import alek_objects as apy
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
from matplotlib import rc
import csv
from numpy import genfromtxt
from glob import glob

def plot_icebridge(xpts, ypts, thickness, year):

	fig = figure(figsize=(4.25,5))
	ax1 = fig.add_axes([0.0, 0.15, 1.0, 0.85])
	minval = 0.
	maxval = 5.

	#im1 = m.scatter(xpts, ypts, c=thickness, edgecolor='none')
	im1=hexbin(xpts, ypts, C=thickness, gridsize=60, vmin=minval, vmax=maxval, cmap=cm.jet)

	#m.drawmapboundary(fill_color='0.4' , zorder=1)
	m.fillcontinents(color='grey',lake_color='grey', zorder=5)
	m.drawcoastlines(linewidth=0.5)
	m.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
	m.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)
	xS, yS = m(230, 62)
	ax1.text(xS, yS, year,fontsize=12, zorder = 11)



	cax = fig.add_axes([0.1, 0.1, 0.8, 0.04])
	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='max',use_gridspec=True)
	#cbar.set_alpha(1)
	#cbar.draw_all()
	cbar.set_label('Ice thickness (m)')
	xticks = np.linspace(minval, maxval, 5)
	#cbar.set_ticks(xticks)

	savefig('/Users/apetty/NOAA/FIGURES/ICEBRIDGE/icebridge_data'+year+'.png', dpi=300)

def plot_icebridge_scatter(xpts, ypts, thickness, year):

	fig = figure(figsize=(4.25,5))
	ax1 = fig.add_axes([0.0, 0.15, 1.0, 0.85])
	minval = 0.
	maxval = 5.

	#im1 = m.scatter(xpts, ypts, c=thickness, edgecolor='none')
	im1=m.scatter(xpts, ypts, c=thickness, s=50, edgecolor='none', vmin=minval, vmax=maxval, cmap=cm.jet)

	#m.drawmapboundary(fill_color='0.4' , zorder=1)
	m.fillcontinents(color='grey',lake_color='grey', zorder=5)
	m.drawcoastlines(linewidth=0.5)
	m.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
	m.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)
	xS, yS = m(230, 62)
	ax1.text(xS, yS, year,fontsize=12, zorder = 11)



	cax = fig.add_axes([0.1, 0.1, 0.8, 0.04])
	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='max',use_gridspec=True)
	#cbar.set_alpha(1)
	#cbar.draw_all()
	cbar.set_label('Ice thickness (m)')
	xticks = np.linspace(minval, maxval, 5)
	#cbar.set_ticks(xticks)

	savefig('/Users/apetty/NOAA/FIGURES/ICEBRIDGE/icebridge_data'+year+'scatter.png', dpi=300)

def read_icebridgeALL(year, mask_hi=1, mask_nonarctic=1):
	lats_total=[] 
	lons_total=[]
	thickness_total=[]
	snow_thickness_total=[]
	if (year>2013):
		files = glob(rawdatapath+'/ICEBRIDGE/ICEBRIDGE_HI/QLOOK/'+str(year)+'*/*.txt')
	else:
		files = glob(rawdatapath+'/ICEBRIDGE/ICEBRIDGE_HI/IDCSI4/'+str(year)+'.*/*.txt')
	
	for x in xrange(size(files)):
		data = genfromtxt(files[x], delimiter=',', skip_header=1, dtype=str)
		# data is a table-like structure (a numpy recarray) in which you can access columns and rows easily
		lats = data[:, 0].astype(float)
		lons = data[:, 1].astype(float)
		thickness = data[:, 2].astype(float)
		snow_thickness = data[:, 7].astype(float)
		lats_total.extend(lats)
		lons_total.extend(lons)
		thickness_total.extend(thickness)
		snow_thickness_total.extend(snow_thickness)

	thickness_total=array(thickness_total)
	snow_thickness_total=array(snow_thickness_total)
	lats_total=array(lats_total)
	lons_total=array(lons_total)

	if (mask_hi==1):
		good_data=where((thickness_total>=0.)&(thickness_total<=20.))
		thickness_total = thickness_total[good_data]
		snow_thickness_total=snow_thickness_total[good_data]
		lats_total = lats_total[good_data]
		lons_total = lons_total[good_data]
	if (mask_nonarctic==1):
		xptsIB, yptsIB = mplot(lons_total, lats_total)
		region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsIB, yptsIB), method='nearest')
		good_data = where((region_maskR==8))
		lats_total = lats_total[good_data]
		lons_total=lons_total[good_data]
		thickness_total=thickness_total[good_data]
		snow_thickness_total=snow_thickness_total[good_data]

	xpts,ypts = mplot(lons_total, lats_total)

	return xpts,ypts, lats_total, lons_total, thickness_total, snow_thickness_total

def read_icebridge(filepath, var_num):

	

    data = genfromtxt(filepath, delimiter=',', dtype=None)
    # data is a table-like structure (a numpy recarray) in which you can access columns and rows easily
    lats = data[1:-1, 0].astype(float)
    lons = data[1:-1, 1].astype(float)
    var = data[1:-1, var_num].astype(float)

    #var = ma.masked_where(var<var_min , var)
    #lats = ma.masked_where(var<var_min , lats)
    #lons = ma.masked_where(var<var_min , lons)
    #var = ma.masked_where(var>var_max , var)
    #lats = ma.masked_where(var>var_max , lats)
    #lons = ma.masked_where(var>var_max , lons)

    return lats, lons, var


mpl.rc("ytick",labelsize=10)
mpl.rc("xtick",labelsize=10)
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
#mpl.rc('text', usetex=True)

#m = Basemap(projection='npstere',boundinglat=66,lon_0=0, resolution='l'  )
m=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)

for x in xrange(2010, 2011):
	year = str(x)

	files = glob('/Users/apetty/NOAA/DATA/ICEBRIDGE/ICEBRIDGE_HI/DATA/'+year+'.*/*.txt')


	lats_total=[] 
	lons_total=[]
	thickness_total=[]

	for x in xrange(size(files)):
		lats, lons, thickness= read_icebridge(files[x], 2)

		lats_total.extend(lats)
		lons_total.extend(lons)
		thickness_total.extend(thickness)
	thickness_total=array(thickness_total)
	lats_total=array(lats_total)
	lons_total=array(lons_total)	
	xpts,ypts = m(lons_total,lats_total)

	#xpts=ma.getdata(xpts)
	#ypts=ma.getdata(ypts)
	# thickness_totalA[where((thickness_totalA>=0.))]
	#thickness_total=thickness_total[where(~(np.isnan(thickness_total)))]

	good_data=where((thickness_total>=0.)&(thickness_total<=20.))
	thickness_total = thickness_total[good_data]
	xpts = xpts[good_data]
	ypts = ypts[good_data]

	#thickness_total=thickness_total[where(~isnan(thickness_total))]
	#xpts=xpts[where(~isnan(thickness_total))]
	#ypts=ypts[where(~isnan(thickness_total))]



plot_icebridge(xpts, ypts, thickness_total, year)
	plot_icebridge_scatter(xpts, ypts, thickness_total, year)

#thickness_total=array(thickness_total)
#thickness_totalA = thickness_total[where((thickness_total>=0)&(thickness_total<8.))]

#a = ma.getdata(thickness_total)
#a=a[where()]



