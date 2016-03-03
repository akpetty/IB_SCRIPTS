import matplotlib
matplotlib.use("AGG")
# basemap import
from mpl_toolkits.basemap import Basemap, shiftgrid
# Numpy import
from osgeo import osr, gdal
import numpy as np
import mpl_toolkits.basemap.pyproj as pyproj
from pylab import *
import sys
sys.path.append('/Users/apetty/GoogleDrive/GPYTHON/RIDGE_SCRIPTS/')
sys.path.append('/Users/apetty/GoogleDrive/GPYTHON/')
import ridge_objects as ro
import alek_objects as apy
import numpy.ma as ma
from matplotlib.patches import Polygon
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import mark_inset
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
#import cubehelix
from netCDF4 import Dataset
from matplotlib import rc
import itertools
from scipy.spatial import cKDTree as KDTree
from scipy import stats
import scipy.optimize as optimization

from glob import glob
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9
rcParams['legend.fontsize']=9
rcParams['font.size']=8
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
my_cmap=apy.perceptual_colormap("Linear_L", reverse=1)

m = Basemap(projection='npstere',boundinglat=66,lon_0=0, resolution='l'  )

def plot_icebridge(xpts, ypts, thickness, year):

	fig = figure(figsize=(4.25,5))
	ax1 = fig.add_axes([0.0, 0.15, 1.0, 0.85])
	minval = 0.
	maxval = 5.

	#im1 = m.scatter(xpts, ypts, c=thickness, edgecolor='none')
	im1=hexbin(xpts, ypts, C=thickness, gridsize=60, vmin=minval, vmax=maxval, cmap=my_cmap)
	#im1=m.scatter(xpts, ypts, c=thickness, s=50, edgecolor='none', vmin=minval, vmax=maxval, cmap=cm.jet)

	#m.drawmapboundary(fill_color='0.4' , zorder=1)
	#m.fillcontinents(color='grey',lake_color='grey', zorder=5)
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

	savefig('/Users/apetty/NOAA/FIGURES/ICEBRIDGE/icebridge_data'+str(year)+'M3.png', dpi=300)

def plot_regions(xptsIB, yptsIB, region_maskR,year):

	fig = figure(figsize=(4.25,5))
	ax1 = fig.add_axes([0.0, 0.15, 1.0, 0.85])
	minval = 0.
	maxval = 11.

	#im1 = m.scatter(xpts, ypts, c=thickness, edgecolor='none')
	im1=hexbin(xptsIB, yptsIB, C=region_maskR, gridsize=60, vmin=minval, vmax=maxval, cmap=my_cmap)
	#im1=m.scatter(xpts, ypts, c=thickness, s=50, edgecolor='none', vmin=minval, vmax=maxval, cmap=cm.jet)

	#m.drawmapboundary(fill_color='0.4' , zorder=1)
	#m.fillcontinents(color='grey',lake_color='grey', zorder=5)
	m.drawcoastlines(linewidth=0.5)
	m.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
	m.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)
	xS, yS = m(230, 62)
	ax1.text(xS, yS, year,fontsize=12, zorder = 11)



	cax = fig.add_axes([0.1, 0.1, 0.8, 0.04])
	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='max',use_gridspec=True)
	#cbar.set_alpha(1)
	#cbar.draw_all()
	cbar.set_label('Region (m)')
	xticks = np.linspace(minval, maxval, 5)
	#cbar.set_ticks(xticks)

	savefig('/Users/apetty/NOAA/FIGURES/ICEBRIDGE/region_grid'+str(year)+'.png', dpi=300)



def read_icebridgeALL(year, mask_hi=1, mask_fboard=1, mask_nonarctic=1):
	#print mask_hi, mask_fboard
	lats_total=[] 
	lons_total=[]
	thickness_total=[]
	snow_thickness_total=[]
	fboard_unc_total=[]
	if (year>2013):
		files = glob('/Users/apetty/NOAA/DATA/ICEBRIDGE/ICEBRIDGE_HI/QLOOK/'+str(year)+'*/*.txt')
	else:
		files = glob('/Users/apetty/NOAA/DATA/ICEBRIDGE/ICEBRIDGE_HI/IDCSI4/'+str(year)+'.*/*.txt')
	
	for x in xrange(size(files)):
		#print x,'/',size(files)
		data = genfromtxt(files[x], delimiter=',', skip_header=1, dtype=str)
		# data is a table-like structure (a numpy recarray) in which you can access columns and rows easily
		lats = data[:, 0].astype(float)
		lons = data[:, 1].astype(float)
		lons[where(lons>180)]=-(360-lons)
		thickness = data[:, 2].astype(float)
		snow_thickness = data[:, 7].astype(float)
		fboard_unc = data[:, 6].astype(float)

		lats_total.extend(lats)
		lons_total.extend(lons)
		thickness_total.extend(thickness)
		snow_thickness_total.extend(snow_thickness)
		fboard_unc_total.extend(fboard_unc)

	thickness_total=array(thickness_total)
	snow_thickness_total=array(snow_thickness_total)
	fboard_unc_total=array(fboard_unc_total)
	lats_total=array(lats_total)
	lons_total=array(lons_total)

	a=size(thickness_total)
	#print a
	if (mask_hi==1):
		good_data=where((thickness_total>=0.)&(thickness_total<=20.))
		thickness_total = thickness_total[good_data]
		snow_thickness_total=snow_thickness_total[good_data]
		fboard_unc_total=fboard_unc_total[good_data]
		lats_total = lats_total[good_data]
		lons_total = lons_total[good_data]

	if (mask_fboard==1):
		good_data=where((fboard_unc_total<=0.1))
		thickness_total = thickness_total[good_data]
		snow_thickness_total=snow_thickness_total[good_data]
		fboard_unc_total=fboard_unc_total[good_data]
		lats_total = lats_total[good_data]
		lons_total = lons_total[good_data]

	if (mask_nonarctic==1):
		xptsIB, yptsIB = m(lons_total, lats_total)
		region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsIB, yptsIB), method='nearest')
		good_data = where((region_maskR==8))
		lats_total = lats_total[good_data]
		lons_total=lons_total[good_data]
		thickness_total=thickness_total[good_data]
		snow_thickness_total=snow_thickness_total[good_data]
		fboard_unc_total=fboard_unc_total[good_data]
		plot_regions(xptsIB, yptsIB, region_maskR, year)

	b = size(thickness_total)
	#print b
	ratio=float(b)/float(a)
	#print ratio
	return lats_total, lons_total, thickness_total, snow_thickness_total, fboard_unc_total


region_mask, xptsM, yptsM = ro.get_region_mask(m)

region = 0

#REMOVE STUFF MAYBE OVER LAND
if (region==0): 
	region_lonlat = [-150, 10, 81, 90]
	region_str='CA'
if (region==1): 
	region_lonlat = [-170, -120, 69, 79]
	region_str='BC'


for year in xrange(2009, 2015):
	latsIB,lonsIB,thicknessIB, snowIB, fboarduncIB= read_icebridgeALL(year, mask_hi=1, mask_fboard=0, mask_nonarctic=1)
	
	mask = where((lonsIB>region_lonlat[0]) & (lonsIB<region_lonlat[1]) & (latsIB>region_lonlat[2]) & (latsIB<region_lonlat[3]))
	mean_thickness = mean(thicknessIB[mask])
	print mean_thickness

	xpts, ypts = m(lonsIB, latsIB)
	plot_icebridge(xpts, ypts, thicknessIB, year)
	

