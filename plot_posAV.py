############################################################## 
# Date: 20/02/16
# Name: plot_posAV.py
# Author: Alek Petty
# Description: Script to plot posAV info of sea ice flights
# Input requirements: posAV data (downscaled)

import matplotlib
matplotlib.use("AGG")
# basemap import
from mpl_toolkits.basemap import Basemap, shiftgrid
# Numpy import
import numpy as np
import mpl_toolkits.basemap.pyproj as pyproj
from pylab import *
import IB_functions as ro
import numpy.ma as ma
from scipy.interpolate import griddata
import os
from glob import glob
#mpl.cm.register_cmap(name='cubehelix1', data=mpl._cm.cubehelix(gamma=1.0, s=1.5, r=-1.0, h=1.5))
#my_cmap=apy.perceptual_colormap("Linear_L", reverse=1)
viridis=ro.get_new_cmaps(cmap_str='viridis')

def get_pos_AV(datapath, year, pole_str='GR'):
	posAVpath = datapath+'/ICEBRIDGE/POSAV/SEA_ICE/'+pole_str+'/'+str(year)+'_'+pole_str+'_NASA/'
	#files = glob()
	files = glob(posAVpath+'/*.txt')
	print year, size(files)
	time=[]
	lats=[]
	lons=[]
	alt=[]
	vel=[]
	pitch=[]
	roll=[]
	flines=[]

	for i in xrange(size(files)):
		dataT = loadtxt(files[i], skiprows=1)
		flineT=i*np.ones(dataT[:, 6].shape)
	    # [mean(xS), mean(yS), ice_area, size(label_nums_ma)-1, ridge_area_all, ridge_area_big, mean_ridge_height_all, mean_ridge_height_big]
		time.extend(dataT[:, 0].astype(float))
		lats.extend(dataT[:, 1].astype(float))
		lons.extend(dataT[:, 2].astype(float))
		alt.extend(dataT[:, 3].astype(float))
		vel.extend(dataT[:, 4].astype(float))
		pitch.extend(dataT[:, 5].astype(float))
		roll.extend(dataT[:, 6].astype(float))
		flines.extend(flineT)
	xpts, ypts=m(lons, lats)
	return xpts, ypts, alt, vel, pitch, roll, flines

rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize']=8
rcParams['ytick.labelsize']=8
rcParams['font.size']=9
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

#m = Basemap(projection='npstere',boundinglat=68,lon_0=0, resolution='l'  )
m=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)
#-------------- ATM AND DMS PATHS------------------

datapath = '../../DATA/'
figpath = './Figures/'



xpts=[]
ypts=[]
alt=[]
vel=[]
pitch=[]
roll=[]
flines=[]
pole_str='GR'
start_year=2009
end_year=2015
num_years=end_year-start_year+1
for year in xrange(start_year, end_year+1):
	#BIG IS IN TERMS OF AREA - WHERE WE REMOVE RIDGES LESS THAN 100m2
	xptsT, yptsT, altT, velT, pitchT, rollT, flinesT= get_pos_AV(datapath, year, pole_str=pole_str)
	xpts.append(xptsT)
	ypts.append(yptsT)
	alt.append(altT)
	vel.append(velT)
	pitch.append(pitchT)
	roll.append(rollT)
	flines.append(flinesT)

option=4
if (option==1):
	minval = 200
	maxval = 800
	var = alt
	out='Alt'
	label='Plane altitude (m)'
if (option==2):
	minval = -5
	maxval = 5
	var = pitch
	out='Pitch'
	label='Pitch (deg)'
if (option==3):
	minval = -5
	maxval = 5
	var = roll
	out='Roll'
	label='Roll (deg)'
if (option==4):
	minval = 0
	maxval = max([max(flines[x]) for x in xrange(num_years)])
	var = flines
	out='Flightlines'
	label='Flight line (num)'

axesname = ['ax1', 'ax2', 'ax3', 'ax4', 'ax5', 'ax6']
aspect = m.ymax/m.xmax
textwidth=5.
fig = figure(figsize=(textwidth,(textwidth*(3./2.)*aspect)+1.2))
subplots_adjust(left = 0.01, right = 0.99, bottom=0.02, top = 0.99, wspace = 0.02, hspace=0.01)
for i in xrange(7):
	vars()['ax'+str(i)] = subplot(4, 2, i+1)

	im1 = vars()['ax'+str(i)].hexbin(xpts[i], ypts[i], C=var[i], gridsize=100, vmin=minval, vmax=maxval, cmap=viridis, zorder=2, rasterized=True)

	#mplot.fillcontinents(color='w',lake_color='grey', zorder=3)
	m.drawcoastlines(linewidth=0.25, zorder=5)
	m.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	m.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)


	xS, yS = m(177, 64.2)
	vars()['ax'+str(i)].text(xS, yS, str(i+2009), zorder = 11)

cax = fig.add_axes([0.52, 0.2, 0.45, 0.03])
cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both',use_gridspec=True)
#cbar.set_alpha(1)
#cbar.draw_all()
cbar.set_label(label, labelpad=1)
xticks = np.linspace(minval, maxval, 6)
cbar.set_ticks(xticks)
cbar.solids.set_rasterized(True)
savefig(figpath+'POSAV'+out+pole_str+'.png', dpi=300)



