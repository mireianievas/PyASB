#!/usr/bin/env python

'''
SkyMap module

Auxiliary functions to plot the SkyMap
____________________________

This module is part of the PyASB project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyASB project"
__credits__ = ["Miguel Nievas"]
__license__ = "GNU GPL v3"
__shortname__ = "PyASB"
__longname__ = "Python All-Sky Brightness pipeline"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


try:
	import sys
	from astrometry import *
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	import matplotlib.colors as mpc
	import matplotlib.patches as mpp
	import matplotlib as mpl
except:
	print(str(sys.argv[0]) + ': One or more modules missing: matplotlib')
	raise SystemExit

class SkyMap():
	'''	SkyMap class '''
	
	def __init__(self,StarCatalog,ImageInfo,FitsImage):
		if ImageInfo.skymap_path==False:
			# Don't draw anything
			print('Skipping Skymap Graph ...')
			return(None)
		else:
			print('Star Map plot ...')
		
		stretched_fits_data = self.stretch_logdata(FitsImage.fits_data,10,99)
		self.define_skymap();
		self.draw_skymap_data(StarCatalog,ImageInfo,stretched_fits_data)
		self.draw_polar_axes(ImageInfo)
		for Star in StarCatalog.StarList:
			self.annotate_skymap(Star)
		#plt.show(self.skyimage)
		self.show_figure(ImageInfo)
	
	def stretch_logdata(self,fits_data,pmin,pmax):
		log_fits_data = np.log(fits_data-np.min(fits_data)+1,dtype="float32")
		valuemin = np.percentile(log_fits_data,pmin)
		valuemax = np.percentile(log_fits_data,pmax)
		stretched_fits_data = log_fits_data.clip(valuemin,valuemax)
		return stretched_fits_data
	
	def define_skymap(self):
		''' Create figure and self.skyimage subplot. '''
		self.skyfigure = plt.figure(figsize=(10,10))
		self.skyimage  = self.skyfigure.add_subplot(111)
	
	def draw_skymap_data(self,StarCatalog,ImageInfo,fits_data):
		''' Draw identified stars with image as background'''
		xpoints_catalog = [Star.Xcoord for Star in StarCatalog.StarList_woPhot]
		ypoints_catalog = [Star.Ycoord for Star in StarCatalog.StarList_woPhot]
		
		xpoints_phot = [Star.Xcoord for Star in StarCatalog.StarList]
		ypoints_phot = [Star.Ycoord for Star in StarCatalog.StarList]
		
		self.skyimage.imshow(fits_data,cmap=mpl.cm.gray)
		
		self.skyimage.scatter(xpoints_catalog,ypoints_catalog,\
			marker='.',c='g',alpha=0.2,label='Star in the catalog')
		
		self.skyimage.scatter(xpoints_phot,ypoints_phot,\
			marker='x',c='r',alpha=0.2,label='Stars found')
		
		self.skyimage.axis([0,ImageInfo.resolution[0],0,ImageInfo.resolution[1]])
		information=str(ImageInfo.date_string)+" UTC\n"+str(ImageInfo.latitude)+5*" "+\
			str(ImageInfo.longitude)+"\n"+ImageInfo.used_filter
		self.skyimage.text(0.005,0.005,information,fontsize='small',color='white',\
			transform = self.skyimage.transAxes,backgroundcolor=(0,0,0,0.75))
		self.skyimage.legend(('In catalog','Detected'),'upper right')
	
	def draw_polar_axes(self,ImageInfo):
		''' Draws meridian and altitude isolines. '''
		
		zenith_xy = zenith_position(ImageInfo)
		
		for each_altitude in np.arange(0,90,15):
			coord_altitude_0 = horiz2xy(0,each_altitude,ImageInfo)
			radius = math.sqrt(\
				(coord_altitude_0[0]-zenith_xy[0])**2 +\
				(coord_altitude_0[1]-zenith_xy[1])**2)
			self.skyimage.add_patch(\
				mpp.Circle((zenith_xy[0],zenith_xy[1]),radius,
				facecolor='k',fill=False, alpha=0.2,label='_nolegend_'))
			self.skyimage.annotate(\
				str(each_altitude),
				xy=(radius+zenith_xy[0],zenith_xy[1]),
				alpha=0.2,
				fontsize=10)
		
		key_azimuths = {0: "N",90: "E", 180: "S", 270: "W"}
		
		for each_azimuth in np.arange(0,360,30):
			coord_azimuth_0 = horiz2xy(each_azimuth,0,ImageInfo)
			self.skyimage.plot(\
				[zenith_xy[0],coord_azimuth_0[0]],
				[zenith_xy[1],coord_azimuth_0[1]],
				color='k',
				alpha=0.2,)
			
			if each_azimuth in key_azimuths:
				azimuth_label = str(key_azimuths[each_azimuth])
			else:
				azimuth_label = str(each_azimuth)
			self.skyimage.annotate(\
				azimuth_label,
				xy=horiz2xy(each_azimuth,ImageInfo.min_altitude,ImageInfo),
				color='k',
				alpha=0.2,
				fontsize=10)
		
	def annotate_skymap(self,Star):
		# Draw identified stars and measuring circ<les.
		# Annotate HD catalog code and Magnitude for each star.
		self.skyimage.scatter(Star.Xcoord,Star.Ycoord,marker='x',c='r',alpha=0.2,label='Identified stars')
		self.skyimage.add_patch(mpp.Circle((Star.Xcoord,Star.Ycoord),Star.R1,facecolor='none',edgecolor=(0,0,0.8),\
			linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
		self.skyimage.add_patch(mpp.Circle((Star.Xcoord,Star.Ycoord),Star.R2,facecolor='none',edgecolor=(0,0.8,0),\
			linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
		self.skyimage.add_patch(mpp.Circle((Star.Xcoord,Star.Ycoord),Star.R3,facecolor='none',edgecolor=(0.8,0,0),\
			linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
		self.skyimage.annotate(Star.name,xy=(Star.Xcoord,Star.Ycoord), xycoords='data',xytext=(0, 3),\
			textcoords='offset points',fontsize=8)
		self.skyimage.annotate(Star.FilterMag,xy=(Star.Xcoord,Star.Ycoord), xycoords='data',xytext=(0,-10),\
			textcoords='offset points',fontsize=8)
			
	def show_figure(self,ImageInfo):
		def skymap_filename(ImageInfo):
			filename = ImageInfo.skymap_path +\
				"/Skymap_"+ImageInfo.obs_name+"_"+ImageInfo.fits_date+"_"+\
				ImageInfo.used_filter+".png"
			return(filename)
		
		if ImageInfo.skymap_path=="screen":
			plt.show()
		else:
			plt.savefig(skymap_filename(ImageInfo),bbox_inches='tight')
		
		plt.clf()
		plt.close('all')
