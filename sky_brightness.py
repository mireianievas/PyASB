#!/usr/bin/env python

'''
Sky Brightness photometry module

Measure Sky brightness from the image using previous 
instrument calibration. 

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
	import numpy as np
	import scipy.interpolate as sint
	import matplotlib as mpl
	import matplotlib.pyplot as plt
	import matplotlib.colors as mpc
	import matplotlib.patches as mpp
except:
	print('One or more modules missing: numpy, matplotlib')
	raise SystemExit

class SkyBrightness():
	'''
	Class with Sky Brightness measure methods.
	Init requires FitsImage, ImageCoordinates, ImageInfo and Bouguerfit objects.
	'''
	
	def __init__(self,FitsImage,ImageInfo,ImageCoordinates,BouguerFit):
		self.measure_in_grid(FitsImage,ImageInfo,ImageCoordinates,BouguerFit)
		self.measure_in_positions()
	
	@staticmethod
	def sky_brightness_region(BouguerFit,ImageInfo,fits_region_values):
		sky_flux = np.median(fits_region_values)
		sky_flux_err = np.std(fits_region_values)
		sky_brightness = BouguerFit.Regression.mean_zeropoint - \
			2.5*np.log10(sky_flux/(ImageInfo.exposure*ImageInfo.pixel_scale))
		sky_brightness_err = np.sqrt(BouguerFit.Regression.error_zeropoint**2 +\
			(2.5*sky_flux_err/(np.log(10)*sky_flux))**2)
		
		return(sky_brightness,sky_brightness_err)
	
	def measure_in_grid(self,FitsImage,ImageInfo,ImageCoordinates,BouguerFit):
		''' Returns sky brightness measures in a grid with a given separation 
		    in degrees and interpolates the result with griddata.'''
		azseparation  = 30
		altseparation = 15
		
		AZdirs  = np.arange(0,360+1,azseparation)
		ZDdirs = np.arange(0,90+1,altseparation)
		
		self.AZgrid,self.ZDgrid = np.meshgrid(AZdirs,ZDdirs)
		self.ALTgrid = 90.-self.ZDgrid
		
		def sky_brightness_point(az,zd):
			alt = 90.-zd
			fits_region_values = FitsImage.fits_data[\
				(np.array(ImageCoordinates.altitude_map>=alt-altseparation/2.)*\
				 np.array(ImageCoordinates.altitude_map<alt+altseparation/2.)\
				)*(\
				(np.array(ImageCoordinates.azimuth_map>=az-azseparation/2.)*\
				 np.array(ImageCoordinates.azimuth_map<az+azseparation/2.)\
				)+\
				(np.array(ImageCoordinates.azimuth_map>=az-azseparation/2.+360)+\
				 np.array(ImageCoordinates.azimuth_map<az+azseparation/2.-360)\
				))]
				# The last two lines must be set to get a smooth transition between 360 and 0 degrees azimuth
			
			return(self.sky_brightness_region(BouguerFit,ImageInfo,fits_region_values))
		
		sky_brightness_list = np.vectorize(sky_brightness_point)
		self.SBgrid,self.SBgrid_errors = np.array(sky_brightness_list(self.AZgrid,self.ZDgrid))
		
		# Once we measured the sky brightness in the image, convert to radians the azimuths
		self.AZgrid = self.AZgrid*np.pi/180.
		
		# Griddata.
		self.AZgridi,self.ZDgridi = np.mgrid[0:2*np.pi:1000j, 0:75:1000j]
		self.ALTgridi = 90. - self.ZDgridi
		
		coord_reshape = [[self.AZgrid[j][k],self.ZDgrid[j][k]] \
			for k in xrange(len(self.AZgrid[0])) for j in xrange(len(self.AZgrid))]
		
		data_reshape = [ self.SBgrid[j][k] \
			for k in xrange(len(self.AZgrid[0])) for j in xrange(len(self.AZgrid))]
		
		self.SBgridi = sint.griddata(coord_reshape,data_reshape, \
			(self.AZgridi,self.ZDgridi), method='linear')
	
	def measure_in_positions(self):
		# Measure Sky Brightness at zenith
		self.SBzenith = np.median(self.SBgrid[self.ZDgrid==0])
		self.SBzenith_err = np.max(self.SBgrid_errors[self.ZDgrid==0])
		

class SkyBrightnessGraph():
	def __init__(self,SkyBrightness,ImageInfo,BouguerFit):
		self.create_plot()
		self.plot_labels(SkyBrightness,ImageInfo,BouguerFit)
		self.define_contours(ImageInfo)
		self.ticks_and_locators()
		self.plot_data(SkyBrightness)
		self.color_bar()
		self.show_map()
	
	def create_plot(self):
		''' Create the figure (empty) with matplotlib '''
		self.SBfigure = plt.figure(figsize=(8,8))
		self.SBgraph  = self.SBfigure.add_subplot(111,projection='polar')
	
	def plot_data(self,SkyBrightness):
		''' Returns the graph with data plotted.'''
		self.SBcontoursf = self.SBgraph.contourf(\
			SkyBrightness.AZgridi, SkyBrightness.ZDgridi, SkyBrightness.SBgridi, cmap=plt.cm.YlGnBu,levels=self.level_list)
		self.SBcontours  = self.SBgraph.contour(\
			SkyBrightness.AZgridi, SkyBrightness.ZDgridi, SkyBrightness.SBgridi,
			colors='k',alpha=0.3,levels=self.coarse_level_list)
		self.SBcontlabel = self.SBgraph.clabel(self.SBcontours,inline=True,fmt='%.1f',fontsize=10)
		# Limit radius
		self.SBgraph.set_ylim(0,75)
	
	def plot_labels(self,SkyBrightness,ImageInfo,BouguerFit):
		''' Set the figure title and add extra information (annotation) '''
		# Image title
		self.SBgraph.text(0,90, unicode(ImageInfo.backgroundmap_title,'utf-8'),\
			horizontalalignment='center',size='xx-large')
		
		# Image information
		image_information = str(ImageInfo.date_string)+" UTC\n"+str(ImageInfo.latitude)+5*" "+\
			str(ImageInfo.longitude)+"\n"+ImageInfo.used_filter+4*" "+\
			"K="+str("%.3f" % float(BouguerFit.Regression.extinction))+"+-"+\
			str("%.3f" % float(BouguerFit.Regression.error_extinction))+"\n"+\
			"SB="+str("%.2f" % float(SkyBrightness.SBzenith))+"+-"+\
			str("%.2f" % float(SkyBrightness.SBzenith_err))+" mag/arcsec2 (zenith)"
	
		self.SBgraph.text(5*np.pi/4,125,unicode(image_information,'utf-8'),fontsize='x-small')
	
	def define_contours(self,ImageInfo):
		''' Calculate optimal contours for pyplot.contour and pyplot.contourf '''
		
		_min_ = float(ImageInfo.background_levels[ImageInfo.used_filter][0])
		_max_ = float(ImageInfo.background_levels[ImageInfo.used_filter][1])
		
		sval = 0.1
		def create_ticks(_min_,_max_,sval):
			return  np.arange(_min_,_max_+sval/10.,sval)
		
		self.level_list = create_ticks(_min_,_max_,0.1)
		self.label_list = create_ticks(_min_,_max_,sval)
		self.coarse_level_list = create_ticks(_min_,_max_,0.2)
		
		while len(self.label_list)>15:
			sval = sval*2.
			self.label_list = create_ticks(_min_,_max_,sval)
		
		if len(self.level_list)<3: self.update_ticks=False
		else: self.update_ticks=True
	
	def ticks_and_locators(self):
		''' Add ticks to the graph '''
		radial_locator = np.arange(10,90+1,10)
		radial_label = ["$80$","$70$","$60$","$50$","$40$","$30$","$20$","$10$","$0$"]
		theta_locator = np.arange(0,360,45)
		theta_label = ["$N$","$NE$","$E$","$SE$","$S$","$SW$","$W$","$NW$"]
		self.SBgraph.set_rgrids(radial_locator,radial_label,size="large",color='k',alpha=0.75)
		self.SBgraph.set_thetagrids(theta_locator,theta_label,size="large")
		# rotate the graph (North up)
		self.SBgraph.set_theta_direction(-1)
		self.SBgraph.set_theta_offset(np.pi/2)
	
	def color_bar(self):
		''' Add the colorbar '''
		# Separation between colour bar and graph
		self.SBfigure.subplots_adjust(right=1)
		# Color bar 
		self.SBcolorbar = plt.colorbar(self.SBcontoursf,orientation='vertical',shrink=0.85)
		self.SBfigure.subplots_adjust(right=0.80) # Restore separation
		self.SBcolorbar.set_ticks(self.label_list,update_ticks=self.update_ticks)
		self.SBcolorbar.set_label("mag/arcsec2",rotation="vertical",size="large")
	
	def show_map(self):
		#plt.show(self.SBfigure)
		self.SBfigure.savefig("/home/minaya/fondo_cielo_pyasb.png")
		
