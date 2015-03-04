
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


import numpy as np
import scipy.interpolate as sint
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
import matplotlib.patches as mpp


class SkyBrightness(object):
	'''
	Class with Sky Brightness measure methods.
	Init requires FitsImage, ImageCoordinates, ImageInfo and Bouguerfit objects.
	'''
	
	def __init__(self,FitsImage,ImageInfo,ImageCoordinates,BouguerFit):
		if ImageInfo.skybrightness_map_path!=False:
			print('Measuring All-Sky Sky Brightness ...')
			#NOTE: This function is very slow, I need to figure how to improve it.
			self.measure_in_grid(FitsImage,ImageInfo,ImageCoordinates,BouguerFit)
			self.grid_text_to_file(ImageInfo)
		else:
			print('Measuring SB only at zenith ...')
		self.measure_in_positions(FitsImage,ImageInfo,ImageCoordinates,BouguerFit)
	
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
		zdseparation = 15
		
		self.AZdirs  = np.arange(0,360+1,azseparation)
		self.ZDdirs = np.arange(0,90+1,zdseparation)
		self.ALTdirs = 90-self.ZDdirs
		
		self.AZgrid,self.ZDgrid = np.meshgrid(self.AZdirs,self.ZDdirs)
		self.ALTgrid = 90.-self.ZDgrid
		
		def sky_brightness_point(az,zd):
			alt = 90.-zd
			fits_region_values = FitsImage.fits_data[\
				(np.array(ImageCoordinates.altitude_map>=alt-zdseparation/2.)*\
				 np.array(ImageCoordinates.altitude_map<alt+zdseparation/2.)\
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
	
	def grid_text_to_file(self,ImageInfo):
		try:
			assert(ImageInfo.skybrightness_table_path!=False)
		except:
			print('Skipping write skybrightness table to file')
		else:
			print('Write skybrightness table to file')
			def sbtable_filename(ImageInfo):
				filename = ImageInfo.skybrightness_table_path +\
					"/SBTable_"+ImageInfo.obs_name+"_"+ImageInfo.fits_date+"_"+\
					ImageInfo.used_filter+".txt"
				return(filename)
			
			photfile = open(sbtable_filename(ImageInfo),'w+')
			
			header = '#Altitude\Azimuth'
			for az_ in  self.AZdirs:
				header += ', '+str(az_)
			header+='\n'
			
			content = [header]
			
			for k,alt_ in enumerate(self.ALTdirs):
				line = str(alt_)
				for j in xrange(len(self.AZdirs)):
					line += ', '+str("%.3f" % float(self.SBgrid[k][j])) + ' +/- ' + str("%.3f" % float(self.SBgrid_errors[k][j]))
				line+='\n'
				content.append(line)
			
			photfile.writelines(content)
			photfile.close()
		
	
	def measure_in_positions(self,FitsImage,ImageInfo,ImageCoordinates,BouguerFit):
		# Measure Sky Brightness at zenith
		try:
			# If previous grid calculus, then extract from that grid
			self.SBzenith = np.median(self.SBgrid[self.ZDgrid==0])
			self.SBzenith_err = np.max(self.SBgrid_errors[self.ZDgrid==0])
		except:
			# If not previous grid, calculate manually.
			zenith_acceptance = 10
			fits_zenith_region_values = FitsImage.fits_data[\
				ImageCoordinates.altitude_map>=90-zenith_acceptance]
			self.SBzenith,self.SBzenith_err = \
				self.sky_brightness_region(BouguerFit,ImageInfo,fits_zenith_region_values)


class SkyBrightnessGraph(object):
	def __init__(self,SkyBrightness,ImageInfo,BouguerFit):
		if ImageInfo.skybrightness_map_path==False:
			# Don't draw anything
			print('Skipping SkyBrightness Graph ...')
			return(None)
		else:
			print('Generating Sky Brightness Map ...')
		
		self.create_plot()
		self.plot_labels(SkyBrightness,ImageInfo,BouguerFit)
		self.define_contours(ImageInfo)
		self.ticks_and_locators()
		self.grid_data(SkyBrightness)
		self.plot_data()
		self.color_bar()
		self.show_map(ImageInfo)
		plt.clf()
		plt.close('all')
	
	def create_plot(self):
		''' Create the figure (empty) with matplotlib '''
		self.SBfigure = plt.figure(figsize=(8,8))
		self.SBgraph  = self.SBfigure.add_subplot(111,projection='polar')
	
	def grid_data(self,SkyBrightness):
		# Griddata.
		self.AZgridi,self.ZDgridi = np.mgrid[0:2*np.pi:1000j, 0:75:1000j]
		self.ALTgridi = 90. - self.ZDgridi
		
		coord_reshape = [[SkyBrightness.AZgrid[j][k],SkyBrightness.ZDgrid[j][k]] \
			for k in xrange(len(SkyBrightness.AZgrid[0])) for j in xrange(len(SkyBrightness.AZgrid))]
		
		data_reshape = [ SkyBrightness.SBgrid[j][k] \
			for k in xrange(len(SkyBrightness.AZgrid[0])) for j in xrange(len(SkyBrightness.AZgrid))]
		
		self.SBgridi = sint.griddata(coord_reshape,data_reshape, \
			(self.AZgridi,self.ZDgridi), method='linear')
	
	def plot_data(self):
		''' Returns the graph with data plotted.'''
		self.SBcontoursf = self.SBgraph.contourf(\
			self.AZgridi, self.ZDgridi, self.SBgridi, cmap=plt.cm.YlGnBu,levels=self.level_list)
		self.SBcontours  = self.SBgraph.contour(\
			self.AZgridi, self.ZDgridi, self.SBgridi,
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
	
	def show_map(self,ImageInfo):
		def skybrightness_filename(ImageInfo):
			filename = ImageInfo.skybrightness_map_path +\
				"/SkyBrightnessMap_"+ImageInfo.obs_name+"_"+ImageInfo.fits_date+"_"+\
				ImageInfo.used_filter+".png"
			return(filename)
		
		#plt.show(self.SBfigure)
		if ImageInfo.skybrightness_map_path=="screen":
			plt.show()
		else:
			plt.savefig(skybrightness_filename(ImageInfo),bbox_inches='tight')
		
		
