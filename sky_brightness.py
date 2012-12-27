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
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


try:
	import matplotlib.pyplot as mpl
	import matplotlib.colors as mpc
	import matplotlib.patches as mpp
except:
	print('One or more modules missing: pyfits,HeaderTest')
	raise SystemExit



class SkyBrightness():
	'''
	Class with Sky Brightness measure methods.
	Init requires fits_data (numpy array with raw data), 
		ImageInfo and Regression.
	'''
	def __init__(self,fits_data,ImageInfo,Regression):
		print('SkyBrightness(): Loading calibration data ...'),
		try:
			self.fits_data = fits_data
			self.exposure    = ImageInfo.Properties.exposure
			self.pixel_scale = ImageInfo.Config.pixel_scale
			self.zp        = Regression.zp
			self.zp_err    = Regression.zp_err
		except:
			print('Failed')
			raise
		else:
			print('OK')
		
	''' Image processing '''
	def _sky_brightness_measure(self,image_region):
		'''
		Return measured sky brightness with its error at a given point in image.
		image_region must contain a pixel list
		'''

		# Measure Sky fluxes
		sky_flux,sky_flux_err = sky_flux_measure(self,image_region)
		# Compute Sky brightness in magnitudes
		sky_brightness,sky_brightness_err = \
			self.zp-2.5*log10(sky_flux/(self.exposure*self.pixel_scale)),\
			sqrt(self.zp_err**2 + (2.5*sky_flux_err/(log(10)*sky_flux))**2)
		
		return sky_brightness,sky_brightness_err
		
	def _sky_brightness_region(self,azimuth,altitude,radius):
		'''
		Return a class with SB measured at a given azimuth/altitude position
		with a fixed integration radius in pixels.
		'''
		Xcenter,Ycenter = horiz2xy(azimuth,altitude,ImageInfo)
		image_region = apphot_pixels(Xcenter,Ycenter,radius,radius,radius,ImageInfo)[0]
		sky_brightness,sky_brightness_err = sky_brightness_measure(self,image_region)
		
		return _sky_brightness_measure(self,image_region)
	
	def measure_in_grid(self):
		'''
		Measure sky brightness only at selected points in the sky
		Return the sky brightness table
		'''
		print('SkyBrightness(): Measuring Sky Brightness on the standard grid ...'),
		try:
			radius = 50 # px
			self.ALTgrid = np.arange(10,90+1,10) # degrees
			self.AZgrid  = np.arange(0,360+1,10) # degrees
			self.SBgrid  = [_sky_brightness_region(self,az,alt,radius) \
				for az in self.SBazlist for alt in self.SBaltlist]
		except:
			print('Failed')
			raise
		else:
			print('OK')
	
	def measure_positions(self,AZlist,ALTlist,RADlist=None):
		'''
		Measure sky brightness only at special selected regions (zenith p.e.)
		AZlist (deg), ALTlist (deg) and RAD (optional, in pixels) are lists
		'''
		print('SkyBrightness(): Measuring Sky Brightness on selected positions ...'),
		try:
			AZlist[len(ALTlist)-1];
			ALTlist[len(AZlist)-1];
			if RADlist==None:
				RADlist = [50]*len(AZlist)
			
			self.AZselected  = AZlist
			self.ALTselected = ALTlist
			self.SBselected  = [sky_brightness_region(self,AZlist[k],ALTlist[k],RADlist[k]) \
				for k in xrange(len(AZlist))]
		except:
			print('Failed')
			raise
		else:
			print('OK')

class SkyBrightnessGraph():
	'''
	Class with Sky Brightness graph methods
	'''
	def __init__(self,SkyBrightness,ImageInfo,ObsPyephem,Regression):
		print('SkyBrightnessGraph(): parameters and data loading ...'),
		try:
			# Measures
			self.AZgrid  = SkyBrightness.AZgrid  *pi/180.                   # degrees
			self.ALTgrid = SkyBrightness.ALTgrid *pi/180.                   # degrees
			self.Radial  = pi/180. - self.ALTgrid                            # degrees
			self.SBgrid  = [SBtuple[0] for SBtuple in SkyBrightness.SBgrid] # mag/arcsec2
			self.Flux    = 10**(-0.4*self.SBgrid)                           # Flux (not calibrated)
			# Sky Brightness at zenith
			self.SBzenith     = SkyBrightness.SBgrid[0]
			self.SBzenith_err = SkyBrightness.SBgrid[1]
			# Graph parameters
			self.Title     = ImageInfo.Config.SBTitle     # Graph title (str)
			self.ImageFilter   = ImageInfo.Properties.used_filter
			try: 
				self.ContourLimits = ImageInfo.Config.ContourLimits[self.ImageFilter]
			except:
				self.ContourLimits = [np.min(self.SBgrid),np.max(self.SBgrid)]
			# Other data	
			try:	
				extinction     = Regression.slope
				extinction_err = Regression.slope_err
				self.extinction_str = "K="+str("%.3f" % float(extinction))+"+-"+\
					str("%.3f" % float(extinction_err))+"\n"
			except:
				self.extinction_str = ""
			# Observatory
			self.ObsPyephem = ObsPyephem
		except:
			print('Failed')
			raise
		else:
			print('OK')
	
	def grid_data(self):
		'''
		Grid scattered data
		'''
		print('SkyBrightnessGraph(): gridding data ...'),
		try:
			self.Radiali  = np.linspace(0*pi/180,75*pi/180,76)
			self.AZgridi  = np.linspace(0,360*pi/180,361)
			self.Fluxi    = griddata((self.AZgrid,self.Radial),self.Flux,\
				(self.azimuthi[None,:],self.radiali[:,None]),method='linear',\
				fill_value=min(self.Flux))
			self.SBgridi  = -2.5*log10(self.Fluxi)
		except:
			print('Failed')
			raise
		else:
			print('OK')
	
	def define_contours(self):
		'''
		Optimize contours for a given data
		'''
		
		_min_ = float(self.ContourLimits[0])
		_max_ = float(self.ContourLimits[1])
		sval = 0.1
		
		def create_ticks(_min_,_max_,sval):
			return  np.arange(_min_,_max_+sval/10.,sval)
		
		self.level_list = create_ticks(_min_,_max_,0.1)
		self.label_list = create_ticks(_min_,_max_,sval)
		
		while len(self.label_list)>15:
			sval = sval*2.
			self.label_list = create_ticks(_min_,_max_,sval)
		
		if len(barra_marcas)<3: self.update_ticks=False
		else: self.update_ticks=True

		
	def create_plot(self,onscreen=False,writefile=False):
		'''
		Interpolate and plot the Sky Brightness map
		'''
		print('SkyBrightnessGraph(): ploting data ...'),
		try:
			SBfigure = plt.figure(figsize=(8,8),dpi=100)
			SBgraph  = SBfigure.add_subplot(111,projection='polar')
			SBgraph.text(0,pi/2, unicode(self.Title,'utf-8'),\
				horizontalalignment='center',size='xx-large')
				
			define_contours(self)
			
			# Contours
			SBcontours = SBgraph.contourf(self.AZgridi,self.Radiali,self.SBgridi,\
				cmap=cm.YlGnBu,levels=self.level_list)
		
			# Radial/azimuthal ticks and locators
			radial_locator = [ (num+1)*pi/18 for num in xrange(7) ]
			radial_label = ["$80$","$70$","$60$","$50$","$40$","$30$","$20$","$10$","$0$"]
			theta_locator = [ 45*num for num in xrange(8)]
			theta_label = ["$N$","$NE$","$E$","$SE$","$S$","$SW$","$W$","$NW$"]
			SBgraph.set_rgrids(radial_locator,radial_label,size="large")
			SBgraph.set_thetagrids(theta_locator,theta_label,size="large")
			SBgraph.set_theta_direction(-1)
			SBgraph.set_theta_offset(pi/2)
			SBgraph.grid(linewidth=0.8)
		
			# Separation between colour bar and graph
			subplots_adjust(right=1)
			
			# Color bar 
			SBcolorbar = colorbar(SBcontours,orientation='vertical',shrink=0.85)
			subplots_adjust(right=0.80)
			SBcolorbar.set_ticks(self.label_list,update_ticks=self.update_ticks)
			SBcolorbar.set_label("mag/arcsec2",rotation="vertical",size="large")

			# Image information
			image_information = str(self.ObsPyephem.date)+" UTC\n"+str(self.ObsPyephem.lat)+5*" "+\
				str(self.ObsPyephem.lon)+"\n"+self.ImageFilter+4*" "+self.extinction_str+\
				str("%.2f" % float(self.SBzenith))+"+-"+str("%.2f" % float(self.SBzenith_err))+\
				"mag/arcsec2 (zenith)"
	
			SBgraph.text(5*pi/4,125*pi/180,unicode(image_information,'utf-8'),fontsize='x-small')
	
			# Show or save the graph
			if onscreen == True:
				show(figurafondocielo)
			if writefile != False:
				savefig(output_file_name(writefile,self.ObsPyephem.date,self.ImageFilter))
			close(SBfigure); 
		except:
			print('Failed')
			raise
		else:
			print('OK')

		
