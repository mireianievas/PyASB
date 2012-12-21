#!/usr/bin/env python

'''
Star detection module

Perform star detection and measure star fluxes
____________________________

This module is part of the PyAstMonUCM project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyAstMonUCM project"
__credits__ = ["Miguel Nievas"]
__license__ = "GPL"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


try:
	import numpy as np
	import matplotlib.pyplot as mpl
	import matplotlib.colors as mpc
	import matplotlib.patches as mpp
	from skymap_plot import * 
except:
	print 'One or more modules missing: pyfits,HeaderTest'
	raise SystemExit

class MeasuredStar(MeasuredCatalog):
	''' Derivated variables used in Bouguer fit'''
	def photometry_bouguervar(self):
		# Calculate parameters used in bouguer law fit
		self.25logF      = 2.5*log10(self.flux_star)
		self.25logF_unc  = (2.5/log(10))*,self.flux_star_unc/self.flux_star
		self.m25logF     = self.FilterMag+Star.25logF
		self.m25logF_unc = self.25logF_unc
	
	''' Background and Star Fluxes '''
	def measure_background(self,fits_region):
		self.skyflux     = np.median(fits_region)
		self.skyflux_err = np.std(fits_region)
		
	def measure_totalflux(self,fits_region):
		self.totalflux    = np.sum(fits_region)
		
	def measure_starflux(self,pixels_star,pixels_background):
		measure_background(self,pixels_background)
		measure_totalflux(self,pixels_star)
		self.starflux = self.totalflux - len(pixel_star)*self.skyflux
		self.starflux_err = len(pixel_star)*self.skyflux_err
	
	''' Centroid determination '''
	def centroid_estimation(self,fits_region):
		exponent = 2 # Valures > 1 intensify the convergence of the method.
		detection_threshold = 1.5
		try:
			xweight = np.array([range(len(fits_region[0]))]*len(fits_region))
			yweight = np.array([range(len(fits_region))]*len(fits_region[0])).transpose()
			self.xcentroid = np.sum(xweight*fits_region**exponent)/np.sum(fits_region**exponent)
			self.ycentroid = np.sum(yweight*fits_region**exponent)/np.sum(fits_region**exponent)
			assert np.std(fits_region)>detection_threshold*np.mean(np.abs(fits_region))
		except:
			raise

	'''Apperture photometry'''		
		
	def star_detection(self,ImageInfo):
		# Define aperture disk for photometry
		R1,R2,R3 = aperture_photometry_radius(Star,ImageInfo)
		pixels1,pixels2,pixels3 = apphot_pixels(Star.Xcoord,Star.Ycoord,R1,R2,R3,ImageInfo)
		
		# Coarse flux measure
		aproximate_fluxes(self.fits_data,pixels1,pixels3,self.ImageInfo)
		# NOTE: flux_sky is absolute, not normalized.
		
		# Filter the fits matrix for future background SB measures
		for pixel in pixels1+pixels2:
			filtered_self.fits_data[pixel[0]][pixel[1]] = Star.flux_sky
		
		# Delete stars below detection limit or low at sky
		Star.max_magnitude = magnitude_limit(Star,self.ImageInfo)
		Star.min_flux      = flux_lowlimit(Star,self.ImageInfo)
		
		if Star.FilterMag < Star.max_magnitude and Star.flux_star > Star.min_flux\
		and Star.altit_appa > self.ImageInfo.Config.min_altitude:
			try:
				# Fine flux measure and star position (centroid estimation)
				Star.Xcoord,Star.Ycoord = estimate_centroid(self.fits_data,pixels1+pixels2,Star)
				pixels1,pixels2,pixels3 = \
					apphot_pixels(Star.Xcoord,Star.Ycoord,R1,R2,R3,self.ImageInfo)
				Star.flux_sky,Star.flux_sky_err   = sky_flux_measure(self.fits_data,pixels3,self.ImageInfo)
				Star.flux_star,Star.flux_star_err = \
					star_flux_measure(self.fits_data,pixels1,Star.flux_sky,Star.flux_sky_err,self.ImageInfo)
				# Log of flux and m+2.5logF, needed for bouguer law fit
				Star=photometry_bouguervar(Star)
				# Filter the fits matrix for future background SB measures
				for pixel in pixels1+pixels2:
					filtered_self.fits_data[pixel[0]][pixel[1]] = Star.flux_sky
				# Append to the measured star list and plot it if neccesary
				StarsMeasured.append(Star)

	
			except:
				''' Continue with next star, we cannot calculate centroid or measure 
				    fluxes for this one. '''
				print 'Cannot measure fluxes on '+str(Star.HDcode)+'. Trying next star...'
				continue
	

class MeasuredCatalog():
	def __init__(self.fits_data,self.StarCatalog,self.ImageInfo,self.ObsPyephem):
		self.fits_data     = fits_data
		self.StarsOriginal = StarCatalog.Stars
		self.ImageInfo     = ImageInfo
		self.ObsPyephem    = ObsPyephem
		self.filtered_fits_data = deepcopy(self.fits_data)
	

	
	




	def stellar_photometry(self):
		'''
		This function measures fluxes for Star in the self.StarCatalog
		and give the neccesary measures to fit a Bouguer extinction law.
		'''
	
		if self.ImageInfo.Config.skymap_file!=False:
			skyfigure,skyimage = create_skymap(self.fits_data,self.StarCatalog,self.ImageInfo,self.ObsPyephem)
	
		StarsMeasured = []
		for star_index in xrange(len(self.StarCatalog)):
			Star = self.StarCatalog[star_index]
			try:
				Star = _star_detection(Star)
			except:
				continue
				
				
				
				
				
			# Define aperture disk for photometry
			R1,R2,R3 = aperture_photometry_radius(Star,self.ImageInfo)
			pixels1,pixels2,pixels3 = apphot_pixels(Star.Xcoord,Star.Ycoord,R1,R2,R3,self.ImageInfo)
			if pixel1==[] or pixel1==[] or pixel3=[]:
				# Continue with next star, this one is not in the FoV.
				self.StarCatalog.pop(star_index)
				continue
			
			# Coarse flux measure
			Star.flux_star,Star.flux_sky = aproximate_fluxes(self.fits_data,pixels1,pixels3,self.ImageInfo)
			# NOTE: flux_sky is absolute, not normalized.
			
			# Filter the fits matrix for future background SB measures
			for pixel in pixels1+pixels2:
				filtered_self.fits_data[pixel[0]][pixel[1]] = Star.flux_sky
			
			# Delete stars below detection limit or low at sky
			Star.max_magnitude = magnitude_limit(Star,self.ImageInfo)
			Star.min_flux      = flux_lowlimit(Star,self.ImageInfo)
			
			if Star.FilterMag < Star.max_magnitude and Star.flux_star > Star.min_flux\
			and Star.altit_appa > self.ImageInfo.Config.min_altitude:
				try:
					# Fine flux measure and star position (centroid estimation)
					Star.Xcoord,Star.Ycoord = estimate_centroid(self.fits_data,pixels1+pixels2,Star)
					pixels1,pixels2,pixels3 = \
						apphot_pixels(Star.Xcoord,Star.Ycoord,R1,R2,R3,self.ImageInfo)
					Star.flux_sky,Star.flux_sky_err   = sky_flux_measure(self.fits_data,pixels3,self.ImageInfo)
					Star.flux_star,Star.flux_star_err = \
						star_flux_measure(self.fits_data,pixels1,Star.flux_sky,Star.flux_sky_err,self.ImageInfo)
					# Log of flux and m+2.5logF, needed for bouguer law fit
					Star=photometry_bouguervar(Star)
					# Filter the fits matrix for future background SB measures
					for pixel in pixels1+pixels2:
						filtered_self.fits_data[pixel[0]][pixel[1]] = Star.flux_sky
					# Append to the measured star list and plot it if neccesary
					StarsMeasured.append(Star)
					if self.ImageInfo.Config.skymap_file!=False:
						skyimage = annotate_skymap(skyimage,Star,R1,R2,R3)
	
				except:
					''' Continue with next star, we cannot calculate centroid or measure 
					    fluxes for this one. '''
					print 'Cannot measure fluxes on '+str(Star.HDcode)+'. Trying next star...'
					continue
	
		if self.ImageInfo.Config.skymap_file!=False:
			# Show or save the map
			show_or_save_skymap(skyfigure,self.ImageInfo,self.ObsPyephem)
		
		return StarsMeasured,filtered_self.fits_data

