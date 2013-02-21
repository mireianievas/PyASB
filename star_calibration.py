#!/usr/bin/env python

'''
Load Catalog file and make PyASB StarCatalog

This module loads the catalog file and returns 
an array of Star objects (StarCatalog) with
their fluxes.

____________________________

This module is part of the PyASB project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

DEBUG = False

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
	import ephem
	import math
	from astrometry import *
	from skymap_plot import *
except:
	print(str(sys.argv[0]) + ': One or more modules missing: pyephem,math,astrometry')
	raise SystemExit

class Star():
	def __init__(self,StarCatalogLine,FitsImage,ImageInfo,ObsPyephem):
		''' Takes StarCatalogLine (line from catalog file) and 
		      FitsImage, ImageInfo and ObsPyephem objects
		    Returns a Star object with photometric and astrometic properties 
		      or a destroy flag if errors ocurred during process'''
		self.destroy=False
		if self.destroy==False: 
			#if DEBUG==True:
			#	print('Extracting info from catalog ...')
			self.from_catalog(StarCatalogLine)
			if self.destroy==True: 
				if DEBUG==True:
					print('Error extracting from catalog')
		
		if self.destroy==False:
			#if DEBUG==True:
			#	print('Reading used filter ...')
			self.magnitude_on_image(ImageInfo)
			if self.destroy==True: 
				if DEBUG==True:
					print(self.name + ' Error reading used filter or magnitude>max')
			
		if self.destroy==False:
			#if DEBUG==True:
			#	print('Star astrometry ...')
			self.star_astrometry(ObsPyephem,ImageInfo)
			if self.destroy==True: 
				if DEBUG==True:
					print(self.name + ' Error performing star astrometry, maybe star isnt visible?')
		
		if self.destroy==False:
			#if DEBUG==True:
			#	print('Photometric radius ...') 
			self.photometric_radius(ImageInfo)
			if self.destroy==True: 
				if DEBUG==True:
					print(self.name + ' Error generating photometric radius')
		
		if self.destroy==False:
			#if DEBUG==True:
			#	print('Estimate fits region ...')
			self.estimate_fits_region_star(FitsImage)
			self.estimate_fits_region_complete(FitsImage)
			if self.destroy==True: 
				if DEBUG==True:
					print(self.name + ' Error creating regions of stars and star+background')
		
		if self.destroy==False:
			#if DEBUG==True:
			#	print('Measure star flux ...')
			self.measure_star_fluxes(FitsImage.fits_data)
			if self.destroy==True: 
				if DEBUG==True:
					print(self.name + ' Error measuring fluxes')
			
		if self.destroy==False:
			#if DEBUG==True:
			#	print('Test if star is saturated ...')
			self.star_is_saturated(ImageInfo)
			if self.destroy==True: 
				if DEBUG==True:
					print(self.name + ' Error, star is saturated')
		
		if self.destroy==False:
			#if DEBUG==True:
			#	print('Test if star is detectable ...')
			self.star_is_detectable(ImageInfo)
			if self.destroy==True: 
				if DEBUG==True:
					print(self.name + ' Error, star isnt detectable')
		
		if self.destroy==False:
			#if DEBUG==True:
			#	print('Estimate centroid ...')
			self.estimate_centroid()
			if self.destroy==True: 
				if DEBUG==True:
					print(self.name + ' Error estimating centroid')
		
		if self.destroy==False:
			#if DEBUG==True:
			#	print('Optimize aperture photometry ...')
			self.optimal_aperture_photometry(ImageInfo,FitsImage.fits_data)
			if self.destroy==True: 
				if DEBUG==True:
					print(self.name + ' Error doing optimal photometry')
		
		if self.destroy==False:
			#if DEBUG==True:
			#	print('Photometric bouguer variables ...')
			self.photometry_bouguervar(ImageInfo)
			if self.destroy==True: 
				if DEBUG==True:
					print(self.name + ' Error calculating bouguer variables')
		
		if DEBUG==True:
			print('Clearing Object')
		self.__clear__()
		
		if self.destroy==False:
			if DEBUG==True:
				print('DONE')
		
	def from_catalog(self,CatalogLine):
		''' Populate class with properties extracted from catalog:
		    recno, HDcode, RA1950, DEC1950, Vmag, U_V, B_V, R_V, I_V '''
		
		def coord_pyephem_format(coord_str):
			# Return coordinate in pyepheem str
			coord_separated = coord_str.split(' ')
			coord_pyephem = str(int(coord_separated[0]))+\
				':'+str(int(coord_separated[1]))+":"+str(float(coord_separated[2]))
			return coord_pyephem
		
		try:
			self.recno   = int(CatalogLine[0])
			self.HDcode  = str(CatalogLine[1]).replace(' ','')
			self.RA1950  = coord_pyephem_format(CatalogLine[2])
			self.DEC1950 = coord_pyephem_format(CatalogLine[3])
			self.Vmag    = float(CatalogLine[4])
			self.U_V     = float(CatalogLine[5])
			self.B_V     = float(CatalogLine[6])
			self.R_V     = float(CatalogLine[7])
			self.I_V     = float(CatalogLine[8])
			try:
				''' Try to find the common name '''
				self.name = str(CatalogLine[9])
			except:
				''' Use the HDcode as name '''
				self.name = self.HDcode
			self.Umag = self.Vmag + self.U_V
			self.Bmag = self.Vmag + self.B_V 
			self.Rmag = self.Vmag + self.R_V
			self.Imag = self.Vmag + self.I_V
		except:
			self.destroy=True
	
	def magnitude_on_image(self,ImageInfo):
		''' Set the magnitude and color (#-V) that match image filter.'''
		if ImageInfo.used_filter=="Johnson_U":
			self.FilterMag = self.Umag
			self.Color     = self.U_V
		elif ImageInfo.used_filter=="Johnson_B":
			self.FilterMag = self.Umag
			self.Color     = self.U_V
		elif ImageInfo.used_filter=="Johnson_V":
			self.FilterMag = self.Vmag
			self.Color     = 0.0
		elif ImageInfo.used_filter=="Johnson_R":
			self.FilterMag = self.Rmag
			self.Color     = self.R_V
		elif ImageInfo.used_filter=="Johnson_I":
			self.FilterMag = self.Imag
			self.Color     = self.I_V
		else: self.destroy=True
		
		try:
			assert(self.FilterMag<ImageInfo.max_magnitude)
		except:
			self.destroy=True
	
	def star_astrometry(self,ObsPyephem,ImageInfo):
		''' Perform astrometry. Returns (if star is visible and well defined) its position on the sky and image'''
		
		def pyephem_declaration(self,ObsPyephem):
			''' Define the star in Pyephem to make astrometric calculations '''
			pyephem_star = ephem.FixedBody()
			pyephem_star = ephem.readdb('"'+str(self.name)+'"'+",f|S|A0,"+str(self.RA1950)+'|0'+\
				","+str(self.DEC1950)+'|0'+","+str(self.Vmag)+',1950,0"')
			pyephem_star.compute(ObsPyephem)
			return pyephem_star
		
		try:
			pyephem_star = pyephem_declaration(self,ObsPyephem)
		except:
			self.destroy=True
			
		if self.destroy==False:
			# Real coordinates (from catalog)
			self.altit_real = float(pyephem_star.alt)*180.0/math.pi
			self.zdist_real = 90.0-self.altit_real
			self.azimuth    = float(pyephem_star.az)*180.0/math.pi
			try:
				assert(self.altit_real)>float(ImageInfo.min_altitude)
			except:
				self.destroy=True
		
		if self.destroy==False:
			# Apparent coordinates in sky. Atmospheric refraction effect.
			self.altit_appa = atmospheric_refraction(self.altit_real,'dir')
			self.zdist_appa = 90.0-self.altit_appa
			self.airmass    = calculate_airmass(self.altit_appa)
			try:
				assert(self.altit_appa)>float(ImageInfo.min_altitude)
			except:
				self.destroy=True
		
		if self.destroy==False:
			# Image coordinates
			XYCoordinates = horiz2xy(self.azimuth,self.altit_appa,ImageInfo)
			self.Xcoord = XYCoordinates[0]
			self.Ycoord = XYCoordinates[1]
			try:
				assert(self.Xcoord>0. and self.Xcoord<ImageInfo.resolution[0])
				assert(self.Ycoord>0. and self.Ycoord<ImageInfo.resolution[1])
			except:
				self.destroy=True
	
	def photometric_radius(self,ImageInfo):
		''' Needs astrometry properties, photometric filter properties and ImageInfo
		    Returns R1,R2 and R3 '''
		try:
			'''Returns R1,R2,R3. Needs photometric properties and astrometry.'''
			MF_magn = 10**(-0.4*self.FilterMag)
			MF_reso = 0.5*(min(ImageInfo.resolution)/2500)
			MF_texp = 0.1*ImageInfo.exposure
			MF_airm = 0.7*self.airmass
			MF_totl = 1+MF_magn+MF_reso+MF_texp+MF_airm
			
			self.R1 = int(ImageInfo.base_radius*MF_totl)
			self.R2 = self.R1*2
			self.R3 = self.R1*3
		except:
			self.destroy=True
	
	def estimate_fits_region_star(self,FitsImage):
		''' Return the region that contains the star (both for calibrated and uncalibrated data)'''
		self.fits_region_star = [[FitsImage.fits_data[y,x] \
			for x in xrange(int(self.Xcoord - self.R1 + 0.5),int(self.Xcoord + self.R1 + 0.5))] \
			for y in xrange(int(self.Ycoord - self.R1 + 0.5),int(self.Ycoord + self.R1 + 0.5))]
		# We will need this to look for saturated pixels.
		self.fits_region_star_uncalibrated = [[FitsImage.fits_data[y,x] \
			for x in xrange(int(self.Xcoord - self.R1 + 0.5),int(self.Xcoord + self.R1 + 0.5))] \
			for y in xrange(int(self.Ycoord - self.R1 + 0.5),int(self.Ycoord + self.R1 + 0.5))]
	
	def estimate_fits_region_complete(self,FitsImage):
		''' Return the region that contains the star+background '''
		self.fits_region_complete = [[FitsImage.fits_data[y,x] \
			for x in xrange(int(self.Xcoord - self.R3 + 0.5),int(self.Xcoord + self.R3 + 0.5))] \
			for y in xrange(int(self.Ycoord - self.R3 + 0.5),int(self.Ycoord + self.R3 + 0.5))]
	
	def measure_star_fluxes(self,fits_data):
		'''Needs self.Xcoord, self.Ycoord and self.R[1-3] defined
		   Returns star fluxes'''
		
		# Pixels in each ring
		def less_distance(Xi,Yi,reference):
			# returns True if distance from pixel to the star center is less than a value.
			# False otherwise
			return (Xi)**2 + (Yi)**2 <= reference**2
		
		try:
			self.pixels1 = [self.fits_region_complete[y][x] \
				for y in xrange(len(self.fits_region_complete))\
				for x in xrange(len(self.fits_region_complete[0])) \
				if less_distance(x-len(self.fits_region_complete)/2.,y-len(self.fits_region_complete[0])/2.,self.R1)]
			
			self.pixels2 = [self.fits_region_complete[y][x] \
				for y in xrange(len(self.fits_region_complete))\
				for x in xrange(len(self.fits_region_complete[0])) \
				if less_distance(x-len(self.fits_region_complete)/2.,y-len(self.fits_region_complete[0])/2.,self.R2) and\
				not less_distance(x-len(self.fits_region_complete)/2.,y-len(self.fits_region_complete[0])/2.,self.R1)]
			
			self.pixels3 = [self.fits_region_complete[y][x] \
				for y in xrange(len(self.fits_region_complete))\
				for x in xrange(len(self.fits_region_complete[0])) \
				if less_distance(x-len(self.fits_region_complete)/2.,y-len(self.fits_region_complete[0])/2.,self.R3) and\
				not less_distance(x-len(self.fits_region_complete)/2.,y-len(self.fits_region_complete[0])/2.,self.R2)]
			
			# Sky background flux
			self.skyflux = 2.5*np.median(self.pixels3)-1.5*np.mean(self.pixels3)
			self.skyflux_err = np.std(self.pixels3)
			# Sky background + Star flux
			totalflux = np.sum(self.pixels1)
			# Only star flux
			self.starflux = totalflux - len(self.pixels1)*self.skyflux
			self.starflux_err = math.sqrt(len(self.pixels1))*self.skyflux_err
			
		except:
			self.destroy=True
	
	def star_is_saturated(self,ImageInfo):
		''' Return true if star has one or more saturated pixels 
		    requires that self.fits_region_star is defined'''
		try:
			assert(np.max(self.fits_region_star_uncalibrated)<(2./3)*2**ImageInfo.ccd_bits)
		except:
			self.destroy=True
	
	def star_is_detectable(self,ImageInfo):
		''' Set a detection limit to remove weak stars'''
		''' Check if star is detectable '''
		try:
			assert(self.starflux>ImageInfo.baseflux_detectable*self.skyflux_err*\
				math.sqrt(ImageInfo.exposure)*1./self.airmass)
		except:
			self.destroy=True
	
	def estimate_centroid(self,iterations=1):
		''' Returns star centroid from a region that contains the star
		    needs self.R1'''
		
		try:
			exponent = 2 # Values > 1 intensify the convergence of the method.
			xweight = np.array([range(1,len(self.fits_region_star[0])+1)]*len(self.fits_region_star))
			yweight = np.array([range(1,len(self.fits_region_star)+1)]*len(self.fits_region_star[0])).transpose()
			xcentroid = np.sum(xweight*np.power(self.fits_region_star-self.skyflux,exponent))\
				/np.sum(np.power(self.fits_region_star-self.skyflux,exponent))
			ycentroid = np.sum(yweight*np.power(self.fits_region_star-self.skyflux,exponent))\
				/np.sum(np.power(self.fits_region_star-self.skyflux,exponent))
			self.Xcoord += xcentroid - len(self.fits_region_star[0])/2
			self.Ycoord += ycentroid - len(self.fits_region_star)/2
		except:
			self.destroy=True
	
	def optimal_aperture_photometry(self,ImageInfo,fits_data):
		'''
		Optimize the aperture to minimize uncertainties and assert
		all flux is contained in R1
		'''
		try:
			radius = (ImageInfo.base_radius+self.R1)/2.
			iterate = True
			num_iterations = 0
			
			self.starflux = 0
			while iterate:
				num_iterations+=1
				old_starflux = self.starflux
				self.R1 = radius
				self.measure_star_fluxes(fits_data)
				if self.starflux < (1+0.01*num_iterations**2)*old_starflux:
					iterate=False
				else:
					radius+=1
				
				assert(radius<self.R2)
		except:
			self.destroy=True
	
	def photometry_bouguervar(self,ImageInfo):
		# Calculate parameters used in bouguer law fit
		try:
			_25logF      = 2.5*math.log10(self.starflux/ImageInfo.exposure)
			_25logF_unc  = (2.5/math.log(10))*self.starflux_err/self.starflux
			color_term = ImageInfo.color_terms[ImageInfo.used_filter][0]
			color_term_err = ImageInfo.color_terms[ImageInfo.used_filter][1]
			self.m25logF     = self.FilterMag+_25logF+color_term*self.Color
			self.m25logF_unc = math.sqrt(_25logF_unc**2 + (color_term_err*self.Color)**2)
		except:
			self.destroy=True
	
	def __clear__(self):
		backup_attributes = [\
			"destroy","name","FilterMag","Color",\
			"RA1950","DEC1950","azimuth","altit_real","airmass","Xcoord","Ycoord",\
			"R1","R2","R3","starflux","starflux_err","m25logF","m25logF_unc"]
		for atribute in list(self.__dict__):
			if atribute[0]!="_" and atribute not in backup_attributes:
				del vars(self)[atribute]
	

class StarCatalog():
	''' This class processes the catalog.
	    Takes FitsImage,ImageInfo,ObsPyephem, returns an object with 
	    the processed Star list'''
	
	def __init__(self,FitsImage,ImageInfo,ObsPyephem):
		print('Creating Star Catalog ...')
		self.load_catalog_file(ImageInfo.catalog_filename)
		print('Star processing ...')
		self.process_catalog(FitsImage,ImageInfo,ObsPyephem)
	
	def load_catalog_file(self,catalog_filename):
		''' Returns Catalog lines from the catalog_filename '''
		try:
			self.catalogfile = open(catalog_filename, 'r')
			CatalogContent = self.catalogfile.readlines()
			MinLine = 31; separator = "\t"
			self.CatalogLines = [ CatalogContent[line][:-1].split(separator) \
				for line in xrange(len(CatalogContent)) \
				if line>=MinLine-1 and CatalogContent[line].replace(" ","")!=""]
		except IOError:
			print('IOError. Error opening file '+catalog_filename+'.')
			#return 1
		except:
			print('Unknown error:')
			raise
			#return 2
		else:
			print('File '+str(catalog_filename)+' opened correctly.')
	
	def process_catalog(self,FitsImage,ImageInfo,ObsPyephem):
		'''Returns the processed catalog. '''
		self.StarList = []
		self.StarList_woPhot = []
		for each_star in self.CatalogLines:
			PhotometricStar = Star(each_star,FitsImage,ImageInfo,ObsPyephem)
			if PhotometricStar.destroy==False:
				self.StarList.append(PhotometricStar)
			try:
				PhotometricStar.Xcoord
				PhotometricStar.Ycoord
			except:
				pass
			else:
				self.StarList_woPhot.append(PhotometricStar)
		
		print(" - Total stars: "+str(len(self.StarList_woPhot)))
		print(" - With photometry: " +str(len(self.StarList)))
		

