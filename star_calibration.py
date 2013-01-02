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

try:
	import ephem
	import math
	from astrometry import *
except:
	print 'One or more modules missing: pyfits,CommonErrors,HeaderTest'
	raise SystemExit

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



''' Basic Star Properties, from Photometric and Astrometric Catalog''' 
class CatalogStar():
	''' Extract information for each Star in the Catalog '''
	def __init__(self,CatalogLine,ObsPyephem,ImageInfo):
		self.basic_properties(CatalogLine)
		self.pyephem_declaration(ObsPyephem)
		self.set_actual_filter(ImageInfo)
		self.imagedep_properties(ObsPyephem,ImageInfo)
		self.destroy = False
	
	def coord_pyephem_format(self,coord_str):
		# Return coordinate in pyepheem str
		coord_separated = coord_str.split(' ')
		coord_pyephem = str(int(coord_separated[0]))+\
			':'+str(int(coord_separated[1]))+str(float(coord_separated[2]))
		return coord_pyephem
			
	def basic_properties(self,CatalogLine):
		''' Object Star. Original catalog values :
			recno, HDcode, RA1950, DEC1950, Vmag, U_V, B_V, R_V, I_V '''		
		self.recno   = int(CatalogLine[0])
		self.HDcode  = str(CatalogLine[1]).replace(' ','')
		self.RA1950  = self.coord_pyephem_format(CatalogLine[2])
		self.DEC1950 = self.coord_pyephem_format(CatalogLine[3])
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
		
	def pyephem_declaration(self,ObsPyephem):
		''' Define the star in Pyephem to make astrometric calculations '''
		pyephem_star = ephem.readdb('"'+str(self.name)+'"'+",f|S|A0,"+str(self.RA1950)+'|0'+\
			","+str(self.DEC1950)+'|0'+","+str(self.Vmag)+',1950,0"')
		pyephem_star.compute(ObsPyephem)	
		return pyephem_star
	
	def set_actual_filter(self,ImageInfo):
		''' Test which filter is in use '''
		used_filter = ImageInfo.used_filter
		if used_filter=="Johnson_U":
			self.FilterMag = self.Umag
			self.Color     = self.U_V
		elif used_filter=="Johnson_B":
			self.FilterMag = self.Umag
			self.Color     = self.U_V
		elif used_filter=="Johnson_V":
			self.FilterMag = self.Vmag
			self.Color     = 0.0
		elif used_filter=="Johnson_R":
			self.FilterMag = self.Rmag
			self.Color     = self.R_V
		elif used_filter=="Johnson_I":
			self.FilterMag = self.Imag
			self.Color     = self.I_V
		else: pass

	def imagedep_properties(self,ObsPyephem,ImageInfo):
		'''
		Calculate star position and dynamic properties.
		Don't proceed if the star altitude < 0.
		Return updated StarCatalog (image dependent)
		'''
		self.destroy = False
		try:
			pyephem_star = self.pyephem_declaration(ObsPyephem)
			if float(pyephem_star.alt)*180.0/math.pi > float(ImageInfo.min_altitude):
				print '1'
				# Real coordinates (from catalog)
				self.altit_real = float(pyephem_star.alt)*180.0/math.pi
				self.zdist_real = 90.0-self.altit_real
				self.azimuth    = float(pyephem_star.az)*180.0/math.pi
				print '2'
				# Apparent coordinates in sky. Atmospheric refraction effect.
				self.altit_appa = atmospheric_refraction(self.altit_real,'dir')
				self.zdist_appa = 90.0-self.altit_appa
				self.airmass    = calculate_airmass(self.altit_appa)
				print '3'
				# Photometric properties
				self.set_actual_filter(ImageInfo)
				print '4'
				# Image coordinates
				XYCoordinates = horiz2xy(self.azimuth,self.altit_appa,ImageInfo)
				self.Xcoord = XYCoordinates[0]
				self.Ycoord = XYCoordinates[1]
				print '5'
				if self.Xcoord<0. or self.Ycoord<0.or self.Xcoord>ImageInfo.resolution[0] \
				or self.Ycoord>ImageInfo.resolution[1] or self.altit_real < ImageInfo.min_altitude:
					# Star doesn't fit in the image
					self.destroy = True
			else:
				self.destroy = True
			
			if self.FilterMag > ImageInfo.max_magnitude:
				self.destroy = True
		except:
			raise
			self.destroy = True
			print 'destroyed'
		else:
			if self.destroy == False:
				print 'ok'
		
	def __del__(self):
		#print('Deleting star ...'),
		del(self)
		#print('OK')

''' Extended Star Properties, from image analysis''' 
class PhotometricStar(CatalogStar):
	''' Derivated variables to be used in Bouguer fit'''
	def photometry_bouguervar(self):
		# Calculate parameters used in bouguer law fit
		self._25logF      = 2.5*log10(self.flux_star)
		self._25logF_unc  = (2.5/log(10))*self.flux_star_unc/self.flux_star
		self.m25logF     = self.FilterMag+self._25logF
		self.m25logF_unc = self._25logF_unc
	
	''' Aperture photometry. Define R1, R2 and R3 radius '''
	def estimate_photometric_radius(self,ImageInfo):
		MF_magn = 10**(-0.4*self.FilterMag)
		MF_reso = 0.5*(min(ImageInfo.resolution)/2500)
		MF_texp = 0.5*ImageInfo.exposure
		MF_airm = 0.7*self.airmass
		MF_totl = 1+MF_magn+MF_reso+MF_texp+MF_airm
		
		self.R1 = int(ImageInfo.base_radius*MF_totl)
		self.R2 = self.R1*3
		self.R3 = self.R1*4
	
	''' Aperture photometry. List pixels in each ring '''
	def apphot_pixels(self,X,Y):
		X = int(X+0.5); Y = int(Y+0.5)
		pixels_region = [[X,Y] for X in xrange(X-self.R3-1,X+self.R3+1)\
			for Y in xrange(Y-self.R3-1,X+self.R3+1)]
		# Pixels in each ring
		def less_distance(Xi,Yi,reference):
			return (Xi-X)**2 + (Yi-Y)**2 <=reference**2
			
		self.pixels1 = [Pixel for Pixel in Pixels_Region if \
			less_distance(Pixel[0],Pixel[1],R1)]
		self.pixels2 = [Pixel for Pixel in Pixels_Region if \
			less_distance(Pixel[0],Pixel[1],R2) and\
			not less_distance(Pixel[0],Pixel[1],R1)]
		self.pixels3 = [Pixel for Pixel in Pixels_Region if \
			less_distance(Pixel[0],Pixel[1],R3) and\
			not less_distance(Pixel[0],Pixel[1],R2)]
	
	''' Turn pixel coordinates into pixel values '''
	@staticmethod
	def pixel_list_values(pixel_list,fits_data):
		return np.array([fits_data[pixel[0],pixel[1]] for pixel in pixel_list])
	
	''' Background and Star Fluxes '''
	def measure_background(self,fits_region):
		self.skyflux     = np.median(fits_region)
		self.skyflux_err = np.std(fits_region)
		
	def measure_totalflux(self,fits_region):
		self.totalflux = np.sum(fits_region)
		
	def measure_starflux(self,pixelvalues_star,pixelvalues_background):
		measure_background(self,pixelvalues_background)
		measure_totalflux(self,pixelvalues_star)
		self.starflux = self.totalflux - len(pixelvalues_star)*self.skyflux
		self.starflux_err = len(pixelvalues_star)*self.skyflux_err
	
	@staticmethod
	def star_is_saturated(fits_region,max_value):
		assert np.max(fits_region)<max_value
	
	
	''' Detectability conditions '''
	def star_min_flux_detectable(self,ImageInfo):
		#NOTE: Changed skyflux -> skyflux_err
		return ImageInfo.baseflux_detectable*self.skyflux_err*\
			math.sqrt(ImageInfo.exposure)*1./self.airmass	
	
	''' Remove the stars in the imagen. NOTE: Actually unnecesary '''
	@staticmethod
	def star_mask_image(fits_data,pixels,value):
		for pixel in pixels:
			fits_data[pixel[0]][pixel[1]]=value
	
	''' 
	Optimize photometry for good stars
	
	Centroid determination
	'''
	def estimate_centroid(self,fits_region):
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
	
	def optimal_aperture_photometry(self,min_radius,max_radius):
		'''
		Optimize the aperture to minimize uncertainties and assert
		all flux is contained in R1
		'''
		radius = min_radius
		continue_iterate = True
		while continue_iterate==True and radius<max_radius:
			radius += 1
			pixels1,pixels2,pixels3 = self.apphot_pixels(self.Xcoord,self.Ycoord)
			old_starflux = self.starflux
			region_background = self.pixel_list_values(pixels3,FitsImage.fits_data)
			region_star = self.pixel_list_values(pixels1,FitsImage.fits_data)
			measure_starflux(self,region_star,region_background)
			if self.starflux < old_starflux*1.02:
				continue_iterate = False
		assert radius<max_radius
	
	'''Aperture photometry. The main function'''	
	def stellar_aperture_photometry(self,FitsImage,ImageInfo):
		# Define aperture disk for photometry
		R1,R2,R3 = self.estimate_photometric_radius(ImageInfo)
		pixels1,pixels2,pixels3 = self.apphot_pixels(self.Xcoord,self.Ycoord)
		region_background = self.pixel_list_values(pixels3,FitsImage.fits_data)
		region_star = self.pixel_list_values(pixels1,FitsImage.fits_data)
		# Measure fluxes
		self.star_is_saturated(region_star)
		self.measure_starflux(region_star,region_background)
		
		# Test if its detectable
		assert self.starflux >= self.star_min_flux_detectable(ImageInfo)
		
		# More accurate photometry
		optimal_aperture_photometry(pixels1,R1/2,R2)
		photometry_bouguervar()


class CatalogFile():
	''' Catalog with Star photometry in U,B,V,R,I bands '''
	def __init__(self,catalog_filename):
		# Open catalog file and return its content
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
	
	def __del__(self):
		print('Closing file ...'),
		self.catalogfile.close()
		print('OK')

class PhotometricCatalog():
	''' This class processes the entire catalog '''
	def __init__(self,ObsPyephem,ImageInfo,catalog_filename='ducati_catalog.tsv'):
		# Load the catalog
		RawCatalog = CatalogFile(catalog_filename)
		
		# Calculate some properties for each star
		self.ProcessedCatalog = []
		for Star in RawCatalog.CatalogLines:
			PhotometricStar_ = PhotometricStar(Star,ObsPyephem,ImageInfo)
			self.ProcessedCatalog.append(PhotometricStar_)
		
		# Filter the catalog
		self.Stars = [Star for Star in self.ProcessedCatalog if Star.destroy==False]
		RawCatalog.__del__()
	
	''' Detect and measure star fluxes in the image. Destroy bad stars'''
	def photometric_measures(self,FitsImage,ImageInfo):
		try:
			for Star in self.ProcessedCatalog:
				Star.stellar_aperture_photometry(FitsImage,ImageInfo)
		except:
			raise
			print 'phot error!'
			Star.destroy=True
		self.Stars = [Star for Star in self.ProcessedCatalog if Star.destroy==False]
	
	def plot_starmap(self,ImageInfo,ObsPyephem):
		TheSkyMap = SkyMap(fits_data,MeasuredCatalog,ImageInfo,ObsPyephem)
		for Star in self.Stars:
			TheSkyMap.annotate_skymap(Star,R1,R2,R3)
		show_or_save_skymap()
	
	def __del__(self):
		print('Deleting catalog ...'),
		del(self)
		print('OK')
