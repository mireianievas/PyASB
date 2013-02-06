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

try:
	import ephem
	import math
	from astrometry import *
except:
	print 'One or more modules missing: pyephem,math,astrometry'
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
		self.destroy = False
		try:
			self.basic_properties(CatalogLine)
			self.pyephem_declaration(ObsPyephem)
			self.set_actual_filter(ImageInfo)
			self.imagedep_properties(ObsPyephem,ImageInfo)
		except:
			if DEBUG==True:
				print('Error reading Star');
				print(CatalogLine)
			
			self.destroy = True	
	
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
		self.destroy = True
		try:
			pyephem_star = self.pyephem_declaration(ObsPyephem)
			assert(float(pyephem_star.alt)*180.0/math.pi > float(ImageInfo.min_altitude))
			#print '1'
			# Real coordinates (from catalog)
			self.altit_real = float(pyephem_star.alt)*180.0/math.pi
			self.zdist_real = 90.0-self.altit_real
			self.azimuth    = float(pyephem_star.az)*180.0/math.pi
			#print '2'
			# Apparent coordinates in sky. Atmospheric refraction effect.
			self.altit_appa = atmospheric_refraction(self.altit_real,'dir')
			self.zdist_appa = 90.0-self.altit_appa
			self.airmass    = calculate_airmass(self.altit_appa)
			#print '3'
			# Photometric properties
			self.set_actual_filter(ImageInfo)
			#print '4'
			# Image coordinates
			XYCoordinates = horiz2xy(self.azimuth,self.altit_appa,ImageInfo)
			self.Xcoord = XYCoordinates[0]
			self.Ycoord = XYCoordinates[1]
			#print '5'
			assert(self.Xcoord>0. and self.Ycoord>0. and self.Xcoord<ImageInfo.resolution[0] \
			and self.Ycoord<ImageInfo.resolution[1] and self.altit_real > ImageInfo.min_altitude)
			# If we get this line, the start should survive.
			self.destroy=False
			#print('Star OK')
		except:
			pass
			#raise
			#print('destroyed')

	def __del__(self):
		#print('Deleting star ...'),
		del(self)
		#print('OK')

''' Extended Star Properties, from image analysis''' 
class PhotometricStar(CatalogStar):
	''' Derivated variables to be used in Bouguer fit'''
	def photometry_bouguervar(self):
		# Calculate parameters used in bouguer law fit
		self._25logF      = 2.5*math.log10(self.starflux)
		self._25logF_unc  = (2.5/math.log(10))*self.starflux_err/self.starflux
		self.m25logF     = self.FilterMag+self._25logF
		self.m25logF_unc = self._25logF_unc
	
	''' Aperture photometry. Define R1, R2 and R3 radius '''
	def estimate_photometric_radius(self,ImageInfo):
		MF_magn = 10**(-0.4*self.FilterMag)
		MF_reso = 0.5*(min(ImageInfo.resolution)/2500)
		MF_texp = 0.1*ImageInfo.exposure
		MF_airm = 0.7*self.airmass
		MF_totl = 1+MF_magn+MF_reso+MF_texp+MF_airm

		self.R1 = int(ImageInfo.base_radius*MF_totl)
		self.R2 = self.R1*3
		self.R3 = self.R1*4
		
		if DEBUG==True:
			print('MF_magn,MF_reso,MF_texp,MF_airm,MF_totl:')
			print(MF_magn,MF_reso,MF_texp,MF_airm,MF_totl)
			print('ImageInfo.base_radius')
			print(ImageInfo.base_radius)
			print('R1,R2,R3')
			print(self.R1,self.R2,self.R3)
	
	''' Aperture photometry. List pixels in each ring '''
	def apphot_pixels(self,X,Y):
		X = int(X+0.5); Y = int(Y+0.5)
		pixels_region = [[Xj,Yj] for Xj in xrange(X-self.R3-1,X+self.R3+1)\
			for Yj in xrange(Y-self.R3-1,Y+self.R3+1)]
		
		if DEBUG==True:
			print('pixels_region:')
			print(len(pixels_region))
			
		# Pixels in each ring
		def less_distance(Xi,Yi,reference):
			return (Xi-X)**2 + (Yi-Y)**2 <= reference**2
		
		if DEBUG==True:
			print('R1,R2,R3')
			print(self.R1,self.R2,self.R3)
		
		self.pixels1 = [Pixel for Pixel in pixels_region if \
			less_distance(Pixel[0],Pixel[1],self.R1)]
		self.pixels2 = [Pixel for Pixel in pixels_region if \
			less_distance(Pixel[0],Pixel[1],self.R2) and\
			not less_distance(Pixel[0],Pixel[1],self.R1)]
		self.pixels3 = [Pixel for Pixel in pixels_region if \
			less_distance(Pixel[0],Pixel[1],self.R3) and\
			not less_distance(Pixel[0],Pixel[1],self.R2)]
		
		if DEBUG==True:
			print('X,Y')
			print(X,Y)
			print('pixels1,pixels2,pixels3')
			print(len(self.pixels1),len(self.pixels2),len(self.pixels3))
		
		try:
			assert(self.pixels1!=[] and self.pixels2!=[] and self.pixels3!=[])
		except:
			self.destroy=True
	
	''' Turn pixel coordinates into pixel values '''
	@staticmethod
	def pixel_list_values(pixel_list,fits_data):
		def valid_pixel(pixel,fits_data):
			try:
				fits_data[pixel[0],pixel[1]]
			except:
				return(False)
			else:
				return(True)
		
		pixel_values = np.array([fits_data[pixel[0],pixel[1]] \
			for pixel in pixel_list if valid_pixel(pixel,fits_data)])
		
		return(pixel_values)
	
	''' Background and Star Fluxes '''
	def measure_background(self,fits_region):
		if DEBUG==True:
			print('Measure background in:')
			print(fits_region)
			assert(fits_region!=[])
		
		self.skyflux     = np.median(fits_region)
		self.skyflux_err = np.std(fits_region)
		
	def measure_totalflux(self,fits_region):
		if DEBUG==True:
			print('Measure total flux')
		self.totalflux = np.sum(fits_region)
		
	def measure_starflux(self,pixelvalues_star,pixelvalues_background):
		if DEBUG==True:
			print('Measure star flux')
		self.measure_background(pixelvalues_background)
		self.measure_totalflux(pixelvalues_star)
		self.starflux = self.totalflux - len(pixelvalues_star)*self.skyflux
		self.starflux_err = len(pixelvalues_star)*self.skyflux_err
	
	@staticmethod
	def star_is_saturated(fits_region,max_value):
		try:
			assert(np.max(fits_region)<max_value)
			return False
		except:
			return True
	
	
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
		if DEBUG==True: print('Estimate centroid')
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
	
	def optimal_aperture_photometry(self,FitsImage,min_radius,max_radius):
		if DEBUG==True: print('Optimal aperture photometry')
		'''
		Optimize the aperture to minimize uncertainties and assert
		all flux is contained in R1
		'''
		radius = min_radius
		continue_iterate = True
		while continue_iterate==True and radius<max_radius:
			radius += 1
			self.apphot_pixels(self.Xcoord,self.Ycoord)
			old_starflux = self.starflux
			region_background = self.pixel_list_values(self.pixels3,FitsImage.fits_data)
			region_star = self.pixel_list_values(self.pixels1,FitsImage.fits_data)
			
			try:
				assert(region_background!=[] and region_star!=[])
			except:
				print('region_background==[] or region_star==[]')
				self.destroy=True
			else:
				if DEBUG==True:
					print('region_background,region_star')
					print(len(region_background),len(region_star))
				
				self.measure_starflux(region_star,region_background)
				if self.starflux < old_starflux*1.02:
					continue_iterate = False
		
		try:
			assert(self.destroy==False)
			assert(radius<max_radius)
		except:
			self.destroy=True
	
	'''Aperture photometry. The main function'''	
	def stellar_aperture_photometry(self,FitsImage,ImageInfo):
		if DEBUG==True: print('Stellar aperture photometry')
		
		# Define aperture disk for photometry
		self.estimate_photometric_radius(ImageInfo)
		if (DEBUG==True):
			print('X,Y')
			print(self.Xcoord,self.Ycoord)
		
		self.apphot_pixels(self.Xcoord,self.Ycoord)
		# Generate regions to measure fluxes.
		if (DEBUG==True):
			print('R1,R2,R3')
			print(self.R1,self.R2,self.R3)
			print('pixels1,pixels2,pixels3')
			print(len(self.pixels1),len(self.pixels2),len(self.pixels3))
		
		if(DEBUG==True and self.destroy==False):
			print('Im here')
		
		region_background = self.pixel_list_values(self.pixels3,FitsImage.fits_data)
		region_star = self.pixel_list_values(self.pixels1,FitsImage.fits_data)
		try:
			assert(self.destroy==False)
			assert(region_background!=[] and region_star!=[])
		except:
			if DEBUG==True:
				print('region_background:')
				print(len(region_background))
				print('region star:')
				print(len(region_star))
			self.destroy=True
		else:
			# Measure fluxes
			if (self.star_is_saturated(region_star,2**(ImageInfo.ccd_bits-1))):
				self.destroy=True
			self.measure_starflux(region_star,region_background)
		
		if(DEBUG==True and self.destroy==False):
			print('Im still alive')
			print('self.starflux,self.star_min_flux_detectable(ImageInfo)')
			print(self.starflux,self.star_min_flux_detectable(ImageInfo))
			print('R1/2,R2')
			print(self.R1/2,self.R2)
		
		# Test if its detectable
		try:
			assert(self.destroy==False)
			assert(self.starflux >= self.star_min_flux_detectable(ImageInfo))
		except:
			self.destroy=True
		else:
			# More accurate photometry
			self.optimal_aperture_photometry(FitsImage,self.R1/2,self.R2)
			if self.destroy==False:
				self.photometry_bouguervar()


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
		if DEBUG==True: print('Closing catalog file ...'),
		self.catalogfile.close()
		if DEBUG==True: print('OK')

class PhotometricCatalog():
	''' This class processes the entire catalog '''
	def __init__(self,ObsPyephem,ImageInfo,catalog_filename='ducati_catalog.tsv'):
		# Load the catalog
		RawCatalog = CatalogFile(catalog_filename)
		
		print('Looking for stars in the image and performing photometry ...')
		
		# Calculate some properties for each star
		self.ProcessedCatalog = []
		if DEBUG==True: print('Loading star properties')
		for Star in RawCatalog.CatalogLines:
			if DEBUG==True: print('Star basic definition')
			PhotometricStar_ = PhotometricStar(Star,ObsPyephem,ImageInfo)
			if DEBUG==True: print('imagedep_properties')
			PhotometricStar_.imagedep_properties(ObsPyephem,ImageInfo)
			try:
				assert(PhotometricStar.destroy==False)
			except:
				self.ProcessedCatalog.append(PhotometricStar_)
			else:
				PhotometricStar.destroy=True
		
		# Filter the catalog
		if DEBUG==True: print('Filtering catalog. Number of remaining stars: '),
		self.Stars = [Star for Star in self.ProcessedCatalog if Star.destroy==False]
		if DEBUG==True: print(len(self.Stars))
		RawCatalog.__del__()
	
	''' Detect and measure star fluxes in the image. Destroy bad stars'''
	def photometric_measures(self,FitsImage,ImageInfo):
		if DEBUG==True:
			print('Len(Stars) before photometric_measures')
			print(len(self.Stars))
			
		for Star in self.Stars:
			if Star.destroy==False:
				try:
					Star.stellar_aperture_photometry(FitsImage,ImageInfo)
					if DEBUG==True: print('Star analized')
				except:
					raise
					self.destroy=True
		
		self.Stars = [Star for Star in self.ProcessedCatalog if Star.destroy==False]
		if DEBUG==True:
			print('Len(Stars) after photometric_measures')
			print(len(self.Stars))
	
	def plot_starmap(self,ImageInfo,ObsPyephem):
		TheSkyMap = SkyMap(fits_data,MeasuredCatalog,ImageInfo,ObsPyephem)
		for Star in self.Stars:
			TheSkyMap.annotate_skymap(Star,R1,R2,R3)
		show_or_save_skymap()
	
	def __del__(self):
		if DEBUG==True: print('Deleting catalog ...'),
		del(self)
		if DEBUG==True: print('OK')
