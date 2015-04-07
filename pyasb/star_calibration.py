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
	import sys,os,inspect
	import ephem
	import scipy.stats
	from astrometry import *
	from skymap_plot import *
except:
	print(str(inspect.stack()[0][2:4][::-1])+': One or more modules missing')
	raise SystemExit


def verbose(function, *args):
	'''
	Run a function in verbose mode
	'''
	try:
		out = function(*args)
	except:
		# Something happened while runing function
		raise
		if DEBUG==True:
			print(str(inspect.stack()[0][2:4][::-1])+' Error')
			raise
	else:
		return(out)

class Star():
	def __init__(self,StarCatalogLine,ImageInfo):
		''' Takes StarCatalogLine (line from catalog file) and 
		      FitsImage, ImageInfo and ObsPyephem objects
		    Returns a Star object with photometric and astrometic properties 
		      or a destroy flag if errors ocurred during process'''
		self.destroy=False
		self.saturated=False
		self.cold_pixels=False
		self.masked=False
		self.to_be_masked=False
		self.camera_independent_astrometry(StarCatalogLine,ImageInfo)
	
	def camera_independent_astrometry(self,StarCatalogLine,ImageInfo):
		# Extract stars from Catalog
		self.verbose_detection(self.from_catalog,StarCatalogLine,\
		 errormsg=' Error extracting from catalog')
		# Estimate magnitude on the given image
		self.verbose_detection(self.magnitude_on_image,ImageInfo,\
		 errormsg=' Error reading used filter or magnitude>max')
		# Astrometry for the current star (sky)
		self.verbose_detection(self.star_astrometry_sky,ImageInfo,\
		 errormsg=' Error performing star astrometry (sky), Star not visible?')
	
	def camera_dependent_astrometry(self,FitsImage,ImageInfo):
		# Astrometry for the current star (image)
		self.verbose_detection(self.star_astrometry_image,ImageInfo,\
		 errormsg=' Error performing star astrometry (image)')
		# Estimate radius to do the aperture photometry
		self.verbose_detection(self.photometric_radius,ImageInfo,\
		 errormsg=' Error generating photometric radius')
		# Create regions of stars and star+background
		self.verbose_detection(self.estimate_fits_region_star,FitsImage,\
		 errormsg=' Cannot create the Star region')
		self.verbose_detection(self.estimate_fits_region_complete,FitsImage,\
		 errormsg=' Cannot create the Star+Background region')
	
	def camera_dependent_photometry(self,FitsImage,ImageInfo):
		# Measure fluxes
		self.verbose_detection(self.measure_star_fluxes,FitsImage.fits_data,\
		 errormsg=' Error measuring fluxes')
		# Check if star is detectable
		self.verbose_detection(self.star_is_detectable,ImageInfo,\
		 errormsg=' Star is not detectable')
		# Check star saturation
		self.verbose_detection(self.star_is_saturated,ImageInfo,\
		 errormsg=' Star is satured or has hot pixels')
		# Check cold pixels
		self.verbose_detection(self.star_has_cold_pixels,ImageInfo,\
		 errormsg=' Star has cold pixels')
		# Estimate centroid
		self.verbose_detection(self.estimate_centroid,ImageInfo,\
		 errormsg=' Star centroid calculated')
		# Update regions with new improved centroid.
		self.verbose_detection(self.estimate_fits_region_star,FitsImage,\
		 errormsg=' Cannot create the Star region')
		self.verbose_detection(self.estimate_fits_region_complete,FitsImage,\
		 errormsg=' Cannot create the Star+Background region')
	
	def check_star_issues(self,FitsImage,ImageInfo):
		# Check if star region is masked
		self.verbose_detection(self.star_region_is_masked,FitsImage,\
		 errormsg=' Star is masked')
		# Optimal aperture photometry
		#self.verbose_detection(\
		# self.optimal_aperture_photometry,ImageInfo,FitsImage.fits_data,\
		# errormsg=' Error doing optimal photometry')
		self.verbose_detection(\
		 self.measure_star_fluxes,FitsImage.fits_data,\
		 errormsg=' Error doing optimal photometry')
		# Check if star is detectable (with optimal astrometry)
		self.verbose_detection(self.star_is_detectable,ImageInfo,\
		 errormsg=' Star is not detectable')
		# Calculate Bouguer variables
		self.verbose_detection(self.photometry_bouguervar,ImageInfo,\
		 errormsg=' Error calculating bouguer variables')
		# Append star to star mask
		#self.verbose_detection(self.append_to_star_mask,FitsImage,\
		# errormsg=' Cannot add star to mask')
		if self.to_be_masked==True:
			try: self.append_to_star_mask(FitsImage)
			except: print(str(inspect.stack()[0][2:4][::-1])+' Cannot add star to mask')
		
	def clear_objects(self):
		if DEBUG==True:
			print('Clearing Object')
		self.__clear__()
		
		if self.destroy==False:
			if DEBUG==True:
				print(self.name + ' DONE. Star detected.')
		
	
	def verbose_detection(self,function, *args, **kwargs):
		'''
		Exec a detection step in verbose mode
		'''
		if self.destroy==False:
			function(*args)
			if self.destroy==True:
				if DEBUG==True:
					print(str(inspect.stack()[0][2:4][::-1])+str(function)+kwargs['errormsg'])
	
	def from_catalog(self,CatalogLine):
		''' Populate class with properties extracted from catalog:
		    recno, HDcode, RA1950, DEC1950, Vmag, U_V, B_V, R_V, I_V '''
		
		def coord_pyephem_format(coord_str):
			# Return coordinate in pyepheem str
			while coord_str[0]==' ':
				coord_str = coord_str[1:]
			
			coord_separated = coord_str.split(' ')
			coord_pyephem = str(int(coord_separated[0]))+\
				':'+str(int(coord_separated[1]))+\
				":"+str(float(coord_separated[2]))
			return coord_pyephem
		
		def get_float(value):
			'''
			Try to get the magnitude of the star,
			If it is missing, then flat it as Incomplete Photometry
			'''
			try:
				return(float(value))
			except:
				self.IncompletePhot = True
				return(0)
		
		def star_is_photometric(self):
			'''
			Flag the star for its photometry usefulness.
			It must have a complete photometric magnitudes
			and not to be double, variable or with 
			[manual flag] bad photometric properties
			'''
			
			self.PhotometricStandard=True
                        # Check few variables
                        if self.isDouble:       self.PhotometricStandard = False
                        if self.isVariab:       self.PhotometricStandard = False
                        if self.isBadPhot:      self.PhotometricStandard = False
                        if self.IncompletePhot: self.PhotometricStandard = False
		
			# Also, if colors are very blue or very red, discard them
                        if self.U_V<0:   self.PhotometricStandard = False
                        if self.B_V>1.5: self.PhotometricStandard = False
		
		self.IncompletePhot = False
		try:
			self.recno     = int(CatalogLine[0])
			self.HDcode    = str(CatalogLine[1]).replace(' ','')
			self.RA2000    = coord_pyephem_format(CatalogLine[2])
			self.DEC2000   = coord_pyephem_format(CatalogLine[3])
			self.RA1950    = coord_pyephem_format(CatalogLine[4])
			self.DEC1950   = coord_pyephem_format(CatalogLine[5])
			self.Vmag      = get_float(CatalogLine[6])
			self.U_V       = get_float(CatalogLine[7])
			self.B_V       = get_float(CatalogLine[8])
			self.R_V       = get_float(CatalogLine[9])
			self.I_V       = get_float(CatalogLine[10])
			self.isDouble  = str(CatalogLine[11]).replace(' ','')=="D"
			self.isVariab  = str(CatalogLine[12]).replace(' ','')=="V"
			self.r_SpTy    = str(CatalogLine[13]).replace(' ','')
			self.SpType    = str(CatalogLine[14]).replace(' ','')
			self.isBadPhot = str(CatalogLine[15]).replace(' ','')=="*"
			
			try:
				# Try to find the common name
				self.name = str(CatalogLine[16])
			except:
				# Use the HDcode as name
				self.name = self.HDcode
			
			#self.name = self.HDcode
			
			star_is_photometric(self)
			
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
			self.FilterMag = self.Bmag
			self.Color     = self.B_V
		elif ImageInfo.used_filter=="Johnson_V":
			self.FilterMag = self.Vmag
			self.Color     = 0.0
		elif ImageInfo.used_filter=="Johnson_R":
			self.FilterMag = self.Rmag
			self.Color     = self.R_V
		elif ImageInfo.used_filter=="Johnson_I":
			self.FilterMag = self.Imag
			self.Color     = self.I_V
		else: 
			self.destroy=True
		
		try:
			assert(self.FilterMag<ImageInfo.max_magnitude)
		except:
			self.destroy=True
	
	def star_astrometry_sky(self,ImageInfo):
		''' Perform astrometry. Returns (if star is visible and well defined) its position on the sky and image'''
		
		ObsPyephem = pyephem_setup_real(ImageInfo)
		
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
			# The catalog is defined for B1950, get the current coordinates
			ra  = float(pyephem_star.a_ra)*12./np.pi
			dec = float(pyephem_star.a_dec)*180./np.pi
			
			# Get the horizontal coordinates
			self.azimuth,self.altit_real = eq2horiz(ra,dec,ImageInfo)
			
			try:
				assert(self.altit_real)>float(ImageInfo.min_altitude)
			except:
				self.destroy=True
			else:
				self.zdist_real = 90.0-self.altit_real
		
		if self.destroy==False:
			# Apparent coordinates in sky. Atmospheric refraction effect.
			self.altit_appa = atmospheric_refraction(self.altit_real,'dir')
			try:
				assert(self.altit_appa)>float(ImageInfo.min_altitude)
			except:
				self.destroy=True
			else:
				self.zdist_appa = 90.0-self.altit_appa
				self.airmass    = calculate_airmass(self.altit_appa)
	
	def star_astrometry_image(self,ImageInfo):
		if self.destroy==False:			
			# Get the X,Y image coordinates
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
			self.R2 = self.R1*2.0
			self.R3 = self.R1*3.5
		except:
			self.destroy=True
	
	def estimate_fits_region_star(self,FitsImage):
		''' Return the region that contains the star 
		(both for calibrated and uncalibrated data)'''
		self.fits_region_star = [[FitsImage.fits_data[y,x] \
			for x in xrange(int(self.Xcoord - self.R1 + 0.5),\
			 int(self.Xcoord + self.R1 + 0.5))] \
			for y in xrange(int(self.Ycoord - self.R1 + 0.5),\
			 int(self.Ycoord + self.R1 + 0.5))]
		# We will need this to look for saturated pixels.
		self.fits_region_star_uncalibrated = [[FitsImage.fits_data_notcalibrated[y,x] \
			for x in xrange(int(self.Xcoord - self.R1 + 0.5),\
			 int(self.Xcoord + self.R1 + 0.5))] \
			for y in xrange(int(self.Ycoord - self.R1 + 0.5),\
			 int(self.Ycoord + self.R1 + 0.5))]
		
		# We have computed the star region. Flag it to be masked
		self.to_be_masked=True
	
	def estimate_fits_region_complete(self,FitsImage):
		''' Return the region that contains the star+background '''
		self.fits_region_complete = [[FitsImage.fits_data[y,x] \
			for x in xrange(int(self.Xcoord - self.R3 + 0.5),\
			 int(self.Xcoord + self.R3 + 0.5))] \
			for y in xrange(int(self.Ycoord - self.R3 + 0.5),\
			 int(self.Ycoord + self.R3 + 0.5))]
	
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
				if less_distance(x-len(self.fits_region_complete)/2.,\
				 y-len(self.fits_region_complete[0])/2.,self.R1)]
			
			self.pixels2 = [self.fits_region_complete[y][x] \
				for y in xrange(len(self.fits_region_complete))\
				for x in xrange(len(self.fits_region_complete[0])) \
				if less_distance(x-len(self.fits_region_complete)/2.,\
				 y-len(self.fits_region_complete[0])/2.,self.R2) and\
				not less_distance(x-len(self.fits_region_complete)/2.,\
				 y-len(self.fits_region_complete[0])/2.,self.R1)]
			
			self.pixels3 = [self.fits_region_complete[y][x] \
				for y in xrange(len(self.fits_region_complete))\
				for x in xrange(len(self.fits_region_complete[0])) \
				if less_distance(x-len(self.fits_region_complete)/2.,\
				 y-len(self.fits_region_complete[0])/2.,self.R3) and\
				not less_distance(x-len(self.fits_region_complete)/2.,\
				 y-len(self.fits_region_complete[0])/2.,self.R2)]
			
			# Sky background flux. t_student 95%.
			t_skyflux = scipy.stats.t.isf(0.025,np.size(self.pixels3))
			#self.skyflux = 2.5*np.median(self.pixels3)-1.5*np.mean(self.pixels3)
			self.skyflux = np.median(self.pixels3)
			self.skyflux_err = t_skyflux*np.std(self.pixels3)/np.sqrt(np.size(self.pixels3))
			# Sky background + Star flux
			totalflux = np.sum(self.pixels1)
			# Only star flux. 
			self.starflux = totalflux - np.size(self.pixels1)*self.skyflux
			self.starflux_err = np.sqrt(2)*np.size(self.pixels1)*self.skyflux_err
		except:
			self.destroy=True
	
	def star_region_is_masked(self,FitsImage):
		''' Check if the star is in the star mask'''
		self.masked = False
		for x in xrange(int(self.Xcoord - self.R1 + 0.5),int(self.Xcoord + self.R1 + 0.5)):
			for y in xrange(int(self.Ycoord - self.R1 + 0.5),int(self.Ycoord + self.R1 + 0.5)):
				if FitsImage.star_mask[y][x] == True:
					self.masked=True
					self.destroy = True
					return(0)

	def star_is_saturated(self,ImageInfo):
		''' Return true if star has one or more saturated pixels 
		    requires a defined self.fits_region_star'''
		try:
			assert(np.max(self.fits_region_star_uncalibrated)<(2./3)*2**ImageInfo.ccd_bits)
		except:
			#self.destroy=True
			self.PhotometricStandard=False
			self.saturated=True
		else:
			self.saturated=False
	
	def star_has_cold_pixels(self,ImageInfo):
		''' Return true if star has one or more cold (0 value) pixels 
		    requires a defined self.fits_region_star'''
		try:
			assert(np.max(self.fits_region_star_uncalibrated)>0)
		except:
			#self.destroy=True
			self.PhotometricStandard=False
			self.cold_pixels=True
		else:
			self.cold_pixels=False
	
	def star_is_detectable(self,ImageInfo):
		''' Set a detection limit to remove weak stars'''
		''' Check if star is detectable '''
		
		try:
			assert(self.starflux>\
			 ImageInfo.baseflux_detectable*self.starflux_err+1e-10)
		except:
			self.destroy=True
	
        def estimate_centroid(self,iterations=1):
		''' Returns star centroid from a region that contains the star
		    needs self.R1'''
		
		try:
			exponent = 2 # Values > 1 intensify the convergence of the method.
			xweight = np.array([range(1,len(self.fits_region_complete[0])+1)]*\
			 len(self.fits_region_complete))
			yweight = np.array([range(1,len(self.fits_region_complete)+1)]*\
			 len(self.fits_region_complete[0])).transpose()
			xcentroid = np.sum(xweight*np.power(\
			 self.fits_region_complete-self.skyflux,exponent))\
			 /np.sum(np.power(self.fits_region_complete-\
			 self.skyflux,exponent))
			ycentroid = np.sum(yweight*np.power(\
			 self.fits_region_complete-self.skyflux,exponent))\
				/np.sum(np.power(self.fits_region_complete-self.skyflux,exponent))
			self.Xcoord += xcentroid - len(self.fits_region_complete[0])/2
			self.Ycoord += ycentroid - len(self.fits_region_complete)/2
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
				if self.starflux < (1+0.002*num_iterations**2)*old_starflux:
					iterate=False
				else:
					radius+=1
				
				assert(radius<self.R2)
		except:
			self.destroy=True
	
	def photometry_bouguervar(self,ImageInfo):
		# Calculate parameters used in bouguer law fit
		try:
			_25logF      = 2.5*np.log10(self.starflux/ImageInfo.exposure)
			_25logF_unc  = (2.5/np.log(10))*self.starflux_err/self.starflux
			color_term = ImageInfo.color_terms[ImageInfo.used_filter][0]
			color_term_err = ImageInfo.color_terms[ImageInfo.used_filter][1]
			self.m25logF     = self.FilterMag+_25logF+color_term*self.Color
			self.m25logF_unc = np.sqrt(_25logF_unc**2 + (color_term_err*self.Color)**2)
		except:
			self.PhotometricStandard=False
			#self.destroy=True

	def append_to_star_mask(self,FitsImage):
		for x in xrange(int(self.Xcoord - self.R1 + 0.5),int(self.Xcoord + self.R1 + 0.5)):
			for y in xrange(int(self.Ycoord - self.R1 + 0.5),int(self.Ycoord + self.R1 + 0.5)):
				FitsImage.star_mask[y][x] = True
	
	def __clear__(self):
		backup_attributes = [\
		 "destroy","PhotometricStandard","HDcode","name","FilterMag",\
		 "Color","saturated","cold_pixels","masked","RA1950","DEC1950",\
		 "azimuth","altit_real","airmass","Xcoord","Ycoord",\
		 "R1","R2","R3","starflux","starflux_err","m25logF","m25logF_unc"\
		 ]
		for atribute in list(self.__dict__):
			#if atribute[0]!="_" and atribute not in backup_attributes:
			if atribute not in backup_attributes:
				del vars(self)[atribute]

class StarCatalog():
	''' This class processes the catalog.
	    Takes FitsImage,ImageInfo,ObsPyephem, returns an object with 
	    the processed Star list'''
	
	def __init__(self,ImageInfo):
		print('Creating Star Catalog ...')
		self.load_catalog_file(ImageInfo.catalog_filename)
		print('Star processing ...')
		self.process_catalog_general(ImageInfo)
		#self.process_catalog(FitsImage,ImageInfo)
		#self.save_to_file(ImageInfo)
		
	
	def load_catalog_file(self,catalog_filename):
		''' Returns Catalog lines from the catalog_filename '''
		
		def line_is_star(textline,sep=';'):
			theline = textline.split(sep)
			try:
				int(theline[0].replace(' ','')[0])
			except:
				return(False)
			else:
				return(True)
		
		def catalog_separation(textline,sep=';'):
			theline = textline
			theline = theline.replace('\r\n','')
			theline = theline.replace('\n','')
			theline = theline.split(sep)
			return(theline)
		
		try:
			self.catalogfile = open(catalog_filename, 'r')
			CatalogContent = self.catalogfile.readlines()
			MinLine = 1; separator = ";"
			self.CatalogLines = [catalog_separation(CatalogContent[line])
			 for line in xrange(len(CatalogContent)) \
			 if line>=MinLine-1 and line_is_star(CatalogContent[line])]
		except IOError:
			print('IOError. Error opening file '+catalog_filename+'.')
			#return 1
		except:
			print('Unknown error:')
			raise
			#return 2
		else:
			print('File '+str(catalog_filename)+' opened correctly.')
	
	def process_catalog_general(self,ImageInfo):
		'''
		Returns the processed catalog with 
		all the starts that should be visible.
		'''
		
		try: ImageInfo.max_star_number
		except: ImageInfo.max_star_number = len(self.CatalogLines)
		
		self.StarList_Tot = []
		for each_star in self.CatalogLines[0:ImageInfo.max_star_number]:
			TheStar = Star(each_star,ImageInfo)
			if (TheStar.destroy==False):
				self.StarList_Tot.append(TheStar)
		
		print(" - Total stars: %d" %len(self.StarList_Tot))
	
	def process_catalog_specific(self,FitsImage,ImageInfo):
		'''
		Returns the processed catalog with 
		all the starts that are detected.
		'''
		
		#Create the masked star matrix
		FitsImage.star_mask = np.zeros(np.array(FitsImage.fits_data).shape,dtype=bool)
		
		self.StarList_Det = []
		self.StarList_Phot = []
		for TheStar in self.StarList_Tot:
			TheStar.camera_dependent_astrometry(FitsImage,ImageInfo)
			TheStar.camera_dependent_photometry(FitsImage,ImageInfo)
			TheStar.check_star_issues(FitsImage,ImageInfo)
			TheStar.clear_objects()
			if (TheStar.destroy==False):
				self.StarList_Det.append(TheStar)
				if TheStar.PhotometricStandard==True:
					self.StarList_Phot.append(TheStar)
		
		print(" - Detected stars: %d" %len(self.StarList_Det))
		print(" - With photometry: %d" %len(self.StarList_Phot))
	
	def save_to_file(self,ImageInfo):
		try:
			assert(ImageInfo.photometry_table_path!=False)
		except:
			print('Skipping write photometric table to file')
		else:
			print('Write photometric table to file')
			
			content = ['#HDcode, CommonName, RA1950, DEC1950, Azimuth, '+\
			 'Altitude, Airmass, Magnitude, Color(#-V), StarFlux, StarFluxErr, '+\
			 'mag+2.5logF, [mag+2.5logF]_Err\n']
			for Star in self.StarList_Phot:
				content.append(str(Star.HDcode)+', '+str(Star.name)+', '+str(Star.RA1950)+', '+\
				 str(Star.DEC1950)+', '+str(Star.azimuth)+', '+str(Star.altit_real)+\
				 ', '+str(Star.airmass)+', '+str(Star.FilterMag)+', '+str(Star.Color)+\
				 ', '+str(Star.starflux)+', '+str(Star.starflux_err)+', '+str(Star.m25logF)+\
				 ', '+str(Star.m25logF_unc)+'\n')
			
			if ImageInfo.photometry_table_path == "screen":
				print(content)
			else:
				def phottable_filename(ImageInfo):
					filename = ImageInfo.photometry_table_path +\
						"/PhotTable_"+ImageInfo.obs_name+"_"+ImageInfo.fits_date+"_"+\
						ImageInfo.used_filter+".txt"
					return(filename)
				
				photfile = open(phottable_filename(ImageInfo),'w+')
				photfile.writelines(content)
				photfile.close()
			
		

