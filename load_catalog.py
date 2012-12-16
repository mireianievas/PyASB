#!/usr/bin/env python

'''
Load Catalog file and make PyAstMon StarCatalog

This module loads the catalog file and returns 
an array of Star objects (StarCatalog).
____________________________

This module is part of the PyAstMonUCM project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

try:
	import ephem
	import CommonErrors, HeaderTest
	from astrometry import atmospheric_refraction,calculate_airmass
except:
	print 'One or more modules missing: pyfits,CommonErrors,HeaderTest'
	raise SystemExit

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyAstMonUCM project"
__credits__ = ["Miguel Nievas"]
__license__ = "GPL"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


class CatalogStar():
	''' Extract information for each Star in the Catalog '''
	
	def __init__(self,CatalogLine,ObsPyephem,ImageInfo):
		self.basic_properties(CatalogLine)
		self.imagedep_properties(ObsPyephem,ImageInfo)
		self.destroy = False
	
	def coord_pyephem_format(self,coord_str):
		# Return coordinate in pyepheem str
		coord_separated = coord_str.split()
		coord_pyephem = str(int(coord_separated[0]))+\
			':'+str(int(coord_separated[1]))+str(int(coord_separated[2]))
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

	def imagedep_properties(self,ObsPyephem,ImageInfo):
		'''
		Calculate star position and dynamic properties.
		Don't proceed if the star altitude < 0.
		Return updated StarCatalog (image dependent)
		'''
	
		def pyephem_declaration(self,ObsPyephem):
			# Define the star in Pyephem to make astrometric calculations
			pyephem_star = ephem.readdb('"'+self.name+'"'+",f|S|A0,"+self.RA1950+'|0'+\
				","+Star.DEC1950+'|0'+","+self.Vmag+',1950,0"')
			pyephem_star.compute(ObsPyephem)	
			return pyephem_star
	
		def set_actual_filter(self,ImageInfo):
			# Test which filter is in use
			used_filter = ImageInfo.Properties.used_filter
			if used_filter=="JohnsonU":
				self.FilterMag = self.Umag
				self.Color     = self.U_V
			elif used_filter=="JohnsonB":
				self.FilterMag = self.Umag
				self.Color     = self.U_V
			elif used_filter=="JohnsonV":
				self.FilterMag = self.Vmag
				self.Color     = 0.0
			elif used_filter=="JohnsonR":
				self.FilterMag = self.Rmag
				self.Color     = self.R_V
			elif used_filter=="JohnsonI":
				self.FilterMag = self.Imag
				self.Color     = self.I_V
			else:
				pass
		
		try:
			pyephem_star = PyephemDeclaration(StarCatalog[line])
			if float(pyephem_star.alt) > 0.0:
				# Real coordinates (from catalog)
				self.altit_real = float(pyephem_star.alt)
				self.zdist_real = 90.0-StarCatalog[line].altit_real
				self.azimuth    = float(pyephem_star.az)
				# Apparent coordinates in sky. Atmospheric refraction effect.
				self.altit_appa = atmospheric_refraction(StarCatalog[line].altit_real,'dir')
				self.zdist_appa = 90.0-StarCatalog[line].altit_apar
				self.airmass    = calculate_airmass(StarCatalog[line].altit_appa)
				# Photometric properties
				set_actual_filter(ImageInfo)
				# Image coordinates
				XYCoordinates = horiz2xy(StarCatalog[line].azimuth,\
					StarCatalog[line].altit_appa,ImageInfo)
				self.Xcoord = XYCoordinates[0]
				self.Ycoord = XYCoordinates[1]
			
			if self.Xcoord<0. or self.Ycoord<0.or self.Xcoord>ImageInfo.Properties.resolution[0] \
			or self.Ycoord>ImageInfo.Properties.resolution[1] or self.altit_real < 0.0:
				# Star doesn't fit in the image
				self.destroy = True
		except:
			self.destroy = True

class CatalogFile():
	''' Catalog with Star photometry in U,B,V,R,I bands '''
	def __init__(self,catalog_filename):
		# Open catalog file and return its content
		try:
			self.file = open(catalog_filename, 'r')
			self.CatalogLines = [ line[:-1] for line in self.file.readlines() ]
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
		self.file.close()
		print('OK')

class PhotometricCatalog():
	''' This class processes the entire catalog '''
	def __init__(self,ObsPyephem,ImageInfo,catalog_filename='ducati_catalog.tsv'):
		RawCatalog = CatalogFile(catalog_filename)
		# Calculate some properties for each star
		ProcessedCatalog = [CatalogStar(Star,ObsPyephem,ImageInfo) \
			for Star in RawCatalog.CatalogLines]
		# Filter the catalog
		self.Stars = [Star for Star in ProcessedCatalog if Star.destroy==False]
		RawCatalog.__del__()
	def __del__(self):
		print('Deleting catalog ...'),
		del(self)
		print('OK')