#!/usr/bin/env python

'''
Load FITS image and header

This module loads the AllSky FITS image and returns both 
the Image binary data and the plain-text header.
____________________________

This module is part of the PyAstMonUCM project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

try:
	import pyfits
except:
	print 'One or more modules missing: pyfits,HeaderTest'
	raise SystemExit

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyAstMonUCM project"
__credits__ = ["Miguel Nievas"]
__license__ = "GPL"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"

class FitsImage():
	def __init__(self,input_file):
		
		self.load_image(input_file)
		
		
	def load_image(self,input_file):
		# Image loading function
		try: 
			self.file_opened = pyfits.open(fichero_imagen)
			self.file_data   = file_opened[0].data
			self.fits_header = file_opened[0].header
		except IOError:	
			print 'IOError. Error opening file '+fichero_imagen+'.'
			return 1
		except:
			print 'Unknown error:'
			raise
			return 2
		else:
			print 'File '+str(input_file)+' opened correctly.'
			return fits_data,fits_header
				
	def __del__(self):
		del(self)
	
class ImageInfo():
	'''
	Extract some data from the image header
	and put it in a class ImageInfo
	'''
	def __init__(self,fits_header):
		# Date and time in different formats
		self.fits_date	 = self.correct_date(fits_header)
		self.date_array	 = [self.fits_date[0:4],self.fits_date[4:6],self.fits_date[6:8],
			self.fits_date[9:11],self.fits_date[11:13],self.fits_date[13:15]]
		self.date_string = date_array[0]+"/"+date_array[1]+"/"+date_array[2]+" "+\
			self.date_array[3]+":"+self.date_array[4]+":"+self.date_array[5]

		# Exposure (float), resolution (2d int array), filter (str)
		self.exposure	 = self.correct_exposure(fits_header)
		self.resolution	 = self.correct_resolution(fits_header)
		self.used_filter = self.correct_filter(fits_header)
	
	def correct_exposure(self,file_header):
		# Exposure
		try: texp = float(file_header['EXPOSURE'])
		except: raise
		else:
			assert texp>0., '0s exposure time detected.'
			return texp
	
	def correct_date(self,file_header):
		# Date and time
		try: date = file_header['DATE']
		except:	raise
		else:
			assert len(date)==6, 'Date format not YYYYMMDD'
			return date	
	
	def correct_resolution(self,file_header):
		# Resolution
		try: resolution = [int(file_header['NAXIS1']),int(file_header['NAXIS2'])]
		except: raise
		else:
			assert resolution[0]>0 and resolution[1]>0, 'Matrix not 2 dimensional'
			return resolution
	
	def correct_filter(self,file_header):
		# Test if there's a known filter
		try: used_filter = file_header['FILTER']
		except: raise
		else:
			# due to an inconsistent format in AstMon, 
			# we found 4 possible formats 'Jonhson_V','JohnsonV','Johnson_V','JonhsonV'
			used_filter = used_filter.replace('_','')
			assert used_filter[0:7] in ['Johnson','Jonhson'], 'Filter type not Johnson'
			assert used_filter[7:] in ['U','B','V','R','I'], 'Filter not U,B,V,R or I'
			return 'Johnson_'+used_filter[7:]

	def __del__(self):
		del(self)
	
	
	