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

import pyfits

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyAstMonUCM project"
__credits__ = ["Miguel Nievas"]
__license__ = "GPL"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"

def load_image(input_file,DataHeader):
	# Image loading function
	try: 
		file_opened = pyfits.open(fichero_imagen)
	except IOError:	
		print 'IOError. Error opening file '+fichero_imagen+'.'
		return 1
	except:
		print 'Unknown error:'
		raise
		return 2
	else:
		print 'File '+str(input_file)+' opened correctly.'
		return file_opened

def image_extract_info(file_header):
	'''
	Extract some data from the image header
	and put it in a class ImageInfo
	'''
	class ImageInfo:
		class Fecha:
			# Date and time in different formats
			fits_date	= HeaderTest.correct_date(file_header)
			date_array	= [fits_date[0:4],fits_date[4:6],fits_date[6:8],
						  fits_date[9:11],fits_date[11:13],fits_date[13:15]]
			date_string	= date_array[0]+"/"+date_array[1]+"/"+date_array[2]+" "+\
						  date_array[3]+":"+date_array[4]+":"+date_array[5]
		class Properties:
			# Exposure, resolution, filter
			exposure	= HeaderTest.correct_exposure(file_header)
			resolution	= HeaderTest.correct_resolution(file_header)
			used_filter = HeaderTest.correct_filter(file_header)
	
	return ImageInfo

class CommonErrors:
	def name_error(variable): 
		print 'Variable '+variable+' not defined.'
	def value_error():
		print 'Value type error.'


class HeaderTest:
	def correct_exposure(file_header):
		# Exposure
		try:
			texp = float(file_header['EXPOSURE'])
		except NameError:
			CommonErrors.name_error('EXPOSURE')
		except ValueError:
			CommonErrors.value_error()
		else:
			assert texp>0.,\
				'0s exposure time detected.'
			return texp
	
	def correct_date(file_header):
		# Date and time
		try:
			date = file_header['DATE']
		except NameError:
			CommonErrors.name_error('DATE')
		else:
			assert len(date)==6,\
				'Date format not YYYYMMDD'
			return date	
	
	def correct_resolution(file_header):
		# Resolution
		try:
			resolution = [int(file_header['NAXIS1']),int(file_header['NAXIS2'])]
		except NameError:
			CommonErrors.name_error('NAXIS1 or NAXIS2')
		else:
			assert resolution[0]>0 and resolution[1]>0,\
				'Matrix not 2 dimensional'
			return resolution
	
	def correct_filter(file_header):
		# Test if there's a known filter
		try:
			used_filter = file_header['FILTER']
		except NameError:
			CommonErrors.name_error('FILTER')
		else:
			# due to an inconsistent format in AstMon, 
			# we found 4 possible formats 'Jonhson_V','JohnsonV','Johnson_V','JonhsonV'
			used_filter = used_filter.replace('_','')
			assert used_filter[0:7] in ['Johnson','Jonhson'],\
				'Filter type not recognized as Johnson system'
			assert used_filter[7:] in ['U','B','V','R','I'],\
				'Filter not U,B,V,R or I'
			return 'Johnson_'+used_filter[7:]






	


