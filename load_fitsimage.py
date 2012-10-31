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

import sys
import pyfits
import numpy as np

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyAstMonUCM project"
__credits__ = ["Miguel Nievas"]
__license__ = "GPL"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Development", or "Production"

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
		# Date and time in different formats
		fits_date    = header_correct_date(file_header)
		date_array   = [fits_date[0:4],fits_date[4:6],fits_date[6:8],
						fits_date[9:11],fits_date[11:13],fits_date[13:15]]
		date_str     =  date_array[0]+"/"+date_array[1]+"/"+date_array[2]+" "+\
						date_array[3]+":"+date_array[4]+":"date_array[5]
		# Exposure
		texp         = header_correct_exposure(file_header)




class CommonErrors:
	def name_error(): 
		print 'Variable not defined.'
	def value_error():
		print 'Value type error.'


def header_correct_exposure(file_header):
	# Exposure
	try:
		texp = float(cabecera_fits['EXPOSURE'])
	except NameError:
		CommonErrors.name_error()
	except ValueError:
		CommonErrors.value_error()
	else:
		assert texp > 0.,
			'0s exposure time detected.'
		return texp

def header_correct_date(file_header):
	# Date and time
	try:
		date = cabecera_fits['DATE']
	except NameError:
		CommonErrors.name_error()
	else:
		assert len(date)==6,
			'Date format not YYYYMMDD'
		return date	

		 



	


