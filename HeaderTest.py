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
	import CommonErrors
except:
	print 'One or more modules missing: CommonErrors'
	raise SystemExit
	
__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyAstMonUCM project"
__credits__ = ["Miguel Nievas"]
__license__ = "GPL"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"

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

