#!/usr/bin/env python
# NOTE: OBSOLETE !

'''
Calibrate FITS image with MasterDark and MasterFlat

This module loads the AllSky FITS image and returns calibrated 
image with MasterDark subtraction and normalized with MasterFlat.
____________________________

This module is part of the PyAstMonUCM project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

try:
	import pyfits
	import numpy as np
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


def scale_darkframe(MasterDark,MasterBias,science_exposure):
	'''
	If not using the same exposure between Dark and Science 
	frames, we need to scale darks.
	MasterBias can be both a Pyfits image or a float value.
	MasterDark is the Pyfits image to scale
	Science_exposure is the exposition time of the science frame in seconds.
	Returns scaled MasterDark
	'''
	
	master_dark_data	= MasterDark[0].data
	master_dark_header	= MasterDark[0].header
	master_dark_exposure = correct_exposure(master_dark_header)

	if type(MasterBias) == "pyfits.hdu.hdulist.HDUList":
		master_bias_data = MasterBias[0].data
	elif type(MasterBias) in ["float","int"]:
		master_bias_data = MasterBias*np.ones((len(master_dark_data),len(master_dark_data[0])))
	else:
		print 'MasterBias not a Pyfits image nor Int/Float value'
		return MasterDark
	
	dark_current = master_dark_data - master_bias_data
	dark_current_normalized = dark_current/master_dark_exposure
	dark_current_scaled = dark_current_normalized*science_exposure
	MasterDark[0].data = dark_current_scaled + master_bias_data
	
	return MasterDark

def subtract_masterdark(ScienceImage,MasterDark):
	'''
	Return DarkFrame (Dark current + Bias) subtracted ScienceImage.
	DarkFrame can be both Pyfits Image or Int/Float value
	'''
	if type(ScienceImage) == "pyfits.hdu.hdulist.HDUList":
		science_data = ScienceImage[0].data
	else:
		print 'ScienceImage not pyfits HDUList Image type'
		return ScienceImage
	
	if type(MasterDark) == "pyfits.hdu.hdulist.HDUList":
		master_dark_data = MasterDark[0].data
	elif type(MasterDark) in ["float","int"]:
		master_dark_data = 	MasterDark*np.ones((len(science_data),len(science_data[0])))
	else:
		print 'MasterDark not a Pyfits image nor Int/Float value'
		return ScienceImage

def flatfield_normalization(FlatField):
	# Return normalized FlatField
	if type(FlatField) == "pyfits.hdu.hdulist.HDUList":
		FlatField[0].data = FlatField[0].data/(np.mean(FlatField[0].data))
	else:
		print 'FlatField not a pyfits HDUList Image type'
	return FlatField

def divide_flatfield(ScienceImage,FlatField):
	'''
	Return Science framed divided with FlatField.
	FlatField must be a Pyfits Image
	'''
	if type(ScienceImage) == "pyfits.hdu.hdulist.HDUList" and \
			type(FlatField) == "pyfits.hdu.hdulist.HDUList":

		FlatField = flatfield_normalization(FlatField)
		
		science_data	= ScienceImage[0].data
		flatfield_data	= FlatField[0].data
		ScienceImage[0].data = science_data/flatfield_data
	else:
		print 'ScienceImage or FlatField not pyfits HDUList Image type'
	return ScienceImage

def calibrate_image(ScienceImage,MasterDark=None,FlatField=None,MasterBias=None,MasterDarkFlat=None):
	# Returns calibrated image with Dark, Flat, Bias
	if MasterDark!=None:
		if MasterBias!=None:
			science_exposure = correct_exposure(ScienceImage[0].data)
			MasterDark = scale_darkframe(MasterDark,MasterBias,science_exposure)
		ScienceImage = subtract_masterdark(ScienceImage,MasterDark)
	if FlatField !=None:
		if MasterDarkFlat!=None:
			if MasterBias!=None:
				flat_exposure = correct_exposure(FlatField[0].data)
				MasterDarkFlat = scale_darkframe(MasterDarkFlat,MasterBias,flat_exposure)
			FlatField = subtract_masterdark(FlatField,MasterDarkFlat)
		ScienceImage = divide_flatfield(ScienceImage,FlatField)
	return ScienceImage
