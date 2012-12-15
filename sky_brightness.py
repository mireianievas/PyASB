#!/usr/bin/env python

'''
Sky Brightness photometry module

Measure Sky brightness from the image using previous 
instrument calibration. 

____________________________

This module is part of the PyAstMonUCM project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyAstMonUCM project"
__credits__ = ["Miguel Nievas"]
__license__ = "GPL"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


try:
	import matplotlib.pyplot as mpl
	import matplotlib.colors as mpc
	import matplotlib.patches as mpp
except:
	print 'One or more modules missing: pyfits,HeaderTest'
	raise SystemExit


def sky_brightness_measure(fits_data,ImageRegion,ImageInfo,Regression):
	'''
	Return measured sky brightness with its error at a given point in image.
	ImageRegion must contain a pixel list
	'''

	try:
		# Measure Sky fluxes
		sky_flux,sky_flux_err = sky_flux_measure(fits_data,ImageRegion,ImageInfo)
		# Compute Sky brightness in magnitudes
		sky_brightness,sky_brightness_err = \
			Regression.zp-2.5*log10(sky_flux/(ImageInfo.Properties.exposure*\
				ImageInfo.Config.pixel_scale)),\
			sqrt(Regression.zp_err**2 + (2.5*sky_flux_err/(log(10)*sky_flux))**2)
		
		return sky_brightness,sky_brightness_err
	except:
		raise

def sky_brightness_grid(fits_data,ImageInfo,Regression):
	'''
	Measure sky brightness at selected points in the sky
	Return the sky brightness table
	'''
	


	

def sky_brightness_map():
	'''
	Interpolate and plot the Sky Brightness map
	'''

