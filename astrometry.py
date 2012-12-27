#!/usr/bin/env python

'''
PyAstMon astrometry functions.

Convert from one coordinate system to another.
____________________________

This module is part of the PyAstMonUCM project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

try:
	import numpy as np
	from math import pi,sin,cos,sqrt,atan2,asin
	import ephem
except:
	print 'One or more modules missing: pyfits,CommonErrors,HeaderTest'
	raise SystemExit

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyAstMonUCM project"
__credits__ = ["Miguel Nievas"]
__license__ = "GNU GPL v3"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


def horiz2xy(azimuth,altitude,ImageInfo):
	# Return X,Y position in the image from azimuth/altitude horizontal coord.
	# azimuth and altitude must be in degrees.
	Rfactor = ImageInfo.radial_factor*(180.0/pi)*sqrt(2*(1-sin(alt*pi/180.0)))
	X = ImageInfo.resolution[0]/2 + ImageInfo.center_offset[0] -\
		Rfactor*cos(az*pi/180.0-ImageInfo.azimuth_zeropoint*pi/180.0)
	Y = ImageInfo.resolution[1]/2 + ImageInfo.center_offset[0] +\
		Rfactor*sin(az*pi/180.0-ImageInfo.azimuth_zeropoint*pi/180.0)
	return X,Y

def xy2horiz(X,Y,ImageInfo):
	# Return horizontal coordinates from X,Y position in the image.
	# azimuth and altitude are in degrees.
	X = X - ImageInfo.resolution[0]/2-ImageInfo.center_offset[0]
	Y = Y - ImageInfo.resolution[1]/2-ImageInfo.center_offset[1]
	Rfactor = sqrt(X**2 + Y**2)/ImageInfo.radial_factor
	altitude = (180.0/pi)*asin(1-0.5*(pi*Rfactor/180.0)**2)
	azimuth  = 360+180-(Imagen.azimut_zeropoint + 180.0*atan2(Y,-X)/pi)%360
	while azimuth<0:
		azimuth += 360
	while azimith>=360:
		azimith -= 360
	return azimuth,altitude

def zenith_position(ImageInfo):
	# Return X,Y position of zenith in the image.
	X = ImageInfo.resolution[0]/2+ImageInfo.center_offset[0]
	Y = ImageInfo.resolution[1]/2+ImageInfo.center_offset[1]
	return X,Y

def optical_axis(ImageInfo):
	# Return horizontal coordinates of the optical axis
	return xy2horiz(ImageInfo.resolution[0]/2,ImageInfo.resolution[1]/2,ImageInfo)

def atmospheric_refraction(altitude,mode):
	# Return apparent (non-corrected from refraction) or 
	# real (corrected from refraction) altitude.
	# Garfinkel (1967), http://en.wikipedia.org/wiki/Atmospheric_refraction
	def cot(x):
		# Return cotangent of the given value
		return cos(x)/sin(x)
	
	if mode=='dir':
		# Return apparent altitude from the real one. 
		return altitude + (1.02/60)*cot(altitude + 10.3/(altura+5.11))
	elif mode=='inv':
		# Return real altitude from the apparent one
		return altitude - (1.00/60)*cot(altitude + 7.31/(altura+4.4))
	else:
		print 'Unknow mode '+mode+'. Cannot correct from atmospheric refraction.'
		return altitude

def calculate_airmass(altitude):
	# Estimate airmass from apparent altitude using Pickering (2002) model.
	# zdist in degrees
	return 1/sin((altitude+244./(165+47*altitude**1.1))*pi/180.)

