#!/usr/bin/env python

'''
Star detection module

Fit fluxes and star data to an extinction law to obtain 
extinction and instrument zeropoint.
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
	from skymap_plot import * 
except:
	print 'One or more modules missing: pyfits,HeaderTest'
	raise SystemExit


def bouguer_fit(StarMeasured, ImageInfo, ObsPyephem):
	''' Fit measured fluxes to an extinction model
		Return regression parameters (ZeroPoint, Extinction) '''
	
	class Regression:
		x    = [Star.airmass for Star in StarMeasured]
		y    = [Star.m25logF for Star in StarMeasured]
		yerr = [Star.m25logF_unc for Star in StarMeasured]
	
	try:
		class fixed_zp:
			fixed_y     = ImageInfo.Config.y
			fixed_y_unc = ImageInfo.Config.y_unc
		Regression = theil_sen(Regression,fixed_zp)
	except:
		Regression = theil_sen(Regression)

			
	


