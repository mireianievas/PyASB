#!/usr/bin/env python

'''
PyASB launcher module

Concatenate processes
____________________________

This module is part of the PyASB project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

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
	from read_config import *
	from program_help import *
	from astrometry import *
	from star_calibration import *
	from load_fitsimage import *
	from bouguer_fit import *
	from sky_brightness import *
	from skymap_plot import *
except:
	print 'One or more modules missing: please check'
	raise SystemExit


config_filename = 'pyasb_config.txt'

class InstrumentCalibration():
	def __init__(Image):
		#ConfigOptions = ConfigOptions(config_name)
		FitsImage = FitsImage(InputFile)
		ImageInfo = ImageInfo(fits_Header,config_filename)
		PlatformHelp = PlatformHelp()
		
		
		
		