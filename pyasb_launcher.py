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
	#from read_config import *
	from program_help import *
	from astrometry import *
	from star_calibration import *
	from load_fitsimage import *
	from bouguer_fit import *
	from sky_brightness import *
	from skymap_plot import *
except:
	print 'One or more modules missing: please check'
	raise


config_filename  = 'pyasb_config.cfg'
catalog_filename = 'ducati_catalog.tsv' 

class InstrumentCalibration():
	def __init__(self,InputFile):
		#ConfigOptions = ConfigOptions(config_name)
		FitsImage_ = FitsImage(InputFile)
		ImageInfo_ = ImageInfo(FitsImage_.fits_Header,config_filename)
		FitsImage_.reduce_science_frame(ImageInfo_.darkframe,ImageInfo_.sel_flatfield,MasterBias=None)
		
		ObsPyephem_ = pyephem_setup(ImageInfo_)
		PlatformHelp_ = PlatformHelp()
		
		PhotometricCatalog_ = PhotometricCatalog(ObsPyephem_,ImageInfo_,catalog_filename)
		PhotometricCatalog_.photometric_measures(FitsImage_,ImageInfo_)
		
		BouguerFit_ = BouguerFit(ImageInfo_,PhotometricCatalog_)
		BouguerFit_.bouguer_plot(ImageInfo_,ObsPyephem_)
		
		SkyBrightness_ = SkyBrightness(FitsImage_.fits_Data,ImageInfo_,BouguerFit_.Regression_)
		SkyBrightnessGraph_ = SkyBrightnessGraph(SkyBrightness_,ImageInfo_,ObsPyephem_,BouguerFit_.Regression_)
		
if __name__ == '__main__':
	Imag = InstrumentCalibration('/home/minaya/facultad/beca_mec/coordenadas/Jonhson_B20100914_000305.fit')
		
		