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
	import sys
	from input_options import *
	from program_help import *
	from astrometry import *
	from star_calibration import *
	from load_fitsimage import *
	from bouguer_fit import *
	from sky_brightness import *
	from skymap_plot import *
	import time
except:
	print(str(sys.argv[0]) + ': One or more modules missing: please check')
	raise SystemExit


config_filename  = 'pyasb_config.cfg'


@profile
class LoadImage():
	def __init__(self,InputFile):
		''' Load fits image '''
		self.FitsImage = FitsImage(InputFile)
		self.ImageInfo = ImageInfo(self.FitsImage.fits_Header,config_filename)
		try:
			self.FitsImage.reduce_science_frame(\
				self.ImageInfo.darkframe,\
				self.ImageInfo.sel_flatfield,\
				MasterBias=None)
		except:
			print('Cannot reduce science frame')

@profile
class ImageAnalysis():
	def __init__(self,Image):
		''' Analize image and perform star astrometry & photometry. 
		    Returns ImageInfo and StarCatalog'''
		
		ObsPyephem_ = pyephem_setup(Image.ImageInfo)
		
		Image.ImageInfo.catalog_filename = 'ducati_catalog.tsv'
		Image.ImageInfo.skymap_file = "test_map.png"
		self.StarCatalog = StarCatalog(Image.FitsImage,Image.ImageInfo,ObsPyephem_)
		
		print('Star Map plot ...'),
		SkyMap_ = SkyMap(self.StarCatalog,Image.ImageInfo,Image.FitsImage)
		print('OK')

@profile
class MultipleImageAnalysis():
	def __init__(self,InputFileList):
		class StarCatalog_():
			StarList = []
			StarList_woPhot = []
		
		for EachFile in InputFileList:
			EachImage = LoadImage(EachFile)
			EachAnalysis = ImageAnalysis(EachImage)
			self.StarCatalog.StarList.append(EachAnalysis.StarCatalog.StarList)
			self.StarCatalog.StarList_woPhot.append(EachAnalysis.StarCatalog.StarList_woPhot)

@profile
class InstrumentCalibration():
	def __init__(self,ImageInfo,StarCatalog):
		print('Calculating Instrument zeropoint and extinction ...'),
		self.BouguerFit = BouguerFit(ImageInfo,StarCatalog)
		self.BouguerFit.bouguer_plot(ImageInfo)
		print('OK')

@profile
class MeasureSkyBrightness():
	def __init__(self,FitsImage,ImageInfo,BouguerFit):
		print('Generating Sky Brightness Map ...'),
		ImageCoordinates_ = ImageCoordinates(ImageInfo)
		SkyBrightness_ = SkyBrightness(FitsImage,ImageInfo,ImageCoordinates_,BouguerFit)
		SkyBrightnessGraph_ = SkyBrightnessGraph(SkyBrightness_,ImageInfo,BouguerFit)
		print('OK')


if __name__ == '__main__':
	PlatformHelp_ = PlatformHelp()
	InputOptions = ReadOptions(sys.argv)
	
	try:
		assert(InputOptions.show_help == False)
	except:
		# Show help and halt
		PlatformHelp_.show_help()
		raise SystemExit
	
	# Load Image into memory & reduce it.
	Image_ = LoadImage(InputFile =InputOptions.fits_filename_list[0])
	
	# Look for stars that appears in the catalog, measure their fluxes. Generate starmap.
	ImageAnalysis_ = ImageAnalysis(Image_)
	
	# Calibrate instrument with image. Generate fit plot.
	InstrumentCalibration_ = InstrumentCalibration(\
		Image_.ImageInfo,
		ImageAnalysis_.StarCatalog)
	
	# Measure sky brightness / background. Generate map.
	ImageSkyBrightness = MeasureSkyBrightness(\
		Image_.FitsImage,
		Image_.ImageInfo,
		InstrumentCalibration_.BouguerFit)
	
	
	