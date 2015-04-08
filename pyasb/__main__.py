#!/usr/bin/env python2

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
	#import gc
	import sys,os,inspect
	import signal
	import time
	
	from input_options import *
	from image_info import *
	from help import *
	from astrometry import *
	from star_calibration import *
	from load_fitsimage import *
	from bouguer_fit import *
	from sky_brightness import *
	from skymap_plot import *
	from cloud_coverage import *
	from write_summary import *
except:
	#raise
	print(str(inspect.stack()[0][2:4][::-1])+\
	 ': One or more modules missing')
	raise# SystemExit

config_file_default  = 'config.cfg'


'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~ Halt handler ~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''


def handler(signum, frame):
        print 'Signal handler called with signal', signum
        print "CTRL-C pressed"
	raise SystemExit
        #sys.exit(0)

signal.signal(signal.SIGTERM, handler)
signal.signal(signal.SIGINT, handler)


#@profile
class LoadImage(object):
	def __init__(self, InputOptions, ImageInfo, ConfigOptions, input_file=None):
		# Load Image file list
		if input_file == None:
			input_file = InputOptions.fits_filename_list[0]
		
		''' Load fits image '''
		self.FitsImage = FitsImage(input_file)
		# Local copy of ImageInfo. We will process it further.
		self.ImageInfo = ImageInfo
		self.ImageInfo.read_header(self.FitsImage.fits_Header)
		self.ImageInfo.config_processing_specificfilter(ConfigOptions)
	
		try:
			self.FitsImage.subtract_corners_background = True
			self.FitsImage.reduce_science_frame(\
				self.ImageInfo.darkframe,\
				self.ImageInfo.sel_flatfield,\
				MasterBias=None,\
				ImageInfo=self.ImageInfo)
		except:
			#raise
			print(inspect.stack()[0][2:4][::-1])
			#raise
			print('Cannot reduce science frame')
		
		self.FitsImage.__clear__()
		self.output_paths(InputOptions)
	

	def output_paths(self,InputOptions):
		# Output file paths (NOTE: should be moved to another file or at least separated function)
		# Photometric table
		
		path_list = [\
			"photometry_table_path", "skymap_path", "bouguerfit_path", \
			"skybrightness_map_path", "skybrightness_table_path", \
			"cloudmap_path", "clouddata_path", "summary_path"]
		
		for path in path_list:
			try: setattr(self.ImageInfo,path,getattr(InputOptions,path))
			except:
				try: getattr(InputOptions,path)
				except: 
					setattr(self.ImageInfo,path,False)
		
#@profile
class ImageAnalysis():
	def __init__(self,Image):
		''' Analize image and perform star astrometry & photometry. 
		    Returns ImageInfo and StarCatalog'''
		Image.ImageInfo.catalog_filename = 'catalog.csv'
		self.StarCatalog = StarCatalog(Image.ImageInfo)
		
		if (Image.ImageInfo.calibrate_astrometry==True):
			TheSkyMap = SkyMap(Image.ImageInfo,Image.FitsImage)
			TheSkyMap.setup_skymap()
			TheSkyMap.set_starcatalog(self.StarCatalog)
			TheSkyMap.astrometry_solver()
		
		self.StarCatalog.process_catalog_specific(Image.FitsImage,Image.ImageInfo)
                self.StarCatalog.save_to_file(Image.ImageInfo)
		TheSkyMap = SkyMap(Image.ImageInfo,Image.FitsImage)
		TheSkyMap.setup_skymap()
		TheSkyMap.set_starcatalog(self.StarCatalog)
		TheSkyMap.complete_skymap()

'''#@profile
class MultipleImageAnalysis():
	def __init__(self,InputOptions):
		class StarCatalog_():
			StarList = []
			StarList_woPhot = []
		
		InputFileList = InputOptions.fits_filename_list
		
		for EachFile in InputFileList:
			EachImage = LoadImage(EachFile)
			EachAnalysis = ImageAnalysis(EachImage)
			self.StarCatalog.StarList.append(EachAnalysis.StarCatalog.StarList)
			self.StarCatalog.StarList_woPhot.append(EachAnalysis.StarCatalog.StarList_woPhot)
'''

#@profile
class InstrumentCalibration():
	def __init__(self,ImageInfo,StarCatalog):
		try:
			self.BouguerFit = BouguerFit(ImageInfo,StarCatalog)
		except Exception as e:
			print(inspect.stack()[0][2:4][::-1])
			print('Cannot perform the Bouguer Fit. Error is: ')
			print type(e)
			print e
			exit(0)
			#raise
		

#@profile
class MeasureSkyBrightness():
	def __init__(self,FitsImage,ImageInfo,BouguerFit):
		ImageCoordinates_ = ImageCoordinates(ImageInfo)
		SkyBrightness_ = SkyBrightness(FitsImage,ImageInfo,ImageCoordinates_,BouguerFit)
		SkyBrightnessGraph_ = SkyBrightnessGraph(SkyBrightness_,ImageInfo,BouguerFit)
		self.SBzenith = SkyBrightness_.SBzenith
		self.SBzenith_err = SkyBrightness_.SBzenith_err

#@profile
def perform_complete_analysis(InputOptions,ImageInfoCommon,ConfigOptions,input_file):
	# Load Image into memory & reduce it.
		# Clean (no leaks)
		Image_ = LoadImage(InputOptions,ImageInfoCommon,ConfigOptions,input_file)
		
		# Look for stars that appears in the catalog, measure their fluxes. Generate starmap.
		# Clean (no leaks)
		ImageAnalysis_ = ImageAnalysis(Image_)
		
		print('Image date: '+str(Image_.ImageInfo.date_string)+\
		 ', Image filter: '+str(Image_.ImageInfo.used_filter))
		
		'Create the needed classes for the summary write'
		class InstrumentCalibration_:
			class BouguerFit:
				class Regression:
					mean_zeropoint = -1
					error_zeropoint = -1
					extinction = -1
					error_extinction = -1
					Nstars_rel = -1
					Nstars_initial = -1
		
		try:
			# Calibrate instrument with image. Generate fit plot.
			# Clean (no leaks)
			InstrumentCalibration_ = InstrumentCalibration(\
				Image_.ImageInfo,
				ImageAnalysis_.StarCatalog)
		except:
			class ImageSkyBrightness:
				SBzenith = '-1'
				SBzenith_err = '-1'
			
		else:
			# Measure sky brightness / background. Generate map.
			ImageSkyBrightness = MeasureSkyBrightness(\
				Image_.FitsImage,
				Image_.ImageInfo,
				InstrumentCalibration_.BouguerFit)
		
		'''
		Even if calibration fails, 
		we will try to determine cloud coverage
		and write the summary
		'''
		
		# Detect clouds on image
		ImageCloudCoverage = CloudCoverage(\
			Image_,
			ImageAnalysis_,
			InstrumentCalibration_.BouguerFit)
		
		Summary_ = Summary(Image_, InputOptions, ImageAnalysis_, \
			InstrumentCalibration_, ImageSkyBrightness, ImageCloudCoverage)
		
		#gc.collect()
		#print(gc.garbage)


def get_config_filename(InputOptions):
	config_file = config_file_default
	try:
		assert(InputOptions.configfile!=False)
	except:
		print(str(inspect.stack()[0][2:4][::-1])+\
		 ': config file not specified, use the default one:')
	else:
		config_file = InputOptions.configfile
	
	return(config_file)


if __name__ == '__main__':
	#gc.set_debug(gc.DEBUG_STATS)
	PlatformHelp_ = PlatformHelp()
	InputOptions = ReadOptions(sys.argv)
	
	config_file = get_config_filename(InputOptions)
	ConfigOptions_ = ConfigOptions(config_file)
	ImageInfoCommon = ImageInfo()
	ImageInfoCommon.config_processing_common(ConfigOptions_,InputOptions)
	try:
		assert(InputOptions.show_help == False)
	except:
		print(inspect.stack()[0][2:4][::-1])
		# Show help and halt
		PlatformHelp_.show_help()
		raise SystemExit
	
	for input_file in InputOptions.fits_filename_list:
		perform_complete_analysis(InputOptions,ImageInfoCommon,ConfigOptions_,input_file)
	
	'''gc.collect()
	
	d = dict()
	for o in gc.get_objects():
		name = type(o).__name__
		if name not in d:
			d[name] = 1
		else:
			d[name] += 1

	items = d.items()
	items.sort(key=lambda x:x[1])
	debug_file = open("debug_objects.txt",'w')
	debug_file.close()
	debug_file = open("debug_objects.txt",'a+')
	for key, value in items:
		print key, value
		debug_file.write(str(key)+",\t"+str(value)+"\n")
	
	debug_file.close()
	'''
	
