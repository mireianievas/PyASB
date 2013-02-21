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
	import gc
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
	raise# SystemExit


config_filename  = 'pyasb_config.cfg'


class LoadImage():
	def __init__(self,InputOptions,input_file=None):
		# Load Image file list
		if input_file == None:
			input_file = InputOptions.fits_filename_list[0]
		
		''' Load fits image '''
		self.FitsImage = FitsImage(input_file)
		self.ImageInfo = ImageInfo(self.FitsImage.fits_Header,config_filename)
		try:
			self.FitsImage.reduce_science_frame(\
				self.ImageInfo.darkframe,\
				self.ImageInfo.sel_flatfield,\
				MasterBias=None)
		except:
			raise
			print('Cannot reduce science frame')
		
		self.FitsImage.__clear__()
		self.output_paths(InputOptions)
		
	def output_paths(self,InputOptions):
		# Output file paths (NOTE: should be moved to another file or at least separated function)
		# Photometric table
		try: self.ImageInfo.photometry_table_path = InputOptions.photometry_table_path
		except:
			try: self.ImageInfo.photometry_table_path
			except: 
				raise
				SystemExit
		
		# Star Map
		try: self.ImageInfo.skymap_path = InputOptions.skymap_path
		except:
			try: self.ImageInfo.skymap_path
			except: 
				raise
				SystemExit
		
		# Bouguer Fit
		try: self.ImageInfo.bouguerfit_path = InputOptions.bouguerfit_path
		except:
			try: self.ImageInfo.bouguerfit_path
			except: 
				raise
				SystemExit
		
		# SkyBrightness
		try: self.ImageInfo.skybrightness_map_path = InputOptions.skybrightness_map_path
		except:
			try: self.ImageInfo.skybrightness_map_path
			except: 
				raise
				SystemExit
		
		try: self.ImageInfo.skybrightness_table_path = InputOptions.skybrightness_table_path
		except:
			try: self.ImageInfo.skybrightness_table_path
			except: 
				raise
				SystemExit
		
		# Summary
		try: self.ImageInfo.summary_path = InputOptions.summary_path
		except:
			try: self.ImageInfo.summary_path
			except: 
				raise
				SystemExit

class ImageAnalysis():
	def __init__(self,Image):
		''' Analize image and perform star astrometry & photometry. 
		    Returns ImageInfo and StarCatalog'''
		
		ObsPyephem_ = pyephem_setup(Image.ImageInfo)
		
		Image.ImageInfo.catalog_filename = 'ducati_catalog.tsv'
		self.StarCatalog = StarCatalog(Image.FitsImage,Image.ImageInfo,ObsPyephem_)
		
		SkyMap_ = SkyMap(self.StarCatalog,Image.ImageInfo,Image.FitsImage)

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

class InstrumentCalibration():
	def __init__(self,ImageInfo,StarCatalog):
		self.BouguerFit = BouguerFit(ImageInfo,StarCatalog)
		self.BouguerFit.bouguer_plot(ImageInfo)


class MeasureSkyBrightness():
	def __init__(self,FitsImage,ImageInfo,BouguerFit):
		ImageCoordinates_ = ImageCoordinates(ImageInfo)
		SkyBrightness_ = SkyBrightness(FitsImage,ImageInfo,ImageCoordinates_,BouguerFit)
		SkyBrightnessGraph_ = SkyBrightnessGraph(SkyBrightness_,ImageInfo,BouguerFit)
		self.SBzenith = SkyBrightness_.SBzenith
		self.SBzenith_err = SkyBrightness_.SBzenith_err

# NOTE: The following 2 functions should be moved to separate file or at least to a new class
# NOTE: Maybe should be rewrite as follows?:
# 1.) Create the file with the header
# 2.) Iterative add lines

def summarize_results(InputOptions, Image, ImageAnalysis, InstrumentCalibration, ImageSkyBrightness):
	sum_date   = str(Image.ImageInfo.fits_date)
	sum_filter = str(Image.ImageInfo.used_filter)
	sum_stars  = str(InstrumentCalibration.BouguerFit.Regression.Nstars_initial)
	sum_gstars = str("%.1f"%float(InstrumentCalibration.BouguerFit.Regression.Nstars_rel))
	sum_zpoint = \
		str("%.3f"%float(InstrumentCalibration.BouguerFit.Regression.mean_zeropoint))+' +/- '+\
		str("%.3f"%float(InstrumentCalibration.BouguerFit.Regression.error_zeropoint))
	sum_extinction = \
		str("%.3f"%float(InstrumentCalibration.BouguerFit.Regression.extinction))+' +/- '+\
		str("%.3f"%float(InstrumentCalibration.BouguerFit.Regression.error_extinction))
	sum_skybrightness = \
		str("%.3f"%float(ImageSkyBrightness.SBzenith))+' +/- '+\
		str("%.3f"%float(ImageSkyBrightness.SBzenith_err))
		
	return(sum_date, sum_filter,sum_stars, sum_gstars, sum_zpoint, sum_extinction, sum_skybrightness)

def save_summary_to_file(ImageInfo,summary_content):
	try:
		assert(ImageInfo.summary_path!=False)
	except:
		print('Skipping write summary to file')
	else:
		print('Write summary to file')
		def summary_filename(ImageInfo):
			filename = ImageInfo.photometry_table_path +\
				"/Summary_"+ImageInfo.obs_name+"_"+ImageInfo.fits_date+"_"+\
				ImageInfo.used_filter+".txt"
			return(filename)
		
		photfile = open(summary_filename(ImageInfo),'w+')
		
		content = ['#Date, Filter, Stars, % Good Stars, ZeroPoint, Extinction, SkyBrightness\n']
		for line in summary_content:
			content_line = ""
			for element in line:
				content_line += element + ", "
			content_line += "\n"
			content.append(content_line)
		
		photfile.writelines(content)
		photfile.close()

if __name__ == '__main__':
	PlatformHelp_ = PlatformHelp()
	InputOptions = ReadOptions(sys.argv)
	
	try:
		assert(InputOptions.show_help == False)
	except:
		# Show help and halt
		PlatformHelp_.show_help()
		raise SystemExit
	
	for input_file in InputOptions.fits_filename_list:
		# Load Image into memory & reduce it.
		Image_ = LoadImage(InputOptions,input_file)
		
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
		
		summary_content = summarize_results(InputOptions, Image_, ImageAnalysis_, \
			InstrumentCalibration_, ImageSkyBrightness)
		
		save_summary_to_file(Image_.ImageInfo, [summary_content])
		gc.collect()
	
