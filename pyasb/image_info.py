#!/usr/bin/env python

'''
Load FITS Info data

This module processes the Image metadata and the program main options
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
	import sys,os,inspect
	import numpy as np
	import astropy.io.fits as pyfits
	from read_config import *
	from load_fitsimage import ImageTest
	from astrometry import pyephem_setup_real, pyephem_setup_common
except:
	print(str(inspect.stack()[0][2:4][::-1])+': One or more modules missing')
	raise SystemExit


class ImageInfo(ImageTest):
	'''
	Extract some data from the image header
	and put it in a class ImageInfo
	'''
	def __init__(self):
		pass
	
	'''
	def __init__(self,fits_header,ConfigOptions):
		self.read_header(fits_header)
		self.config_processing(ConfigOptions)
	'''
	
	def read_header(self,fits_header):
		# Date and time in different formats
		self.fits_date	 = ImageTest.correct_date(fits_header)
		self.date_array	 = [self.fits_date[0:4],self.fits_date[4:6],self.fits_date[6:8],
			self.fits_date[9:11],self.fits_date[11:13],self.fits_date[13:15]]
		self.date_string = self.date_array[0]+"/"+self.date_array[1]+"/"+self.date_array[2]+" "+\
			self.date_array[3]+":"+self.date_array[4]+":"+self.date_array[5]
		
		ObsPyephem = pyephem_setup_real(self)
		self.local_sidereal_time =  ObsPyephem.sidereal_time()*12./np.pi
		ObsPyephem.lon=0
		self.sidereal_time = ObsPyephem.sidereal_time()*12./np.pi

		# Exposure (float), resolution (2d int array), filter (str)
		self.exposure	 = ImageTest.correct_exposure(fits_header)
		self.resolution	 = ImageTest.correct_resolution(fits_header)
		self.used_filter = ImageTest.correct_filter(fits_header)
	
	def config_processing_common(self,ConfigOptions,InputOptions):
		# Default values
		self.darkframe = False
		self.biasframe = False
		self.latitude_offset = 0.0
		self.longitude_offset = 0.0
		
		for atribute in list(InputOptions.__dict__):
			ConfigOptions.FileOptions.append([atribute,vars(InputOptions)[atribute]])
		
		# Config processing
		for option in ConfigOptions.FileOptions:
			if   option[0]=="obs_latitude":             self.latitude=float(option[1])
			elif option[0]=="obs_longitude":            self.longitude=float(option[1])
			elif option[0]=="obs_name":                 self.obs_name=str(option[1]).replace(" ","")
			elif option[0]=="latitude_offset":	    self.latitude_offset=float(option[1])
			elif option[0]=="longitude_offset":         self.longitude_offset=float(option[1])
			elif option[0]=="delta_x":                  self.delta_x=float(option[1])
			elif option[0]=="delta_y":                  self.delta_y=float(option[1])
			elif option[0]=="radial_factor":            self.radial_factor=float(option[1])
			elif option[0]=="azimuth_zeropoint":        self.azimuth_zeropoint=float(option[1])
			elif option[0]=="min_altitude":             self.min_altitude=float(option[1])
			elif option[0]=="base_radius":              self.base_radius=float(option[1])
			elif option[0]=="baseflux_detectable":      self.baseflux_detectable=float(option[1])
			elif option[0]=="lim_Kendall_tau":          self.lim_Kendall_tau=float(option[1])
			elif option[0]=="ccd_bits":                 self.ccd_bits=float(option[1])
			elif option[0]=="ccd_gain":                 self.ccd_gain=float(option[1])

			elif option[0]=="perc_low":                 self.perc_low=float(option[1])
			elif option[0]=="perc_high":                self.perc_high=float(option[1])
			
                        elif option[0]=="read_noise":               self.read_noise=float(option[1])
			elif option[0]=="thermal_noise":            self.thermal_noise=float(option[1])
			elif option[0]=="max_magnitude":            self.max_magnitude = float(option[1])
			elif option[0]=="max_star_number":          self.max_star_number = int(option[1])
                        
                        elif option[0]=="backgroundmap_title":      self.backgroundmap_title = str(option[1])
			elif option[0]=="cloudmap_title":           self.cloudmap_title = str(option[1])
			elif option[0]=="skymap_path":              self.skymap_path = str(option[1]).replace(" ","")
			elif option[0]=="photometry_table_path":    self.photometry_table_path = str(option[1]).replace(" ","")
			elif option[0]=="bouguerfit_path":          self.bouguerfit_path = str(option[1]).replace(" ","")
			elif option[0]=="skybrightness_map_path":   self.skybrightness_map_path = str(option[1]).replace(" ","")
			elif option[0]=="skybrightness_table_path": self.skybrightness_table_path = str(option[1]).replace(" ","")
			elif option[0]=="cloudmap_path":            self.cloudmap_path = str(option[1]).replace(" ","")
			elif option[0]=="clouddata_path":           self.clouddata_path = str(option[1]).replace(" ","")
			elif option[0]=="summary_path":             self.summary_path = str(option[1]).replace(" ","")
			elif option[0]=="darkframe":                self.darkframe = str(option[1]).replace(" ","")
			elif option[0]=="biasframe":                self.biasframe = str(option[1]).replace(" ","")
	
	def config_processing_specificfilter(self,ConfigOptions):
		filters=["U","B","V","R","I"]
		
		self.zero_points = {}
		self.color_terms = {}
		self.background_levels = {}
		self.flatfield = {}
		
		for the_filter in filters:
			self.zero_points["Johnson_"+the_filter] = [False,False]
			self.color_terms["Johnson_"+the_filter] = [0,0]
			self.background_levels["Johnson_"+the_filter] = [False,False]
			self.flatfield["Johnson_"+the_filter] = False
		
		# Options that depends on the filter
		for option in ConfigOptions.FileOptions:
			for the_filter in xrange(len(filters)):
				filter_name = "Johnson_"+filters[the_filter]
				if   option[0]=="zero_point_"+filters[the_filter]:
					self.zero_points[filter_name] = \
						[float(option[1].split(",")[0]), float(option[1].split(",")[1])]
					if "Johnson_"+filters[the_filter] == self.used_filter:
						self.used_zero_point = self.zero_points[filter_name]
				elif option[0]=="color_term_"+filters[the_filter]:
					self.color_terms[filter_name] = \
						[float(option[1].split(",")[0]), float(option[1].split(",")[1])]
					if "Johnson_"+filters[the_filter] == self.used_filter:
						self.sel_color_terms = self.color_terms[filter_name]
				elif option[0]=="bkgnd_minmax_"+filters[the_filter]:
					self.background_levels[filter_name] = [float(option[1].split(",")[0]), \
					float(option[1].split(",")[1])]
					if "Johnson_"+filters[the_filter] == self.used_filter:
						self.sel_background_levels = self.background_levels[filter_name]
				elif option[0]=="flatfield_"+filters[the_filter]:
					self.flatfield[filter_name] = str(option[1]).replace(" ","")
					if "Johnson_"+filters[the_filter] == self.used_filter:
						self.sel_flatfield = self.flatfield[filter_name]
		
		# Some hacks to improve detectability in specific filters
		
		if self.used_filter == 'Johnson_B':
			self.base_radius = self.base_radius*1.2
			self.baseflux_detectable = self.baseflux_detectable*1.2
		
		if self.used_filter == 'Johnson_U':
			self.base_radius = self.base_radius*2
			self.baseflux_detectable = self.baseflux_detectable*2
			self.max_magnitude = self.max_magnitude-0.5
