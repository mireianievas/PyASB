#!/usr/bin/env python

'''
PyASB help module

Build and show program help
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
__longname__ = "Python All-Sky Brightness measuring pipeline"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"

class PlatformHelp():
	def __init__():
		self.make_title()
		self.make_welcome()
		self.make_requisites()
		self.make_options()
	
	def make_title(self):
		nameversion = 5*'#'+__shortname__+' v'+__version__+5*'#'
		self.separator = len(nameversion)*'#'
		self.title = +'\n'+nameversion+'\n'+self.separator+'\n'
	
	def make_welcome(self):
		self.welcome = 'Welcome to '+__shortname__+' ('+longname+')\n'
	
	def make_requisites(self):
		self.requisites = 'Requisites: Python 2.7; Scipy; Numpy; Matplotlib; Pyfits; Pyephem\n'
	
	def make_options(self):
		self.options = \
			'
			-h: print this help message\n\n
			-i input_allsky_image: All Sky image you want to be analyzed\n\n
			-f [year,month,day]: Date to be analyzed (AstMon-UCM only)\n
			  month and day are optional\n\n
			-om output_map_image: Output star map image, full or relative path\n
			  if no output file, show the map on screen\n\n
			-ot output_photometric_table: Output photometric table full or relative path\n
			  if no output file, show the table on screen\n\n
			-or output_results_summary: Summary of analysis, fit parameters and zenith SB\n
			  full or relative path. If no output file, show the table on screen\n\n
			-ob output_bouguerfit_graph: Output bouguer-fit graph, full or relative path\n
			  if no output file, show the graph on screen\n\n
			-os output_skybrightness_graph: Output Sky Brightness graph, full or relative path\n
			  if no output file, show the graph on screen\n\n
			-ost output_skybrightness_table: Output Sky Brightness table, full or relative path\n
			  if no output file, show the graph on screen\n\n
			'
	
	def show_help():
		print(\
			self.separator+'\n'+self.title+self.separator+'\n'+self.welcome+\
			self.requisites+self.options+self.separator)
			
	def incorrect_parameter(self,parameter):
		print('ERROR. Incorrect parameter: '+str(parameter))
		self.show_help()
	
	def date_or_file_input_error(self):
		print('ERROR. Date or file input')
		self.show_help()
	
	def no_parameters_error(self):
		print('ERROR. No input parameters especified')
		self.show_help()
	
	