#!/usr/bin/env python

'''
Load FITS image and header

This module loads the AllSky FITS image and returns both 
the Image binary data and the plain-text header.
____________________________

This module is part of the PyASB project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

try:
	import pyfits
except:
	print 'One or more modules missing: pyfits,HeaderTest'
	raise SystemExit

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


class ImageTest:
	pass
class ConfigOptions:
	pass


class FitsImage(ImageTest):
	def __init__(self,input_file):
		self.load_science(input_file)
		
	def load_science(self,input_file):
		print('Loading ScienceFrame ...'),
		try: 
			file_opened = pyfits.open(input_file)
			self.fits_Data   = file_opened[0].data
			self.fits_Header = file_opened[0].header
			self.fits_Texp   = float(ImageTest.correct_exposure(self.fits_Header))
		except:	raise
		else: print('OK. File '+str(input_file)+' opened correctly.')
		
	def load_dark(self,MasterDark):
		print('Loading MasterDark ...'),
		try:
			MasterDark_HDU    = pyfits.open(MasterDark)
			self.MasterDark_Data   = MasterDark_HDU[0].data
			self.MasterDark_Header = MasterDark_HDU[0].header
			self.MasterDark_Texp   = float(ImageTest.correct_exposure(MasterDark_Header))
		except: raise
		else: print('OK')
	
	def load_flat(self,MasterFlat):
		print('Loading MasterFlat ...'),
		try:
			MasterFlat_HDU    = pyfits.open(MasterDark)
			self.MasterFlat_Data   = MasterDark_HDU[0].data
			self.MasterFlat_Header = MasterDark_HDU[0].header
			self.MasterFlat_Texp   = float(ImageTest.correct_exposure(MasterFlat_Header))
		except: raise
		else: print('OK')
		
	def load_bias(self,MasterBias):
		print('Loading MasterBias ...'),
		try:
			MasterBias_HDU    = pyfits.open(MasterDark)
			self.MasterBias_Data   = MasterDark_HDU[0].data
			self.MasterBias_Header = MasterBias_HDU[0].header
			self.MasterBias_Texp   = float(ImageTest.correct_exposure(MasterBias_Header))
		except: raise
		else: print('OK')
			
	def reduce_science_frame(self,MasterDark,MasterFlat,MasterBias=None):
		'''
		Load MasterDark and MasterFlat. MasterBias is neccesary only if working 
		with different exposures between Dark and Science frames
		'''
		self.load_dark(MasterDark)
		self.load_flat(MasterFlat)
		
		if self.MasterDark_Texp == self.fits_Texp or MasterBias==None:
			if MasterBias==None:
				print('WARNING: Science and Dark dont have the same exposure ! ')
			self.SyntDark_Data   = self.MasterDark_Data
			self.SyntDark_Texp   = self.MasterDark_Texp
			self.SyntDark_Header = self.MasterDark_Header
		else:
			self.load_bias(MasterBias)		
			print('Creating synthetic Dark ...'),
			try:
				self.SyntDark_Data = (self.MasterDark_Data-self.MasterBias_Data)/ \
					(self.MasterDark_Texp-self.MasterBias_Data) *\
					(self.ScienceFrame_Texp-self.MasterBias_Texp)+\
					self.MasterBias_Data
				self.SyntDark_Texp   = self.fits_Texp
				self.SyntDark_Header = self.MasterDark_Header
				self.SyntDark_Header['EXPOSURE'] = self.SyntDark_Texp
			except: raise
			else: print('OK')
		
		print('Calibrating image with MasterFlat and MasterDark ...')
		try: self.fits_Data = (self.fits_Data-self.SyntDark_Data)/self.MasterFlat_Data
		except: raise
		else: print('OK')
		
	def __del__(self):
		del(self)

class ImageInfo(ImageTest,ConfigOptions):
	'''
	Extract some data from the image header
	and put it in a class ImageInfo
	'''
	def __init__(self,fits_header,config_file):
		ConfigOptions.__init__(self,config_file)
		
	def read_header(self,fits_header):
		# Date and time in different formats
		self.fits_date	 = ImageTest.correct_date(fits_header)
		self.date_array	 = [self.fits_date[0:4],self.fits_date[4:6],self.fits_date[6:8],
			self.fits_date[9:11],self.fits_date[11:13],self.fits_date[13:15]]
		self.date_string = date_array[0]+"/"+date_array[1]+"/"+date_array[2]+" "+\
			self.date_array[3]+":"+self.date_array[4]+":"+self.date_array[5]

		# Exposure (float), resolution (2d int array), filter (str)
		self.exposure	 = ImageTest.correct_exposure(fits_header)
		self.resolution	 = ImageTest.correct_resolution(fits_header)
		self.used_filter = ImageTest.correct_filter(fits_header)
	
	def config_processing(self):
		filters=["U","B","V","R","I"]
		
		# Default values
		self.zero_points = {"Johnson_"+filtro:[False,False] for filtro in filtros}
		self.color_terms = {"Johnson_"+filtro:[False,False] for filtro in filtros}
		self.background_levels = {"Johnson_"+filtro:[False,False] for filtro in filtros}
		self.flatfield = {"Johnson_"+filtro:False for filtro in filtros}
		self.darkframe = False
		self.biasframe = False
		
		# Config processing
		for option in self.FileOptions:
			if   option[0]=="obs_latitude":        self.latitude=float(option[1])
			elif option[0]=="obs_longitude":       self.longitude=float(option[1])
			elif option[0]=="delta_x":             self.despl_x_ptos=float(option[1])
			elif option[0]=="delta_y":             self.despl_y_ptos=float(option[1])
			elif option[0]=="radial_factor":       self.radial_factor=float(option[1])
			elif option[0]=="azimuth_zeropoint":   self.azimuth_zeropoint=float(option[1])
			elif option[0]=="min_altitude":        self.min_altitude=float(option[1])
			elif option[0]=="base_radius":         self.base_radius=float(option[1])
			elif option[0]=="baseflux_detectable": self.baseflux_detectable=float(option[1])
			elif option[0]=="lim_Kendall_tau":     self.lim_Kendall_tau=float(option[1])
			elif option[0]=="ccd_bits":            self.ccd_bits=float(option[1])
			elif option[0]=="ccd_gain":            self.ccd_gain=float(option[1])
			elif option[0]=="read_noise":          self.read_noise=float(option[1])
			elif option[0]=="thermal_noise":       self.thermal_noise=float(option[1])
			elif option[0]=="max_magnitude":       self.max_magnitude = float(option[1])
			elif option[0]=="max_star_number":     self.max_star_number = int(option[1])
			elif option[0]=="pixel_scale":         self.pixel_scale = float(option[1])
			elif option[0]=="backgroundmap_title": self.backgroundmap_title = str(option[1])
			elif option[0]=="darkframe":           self.darkframe=option[1]
			elif option[0]=="biasframe":           self.biasframe=option[1]
			else:
				# Options that depends on the filter
				for the_filter in filters:
					filter_name = "Johnson_"+filters[the_filter]
					if   option[0]=="zero_point_"+filters[the_filter]:
						self.zero_points[filter_name] = \
							[float(option[1].split(",")[0]), float(option[1].split(",")[1])]
					elif option[0]=="color_term_"+filters[the_filter]:
						self.color_terms[filter_name] = \
							[float(option[1].split(",")[0]), float(option[1].split(",")[1])]
					elif option[0]=="bkgnd_minmax_"+filters[the_filter]:
						self.background_levels[filter_name] = [float(option[1].split(",")[0]), \
						float(option[1].split(",")[1])]
					elif option[0]=="flatfield_"+filters[the_filter]:
						self.flatfield[filter_name] = option[1]	

	def __del__(self):
		del(self)
		
class ImageTest():
	'''Perform some test on the image header and extract information'''
	
	@staticmethod
	def correct_exposure(file_header):
		# Exposure
		try: texp = float(file_header['EXPOSURE'])
		except: raise
		else:
			assert texp>0., '0s exposure time detected.'
			return texp
	
	@staticmethod
	def correct_date(file_header):
		# Date and time
		try: date = file_header['DATE']
		except:	raise
		else:
			assert len(date)==6, 'Date format not YYYYMMDD'
			return date	
	
	@staticmethod
	def correct_resolution(file_header):
		# Resolution
		try: resolution = [int(file_header['NAXIS1']),int(file_header['NAXIS2'])]
		except: raise
		else:
			assert resolution[0]>0 and resolution[1]>0, 'Matrix not 2 dimensional'
			return resolution
	
	@staticmethod
	def correct_filter(file_header):
		# Test if there's a known filter
		try: used_filter = file_header['FILTER']
		except: raise
		else:
			# due to an inconsistent format in AstMon, 
			# we found 4 possible formats 'Jonhson_V','JohnsonV','Johnson_V','JonhsonV'
			used_filter = used_filter.replace('_','')
			assert used_filter[0:7] in ['Johnson','Jonhson'], 'Filter type not Johnson'
			assert used_filter[7:] in ['U','B','V','R','I'], 'Filter not U,B,V,R or I'
			return 'Johnson_'+used_filter[7:]
	
	
	
