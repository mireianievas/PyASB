
'''
PyASB input options module

Read input options
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

class ReadOptions(object):
	def __init__(self,input_options):
		# Lambda: http://docs.python.org/release/2.5.2/tut/node6.html [4.7.5]
		self.options = { '-h': lambda: 'show_help', '-d': lambda : 'use_date', '-i': lambda : 'use_file', 
			'-om': lambda : 'use_skymap', '-or': lambda : 'use_results', '-ot': lambda : 'use_phottable', 
			'-ob': lambda: 'use_bouguerfit', '-os': lambda: 'use_sb' , '-ost': lambda: 'use_sbtable'}
		
		# By default, we wont show on screen nor save on disk.
		self.show_help = False;
		self.photometry_table_path=False; 
		self.skymap_table_path=False; 
		self.bouguerfit_path=False;
		self.skybrightness_map_path=False; 
		self.skybrightness_table_path=False; 
		self.summary_path=False; 
		
		print('Input Options: '+str(input_options))
		
		self.input_options = input_options
		try: self.input_options[1]
		except: self.no_parameters()
		else:
			while len(self.input_options)>1:
				input_option = self.options.get(self.input_options[1], lambda : None)()
				if input_option == 'show_help':
					self.show_help = True
					# Stop reading options. Program will halt
					self.input_options = []
				elif input_option == 'use_date':
					self.dates=self.input_date()
				elif input_option == 'use_file':
					self.fits_filename_list = self.input_file()
				elif input_option == 'use_phottable':
					self.photometry_table_path=self.output_file()
				elif input_option == 'use_skymap':
					self.skymap_path=self.output_file()
				elif input_option == 'use_bouguerfit':
					self.bouguerfit_path=self.output_file()
				elif input_option == 'use_sb':
					self.skybrightness_map_path=self.output_file()
				elif input_option == 'use_sbtable':
					self.skybrightness_table_path=self.output_file()
				elif input_option == 'use_results':
					self.summary_path=self.output_file()
				else:
					self.input_options.remove(self.input_options[1])
					continue
					#self.incorrect_parameter()
		
		if self.show_help==False:
			self.date_set=True 
			self.inputfile_set=True
			
			try: 
				assert(len(dates)>=1)
			except: self.date_set=False
			
			try: 
				assert(len(self.fits_filename_list)>=1)
			except: self.inputfile_set=False
			
			if not (self.date_set or self.inputfile_set):
				self.need_date_or_file()
	
	def incorrect_parameter(self):
		print '\nERR: Incorrect parameter '+str(self.input_options[1])
		self.input_options = []
		self.show_help = True
	
	def need_date_or_file(self):
		print '\nERR: Need date or input file to proceed'
		self.input_options = []
		self.show_help = True
	
	def no_parameters(self):
		print '\nERR: Need more than one parameter'
		self.input_options = []
		self.show_help = True
	
	def input_file(self):
		# Input could be a list of files, comma separated
		try: self.input_options[2]
		except: self.need_date_or_file()
		else:
			if self.options.get(self.input_options[2], lambda : None)():
				self.need_date_or_file()
			else:
				iterate = True
				list_files = []
				while(iterate == True):
					list_files.append(self.input_options[2])
					self.input_options.remove(self.input_options[2]); 
					if self.input_options[2].replace(" ","")[0]=="-":
						iterate=False
				
				self.input_options.remove(self.input_options[1])
				return list_files
	
	def output_file(self):
		# If output is not disabled, then show on screen or save to file
		try: self.input_options[2]
		except:
			self.input_options.remove(self.input_options[1])
			return 'screen'
		else:
			if self.options.get(self.input_options[2], lambda : None)():
				self.input_options.remove(self.input_options[1])
				file_output='screen'
			else:
				file_output=self.input_options[2]
				self.input_options.remove(self.input_options[2]); self.input_options.remove(self.input_options[1])
			return(file_output)
	
	def input_date(self):
		# Format input date (takes into account short format)
		# Date should be something like 
		# YYYY-MM-DD (long format, complete), YY-M-D (short format, complete), 
		# YYYY-MM (lf, entire month), YY (sf, entire year) ...
		try: self.input_options[2]
		except: self.need_date_or_file()
		else:
			try: date=self.input_options[2].split("-")
			except: 
				self.need_date_or_file()
			else:
				if self.input_options[2]=='todo' or self.input_options[2]=="todo":
					years=range(2010,2012+1,1)
					months=range(1,12+1,1)
				else:
					if len(date)==1:
						years=[int(date[0])]
						months=range(1,12+1,1)
					elif len(date)==2:
						years=[int(date[0])]
						months=[int(date[1])]
					elif len(date)==3:
						years=[int(date[0])]
						months=[int(date[1])]
						days=[int(date[2])]
					else:
						self.need_date_or_file()
			days_month={1:31, 2:28, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30, 10:31, 11:30, 12:31}
			dates=[]
			today_yearmonthday=str(datetime.now()).split(' ')[0].split('-')
			for year in years:
				# Short format processing
				year_ = year%2000 + 2000
				for month in months:
					try: days
					except: 
						# Leap year?
						if month==2 and (year_%4==0 and ((year_ %100!=0) or (year_%400==0))):
							days=range(1,29+1,1)
						else:
							days=range(1,days_month[month]+1,1)
		
					for day in days:
						use_date=True
						if year_>int(today_yearmonthday[0]):
							use_date=False
						elif year_==int(today_yearmonthday[0]):
							if month>int(today_yearmonthday[1]):
								use_date=False
							elif month==int(today_yearmonthday[1]):
								if day>=int(today_yearmonthday[2]):
									use_date=False
						if use_date==True:
							dates.append([year_,month,day])
			self.input_options.remove(self.input_options[2]); self.input_options.remove(self.input_options[1])
			return dates
	
	

