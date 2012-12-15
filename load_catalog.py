#!/usr/bin/env python

'''
Load Catalog file and make PyAstMon StarCatalog

This module loads the catalog file and returns 
an array of Star objects (StarCatalog).
____________________________

This module is part of the PyAstMonUCM project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

try:
	import ephem
	import CommonErrors, HeaderTest
	from astrometry import atmospheric_refraction,calculate_airmass
except:
	print 'One or more modules missing: pyfits,CommonErrors,HeaderTest'
	raise SystemExit

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyAstMonUCM project"
__credits__ = ["Miguel Nievas"]
__license__ = "GPL"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


def open_catalog_file(catalog_filename='ducati_catalog.tsv'):
	# Open catalog file and return its content
	try: 
		tabseparated_catalog = [ line[:-1] for line in \
			open(catalog_filename, 'r').readlines() ]
	except IOError:
		print 'IOError. Error opening file '+catalog_filename+'.'
		return 1
	except:
		print 'Unknown error:'
		raise
		return 2
	else:
		print 'File '+str(catalog_filename)+' opened correctly.'
		return tabseparated_catalog

class StarObject:
	''' # Object Star. Original catalog values :
	    # recno, HDcode, RA1950, DEC1950, Vmag, U_V, B_V, R_V, I_V '''
	pass

def coord_pyephem_format(coord_str):
	# Return coordinate in pyepheem str
	coord_separated = coord_str.split()
	coord_pyephem = str(int(coord_separated[0]))+\
		':'+str(int(coord_separated[1]))+str(int(coord_separated[2]))
	return coord_pyephem

def build_pyastmon_catalog(tabseparated_catalog):
	'''
	Read catalog lines, extract Star properties and put them
	in Object array. 
	Properties are image independent
	'''
	
	print 'Reading catalog ...'
	# Catalog starts when -------- detected at the begining of the line
	startline = 0
	while tabseparated_catalog[startline][0:8] != 8*'-':
		startline+=1
	
	# Build catalog matrix without header
	tabseparated_catalog = [ tabseparated_catalog[line] for line in \
		range(len(tabseparated_catalog)) if line>startline ]

	# Build StarCatalog (array of StarObject)
	StarCatalog = [ StarObject for line in tabseparated_catalog]

	for line in range(len(tabseparated_catalog)-1,0-1,-1):
		try:
			tabseparated_catalog_split = tabseparated_catalog.split('\t')
			StarCatalog[line].recno   = int(tabseparated_catalog_split[0])
			StarCatalog[line].HDcode  = str(tabseparated_catalog_split[1]).replace(' ','')
			StarCatalog[line].RA1950  = coord_pyephem_format(tabseparated_catalog_split[2])
			StarCatalog[line].DEC1950 = coord_pyephem_format(tabseparated_catalog_split[3])
			StarCatalog[line].Vmag    = float(tabseparated_catalog_split[4])
			StarCatalog[line].U_V     = float(tabseparated_catalog_split[5])
			StarCatalog[line].B_V     = float(tabseparated_catalog_split[6])
			StarCatalog[line].R_V     = float(tabseparated_catalog_split[7])
			StarCatalog[line].I_V     = float(tabseparated_catalog_split[8])

			try: 
				StarCatalog[line].name = str(tabseparated_catalog_split[9])
			except:
				StarCatalog[line].name = str(tabseparated_catalog_split[1])

			StarCatalog[line].Umag = StarCatalog[line].Vmag + StarCatalog[line].U_V
			StarCatalog[line].Bmag = StarCatalog[line].Vmag + StarCatalog[line].B_V 
			StarCatalog[line].Rmag = StarCatalog[line].Vmag + StarCatalog[line].R_V
			StarCatalog[line].Imag = StarCatalog[line].Vmag + StarCatalog[line].I_V
		except:
			print 'Error loading star '+line+', deleting from catalog... '
			raise
			StarCatalog.pop(line)
	
	return StarCatalog

def imagedep_properties_catalog(StarCatalog,ObsPyephem,ImageInfo):
	'''
	Calculate star position and dynamic properties.
	Don't proceed if the star altitude < 0.
	Return updated StarCatalog (image dependent)
	'''

	def pyephem_declaration(Star):
		# Define the star in Pyephem to make astrometric calculations
		pyephem_star = ephem.readdb('"'+Star.name+'"'+",f|S|A0,"+Star.RA1950+'|0'+\
			","+Star.DEC1950+'|0'+","+Star.Vmag+',1950,0"')
		pyephem_star.compute(ObsPyephem)	
		return pyephem_star

	def actual_filter(Star,ImageInfo):
		# Test which filter is in use
		used_filter = ImageInfo.Properties.used_filter
		class StarFilter:
			if used_filter=="JohnsonU":
				Mag   = Star.Umag
				Color = Star.U_V
			elif used_filter=="JohnsonB":
				Mag   = Star.Umag
				Color = Star.U_V
			elif used_filter=="JohnsonV":
				Mag   = Star.Vmag
				Color = 0.0
			elif used_filter=="JohnsonR":
				Mag   = Star.Rmag
				Color = Star.R_V
			elif used_filter=="JohnsonI":
				Mag   = Star.Imag
				Color = Star.I_V
			else:
				pass
		return StarFilter
	
	for line in range(len(StarCatalog)-1,0-1,-1):
		try:
			pyephem_star = PyephemDeclaration(StarCatalog[line])
			if float(pyephem_star.alt) > 0.0:
				# Real coordinates (from catalog)
				StarCatalog[line].altit_real = float(pyephem_star.alt)
				StarCatalog[line].zdist_real = 90.0-StarCatalog[line].altit_real
				StarCatalog[line].azimuth    = float(pyephem_star.az)
				# Apparent coordinates in sky. Atmospheric refraction effect.
				StarCatalog[line].altit_appa = atmospheric_refraction(StarCatalog[line].altit_real,'dir')
				StarCatalog[line].zdist_appa = 90.0-StarCatalog[line].altit_apar
				StarCatalog[line].airmass    = calculate_airmass(StarCatalog[line].altit_appa)
				# Photometric properties
				Starfilter = actual_filter(StarCatalog[line],ImageInfo)
				StarCatalog[line].FilterMag = Starfilter.Mag
				StarCatalog[line].Color     = Starfilter.Color
				# Image coordinates
				XYCoordinates = horiz2xy(StarCatalog[line].azimuth,\
					StarCatalog[line].altit_appa,ImageInfo)
				StarCatalog[line].Xcoord = XYCoordinates[0]
				StarCatalog[line].Ycoord = XYCoordinates[1]
				if StarCatalog[line].Xcoord<0. or StarCatalog[line].Ycoord<0.\
				or StarCatalog[line].Xcoord>ImageInfo.Properties.resolution[0] \
				or StarCatalog[line].Ycoord>ImageInfo.Properties.resolution[1]:
					# Star doesn't fit in the image
					StarCatalog.pop(line)
			else:
				StarCatalog.pop(line)
		except:
			StarCatalog.pop(line)
	

	return StarCatalog

