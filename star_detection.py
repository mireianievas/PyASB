#!/usr/bin/env python

'''
Star detection module

Perform star detection and measure star fluxes
____________________________

This module is part of the PyAstMonUCM project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyAstMonUCM project"
__credits__ = ["Miguel Nievas"]
__license__ = "GPL"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


try:
	import matplotlib.pyplot as mpl
	import matplotlib.colors as mpc
	import matplotlib.patches as mpp
	from skymap_plot import * 
except:
	print 'One or more modules missing: pyfits,HeaderTest'
	raise SystemExit


def photometry_bouguervar(Star):
	# Calculate parameters used in bouguer law fit
	Star.25logF      = 2.5*log10(Star.flux_star)
	Star.25logF_unc  = (2.5/log(10))*,Star.flux_star_unc/Star.flux_star
	Star.m25logF     = Star.FilterMag+Star.25logF
	Star.m25logF_unc = Star.25logF_unc
	return Star

def star_photometry(fits_data,StarCatalog,ImageInfo,ObsPyephem):
	'''
	This function measures fluxes for Star in the StarCatalog
	and give the neccesary measures to fit a Bouguer extinction law.
	'''
	
	if ImageInfo.Config.skymap_file!=False:
		skyfigure,skyimage = create_skymap(fits_data,StarCatalog,ImageInfo,ObsPyephem)

	StarsMeasured = []
	for star_index in range(StarCatalog):
		Star = StarCatalog[star_index]
		# Define aperture disk for photometry
		R1,R2,R3 = aperture_photometry_radius(Star,ImageInfo)
		pixels1,pixels2,pixels3 = apphot_pixels(Star,R1,R2,R3,ImageInfo)
		if pixel1==[] or pixel1==[] or pixel3=[]:
			# Continue with next star, this one is not in the FoV.
			StarCatalog.pop(star_index)
			continue
		
		# Coarse flux measure
		Star.flux_star,Star.flux_background = aproximate_fluxes(fits_data,\
			pixels1,pixels3,ImageInfo)
		# NOTE: flux_background is absolute, not normalized.
		
		# Blacklist star positions for future background SB measures
		for pixel in pixels1+pixels2:
			starmask_pixel.append(pixel)
			starmask_background.append(Star.flux_background)
		
		# Delete stars below detection limit or low at sky
		Star.max_magnitude = magnitude_limit(Star,ImageInfo)
		Star.min_flux      = flux_lowlimit(Star,ImageInfo)
		
		if Star.FilterMag < Star.max_magnitude and Star.flux_star > Star.min_flux\
		and Star.altit_appa > ImageInfo.Config.min_altitude:
			# Fine flux measure and star position (centroid estimation)
			try:
				Star.Xcoord,Star.Ycoord = estimate_centroid(fits_data,pixels1+pixels2,Star)
				pixels1,pixels2,pixels3 = apphot_pixels(Star,R1,R2,R3,ImageInfo)
				Star.flux_star,Star.flux_star_unc = \
					precise_star_fluxes(fits_data,pixels1,pixels2,pixels3,R1,R2,ImageInfo)
				Star=photometry_bouguervar(Star)
			except:
				''' Continue with next star, we cannot calculate centroid or measure 
				    fluxes for this one. '''
				raise
				continue

			for pixel in pixels1:
				starmask_pixel.append(pixel)
				starmask_background.append(Star.flux_background)

			StarsMeasured.append(Star)
			
			# If users wants the sky map with annotated stars
			if ImageInfo.Config.skymap_file!=False:
				skyimage = annotate_skymap(skyimage,Star,R1,R2,R3)
	
	if ImageInfo.Config.skymap_file!=False:
		# Show or save the map
		show_or_save_skymap(skyfigure,ImageInfo,ObsPyephem)
			
	return StarsMeasured

		


