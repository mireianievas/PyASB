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

def stellar_photometry(fits_data,StarCatalog,ImageInfo,ObsPyephem):
	'''
	This function measures fluxes for Star in the StarCatalog
	and give the neccesary measures to fit a Bouguer extinction law.
	'''
	filtered_fits_data = deepcopy(fits_data)

	if ImageInfo.Config.skymap_file!=False:
		skyfigure,skyimage = create_skymap(fits_data,StarCatalog,ImageInfo,ObsPyephem)

	StarsMeasured = []
	for star_index in range(StarCatalog):
		Star = StarCatalog[star_index]
		# Define aperture disk for photometry
		R1,R2,R3 = aperture_photometry_radius(Star,ImageInfo)
		pixels1,pixels2,pixels3 = apphot_pixels(Star.Xcoord,Star.Ycoord,R1,R2,R3,ImageInfo)
		if pixel1==[] or pixel1==[] or pixel3=[]:
			# Continue with next star, this one is not in the FoV.
			StarCatalog.pop(star_index)
			continue
		
		# Coarse flux measure
		Star.flux_star,Star.flux_sky = aproximate_fluxes(fits_data,pixels1,pixels3,ImageInfo)
		# NOTE: flux_sky is absolute, not normalized.
		
		# Filter the fits matrix for future background SB measures
		for pixel in pixels1+pixels2:
			filtered_fits_data[pixel[0]][pixel[1]] = Star.flux_sky
		
		# Delete stars below detection limit or low at sky
		Star.max_magnitude = magnitude_limit(Star,ImageInfo)
		Star.min_flux      = flux_lowlimit(Star,ImageInfo)
		
		if Star.FilterMag < Star.max_magnitude and Star.flux_star > Star.min_flux\
		and Star.altit_appa > ImageInfo.Config.min_altitude:
			try:
				# Fine flux measure and star position (centroid estimation)
				Star.Xcoord,Star.Ycoord = estimate_centroid(fits_data,pixels1+pixels2,Star)
				pixels1,pixels2,pixels3 = \
					apphot_pixels(Star.Xcoord,Star.Ycoord,R1,R2,R3,ImageInfo)
				Star.flux_sky,Star.flux_sky_err   = sky_flux_measure(fits_data,pixels3,ImageInfo)
				Star.flux_star,Star.flux_star_err = \
					star_flux_measure(fits_data,pixels1,Star.flux_sky,Star.flux_sky_err,ImageInfo)
				# Log of flux and m+2.5logF, needed for bouguer law fit
				Star=photometry_bouguervar(Star)
				# Filter the fits matrix for future background SB measures
				for pixel in pixels1+pixels2:
					filtered_fits_data[pixel[0]][pixel[1]] = Star.flux_sky
				# Append to the measured star list and plot it if neccesary
				StarsMeasured.append(Star)
				if ImageInfo.Config.skymap_file!=False:
					skyimage = annotate_skymap(skyimage,Star,R1,R2,R3)

			except:
				''' Continue with next star, we cannot calculate centroid or measure 
				    fluxes for this one. '''
				print 'Cannot measure fluxes on '+str(Star.HDcode)+'. Trying next star...'
				continue

	if ImageInfo.Config.skymap_file!=False:
		# Show or save the map
		show_or_save_skymap(skyfigure,ImageInfo,ObsPyephem)
	
	return StarsMeasured,filtered_fits_data

		


