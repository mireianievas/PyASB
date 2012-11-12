#!/usr/bin/env python

'''
Star detection module

Perform star detection and measure star fluxes
____________________________

This module is part of the PyAstMonUCM project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

try:
	import matplotlib.pyplot as mpl
	import matplotlib.colors as mpc
	import matplotlib.patches as mpp
	from draw_functions import draw_polar_axes
except:
	print 'One or more modules missing: pyfits,HeaderTest'
	raise SystemExit




def create_skymap(fits_data,StarCatalog,ImageInfo,ObsPyephem):
	# Create figure and skyimage subplot. 
	# Set axis, labels, info text and draw stars in catalog.
	# Return th
	skyfigure = mpl.figure(figsize=(10,10),dpi=100)
	skyimage  = skyfigure.add_subplot(111)

	skyimage.set_title('Stars in the catalog and identified stars',size="xx-large")
	skyimage.imshow(fits_data,norm=mpc.LogNorm(),cmap=cm.gray)

	xpoints = [Star.Xcoord for Star in StarCatalog]
	ypoints = [Star.Ycoord for Star in StarCatalog]

	skyimage.scatter(xpoints,ypoints,marker='.',c='g',alpha=0.2,label='Star in the catalog')
	skyimage.axis([0,ImageInfo.Properties.resolution[0],0,	
		ImageInfo.Properties.resolution[1]])
	information_skyimage=str(ObsPyephem.date)+" UTC\n"+str(ObsPyephem.lat)+5*" "+\
		str(ObsPyephem.lon)+"\n"+ImageInfo.Properties.used_filter
	skyimage.text(0.005,0.005,information_skyimage,fontsize='x-small',color='white',\
		transform = skyimage.transAxes,backgroundcolor=(0,0,0,0.75))
	skyimage = draw_polar_axes(skyimage,ImageInfo)
	skyimage.legend(('In catalog','Detected'),'upper right')
	return skyfigure,skyimage

def annotate_skymap(skyimage,Star,R1,R2,R3):
	# Draw identified stars and measuring circles.
	# Annotate HD catalog code and Magnitude for each star.
	skyimage.scatter(Star.Xcoord,Star.Ycoord,marker='x',c='r',alpha=0.2,label='Identified stars')
	skyimage.add_patch(mpp.Circle((Star.Xcoord,Star.Ycoord),R1,facecolor='none',edgecolor=(0,0,0.8),\
		linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
	skyimage.add_patch(mpp.Circle((Star.Xcoord,Star.Ycoord),R2,facecolor='none',edgecolor=(0,0.8,0),\
		linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
	skyimage.add_patch(mpp.Circle((Star.Xcoord,Star.Ycoord),R3,facecolor='none',edgecolor=(0.8,0,0),\
		linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
	skyimage.annotate(Star.HDcode,xy=(Star.Xcoord,Star.Ycoord), xycoords='data',xytext=(0, 3),\
		textcoords='offset points',fontsize=6)
	skyimage.annotate(Star.FilterMag,xy=(Estrella.X,Estrella.Y), xycoords='data',xytext=(0,-10),\
		textcoords='offset points',fontsize=6)
	return skyimage

def show_or_save_skymap(skyfigure,ImageInfo,ObsPyephem):
	# Show or save the skymap
	if ImageInfo.Config.skymap_file=='on_screen':
		mpl.show(skyfigure)
	else:
		mpl.savefig(complete_file_name(ImageInfo))

def star_photometry(fits_data,StarCatalog,ImageInfo,ObsPyephem):
	'''
	This function measures fluxes for Star in the StarCatalog
	and give the neccesary measures to fit a Bouguer extinction law.
	'''
	
	if ImageInfo.Config.skymap_file!=False:
		skyfigure,skyimage = create_skymap(fits_data,StarCatalog,ImageInfo,ObsPyephem)

	StarMeasured = []
	for star_index in range(StarCatalog):
		Star = StarCatalog[star_index]
		# Define aperture disk for photometry
		R1,R2,R3 = aperture_photometry_radius(Star,ImageInfo)
		pixels1,pixels2,pixels3 = measure_pixels(Star,R1,R2,R3,ImageInfo)
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
			except:
				# Continue with next star, we cannot calculate centroid for this one.
				raise
				continue

			# Saturation photometry
			try:
				pixels1,pixels2,pixels3 = measure_pixels(Star,R1,R2,R3,ImageInfo)
				Star.flux_star,Star.flux_star_uncertainty = \
					precise_star_fluxes(fits_data,pixels1,pixels2,pixels3,R1,R2,ImageInfo)
			except:
				# Continue with next star, we cannot measure fluxes on this.
				raise
				continue

			for pixel in pixels1:
				starmask_pixel.append(pixel)
				starmask_background.append(Star.flux_background)

			StarMeasured.append(Star)
			
			# If users wants the sky map with annotated stars
			if ImageInfo.Config.skymap_file!=False:
				skyimage = annotate_skymap(skyimage,Star,R1,R2,R3)
	
	if ImageInfo.Config.skymap_file!=False:
		# Show or save the map
		show_or_save_skymap(skyfigure,ImageInfo,ObsPyephem)
			


		


