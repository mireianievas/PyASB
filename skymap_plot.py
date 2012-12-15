#!/usr/bin/env python

'''
Star detection module

Auxiliary functions to plot the SkyMap
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
	from draw_functions import draw_polar_axes
except:
	print 'One or more modules missing: pyfits,HeaderTest'
	raise SystemExit


def draw_polar_axes(skyimage,ImageInfo):
	''' Draw polar axes on the image '''
	num_cir=9,num_rad=12
	div_360 = [ k for k in range(1,360+1,1) if 360%k == 0 ]
	phase = 360/([ k for k in div_360 if num_rad-k>=0][-1])
	# Minimal plotted altitude
	zenith_xy = zenith_position(ImageInfo)
	try:
		min_altitude = max(
			10*int(xy2horiz(0,ImageInfo.Properties.resolution[1]/2,ImageInfo)[1]/10 +1),
			10*int(xy2horiz(ImageInfo.Properties.resolution[0],\
				ImageInfo.resolution[1]/2,ImageInfo)[1]/10 +1),
			10*int(xy2horiz(ImageInfo.Properties.resolution[0]/2,0,ImageInfo)[1]/10 +1),
			10*int(xy2horiz(ImageInfo.Properties.resolution[0]/2,\
				ImageInfo.resolution[1],ImageInfo)[1]/10 +1) )
	except:
		raise
		min_altitude = 0
	
	x_minalt,y_minalt = horiz2xy(0,min_altitude,ImageInfo)
	max_radius = sqrt((x_minalt-ImageInfo.resolution[0]/2)**2+\
		y_minalt-ImageInfo.resolution[1]/2)**2)

	def is_cardinal_point(angle):
		''' this function draws a cardinal point letter instead of number
		    if its a known angle '''
		try: cardinal={0:"N",45:"NE",90:"E",135:"SE",180:"S",225:"SW",270:"W",315:"NW"}[angle]
		except: return "$\hbox{"+str(angle)+"}^\\circ $"
		else: return "$\hbox{"+cardinal+"}$"
	
	for azimuth in range(0,360,phase):
		dx = -max_radius*cos((azimuth-ImageInfo.Astrometry.azimuth_zeropoint)*pi/180)
		dy =  max_radius*sin((azimuth-ImageInfo.Astrometry.azimuth_zeropoint)*pi/180)
		xlabel = zenith_xy[0]+dx/2
		ylabel = zenith_xy[1]+dy/2
		skyimage.add_patch(Arrow(zenith_xy[0],zenith_xy[1],dx,dy,\
			ls='dashed',alpha=0.15,lw=1,label='_nolegend_'))
		imagen.annotate(is_cardinal_point(azimuth),xy=(xlabel,ylabel),xycoords='data',\
			xytext=(0,3),textcoords='offset points',fontsize='small',label='_nolegend_',\
			alpha=0.75)
	
	altitude_labels = range(90,min_altitude-5,5*int((min_altitude-90)/(5*num_cir)))
	if min_altitude not in altitude_labels:
		altitude_labels.append(min_altitude)

	for altitude in altitude_labels:
		x,y = horiz2xy(0,altitude,ImageInfo)
		radius = sqrt((x-zenith_xy[0])**2 + (y-zenith_xy[1])**2)
		skyimage.add_patch(Circle((zenith_xy[0],zenith_xy[1]),radius,facecolor='none',\
			ls='dashed',lw=1,fill=False, alpha=0.15,label='_nolegend_'))
		skyimage.annotate("$"+str(altitude)+"^\\circ$",xy=(x,y), xycoords='data',\
			label='_nolegend_',xytext=(0, 3), textcoords='offset points',\
			fontsize='small',alpha=0.75)
	
	return skyimage

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
