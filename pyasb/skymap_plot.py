#!/usr/bin/env python

'''
SkyMap module

Auxiliary functions to plot the SkyMap
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
	from astrometry import *
	from scipy.ndimage import uniform_filter as denoise
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	import matplotlib.colors as mpc
	import matplotlib.patches as mpp
	import matplotlib as mpl
except:
	print(str(inspect.stack()[0][2:4][::-1])+': One or more modules missing')
	raise SystemExit

class SkyMap():
	'''	SkyMap class '''
	
	def __init__(self,StarCatalog,ImageInfo,FitsImage):
		self.StarCatalog = StarCatalog
                self.ImageInfo   = ImageInfo

		if self.ImageInfo.skymap_path==False:
			# Don't draw anything
			print('Skipping Skymap Graph ...')
			return(None)
		else:
			print('Star Map plot ...')
		
		stretched_fits_data = self.stretch_logdata(FitsImage.fits_data_notcalibrated,40,99)
		self.define_skymap();
		self.draw_skymap_data(stretched_fits_data)
		self.draw_polar_axes()
		for Star in self.StarCatalog.StarList_Det:
			self.annotate_skymap(Star)
		
                self.astrometry_solver()
		self.show_figure()
	
	def stretch_logdata(self,fits_data,pmin,pmax):
		log_fits_data = np.log(fits_data-np.min(fits_data)+1,dtype="float32")
		valuemin = np.percentile(log_fits_data,pmin)
		valuemax = np.percentile(log_fits_data,pmax)
		stretched_fits_data = log_fits_data.clip(valuemin,valuemax)
		return stretched_fits_data
	
	def define_skymap(self):
		''' Create figure and self.skyimage subplot. '''
		self.skyfigure = plt.figure(figsize=(10,10))
		self.skyimage  = self.skyfigure.add_subplot(111)
	
	def mouse_press_callback(self,event):
		''' Coordinate input '''
		if event.button == 3:
			try:
				self.name
				self.azim
				self.alti
			except:
				self.identified_stars = []
				self.star_index = 0
				print('Please, right click on the following stars')
			else:
				ix, iy = event.xdata, event.ydata
				print('x = %d, y = %d'%(ix, iy))
				self.identified_stars.append([self.name,self.azim,self.alti,ix,iy])
				self.star_index +=1

			self.name=self.StarCatalog.StarList_Tot[self.star_index].name
			self.azim=self.StarCatalog.StarList_Tot[self.star_index].azimuth
			self.alti=self.StarCatalog.StarList_Tot[self.star_index].altit_real
			print('Name: %s, Az: %s, Alt: %s' %(self.name,self.azim,self.alti))
		
		return(None)
	
	def key_press_callback(self, event):
	        'whenever a key is pressed'
        	if not event.inaxes: return(None)
		
		if event.key=='n':
			print('Next star')
			self.star_index += 1
			self.name=self.StarCatalog.StarList_Tot[self.star_index].name
			self.azim=self.StarCatalog.StarList_Tot[self.star_index].azimuth
			self.alti=self.StarCatalog.StarList_Tot[self.star_index].altit_real
			print('Name: %s, Az: %s, Alt: %s' %(self.name,self.azim,self.alti))

		if event.key=='p':
			print('Previous star')
			self.identified_stars.pop()
			self.star_index -= 1
			
			self.name=self.StarCatalog.StarList_Tot[self.star_index].name
			self.azim=self.StarCatalog.StarList_Tot[self.star_index].azimuth
			self.alti=self.StarCatalog.StarList_Tot[self.star_index].altit_real
			print('Name: %s, Az: %s, Alt: %s' %(self.name,self.azim,self.alti))

	        if event.key=='q':
			print('End')
			self.skyfigure.canvas.mpl_disconnect(self.cid_mouse)
			self.skyfigure.canvas.mpl_disconnect(self.cid_keyboard)
                        print(self.identified_stars)
                        self.astrometry_optimizer()

		return(None)

        def astrometry_optimizer(self):
                from scipy.optimize import minimize
                from astrometry import horiz2xy

                ImageInfo = self.ImageInfo

                def horiz2xy_chi2(sol,az,alt,x,y):
                    ImageInfo.radial_factor     = sol[0]
                    ImageInfo.azimuth_zeropoint = sol[1]
                    ImageInfo.delta_x           = sol[2]
                    ImageInfo.delta_y           = sol[3]
                    ImageInfo.latitude_offset   = sol[4]
                    ImageInfo.longitude_offset  = sol[5]
                    xf,yf = horiz2xy(az,alt,ImageInfo,derotate=True)
                    return(np.sum((xf-x)**2 + (yf-y)**2))

                coords = np.array(self.identified_stars)[:,1:] # Remove star name
                coords = np.array(coords,dtype=float)          # Convert to float
                [_az,_alt,_x,_y] = np.transpose(coords)        # Transpose and split
                print('Solving equation system')
                res = minimize(horiz2xy_chi2,[10,0,0,0,0,0],args = (_az,_alt,_x,_y))
                print(res.x)
                print(res)

	def astrometry_solver(self):
		print(\
			'*** Star select tool. Press right-click to begin. *** \n'+\
			'Right-click: assign star coords. \n'+\
			'n:           next star (skip current). \n'+\
			'p:           previous star (remove last entry). \n'+\
			'q:           quit star select tool. \n')
		
		self.completed=0
		self.cid_mouse = self.skyfigure.canvas.mpl_connect('button_press_event', self.mouse_press_callback)
		self.cid_keyboard = self.skyfigure.canvas.mpl_connect('key_press_event', self.key_press_callback)

	
	def draw_skymap_data(self,fits_data):
		''' Draw identified stars with image as background'''
		xpoints_catalog = [Star.Xcoord for Star in self.StarCatalog.StarList_Tot]
		ypoints_catalog = [Star.Ycoord for Star in self.StarCatalog.StarList_Tot]
		
		xpoints_phot = [Star.Xcoord for Star in self.StarCatalog.StarList_Det]
		ypoints_phot = [Star.Ycoord for Star in self.StarCatalog.StarList_Det]
		
		self.skyimage.imshow(denoise(fits_data, 5),cmap=mpl.cm.gray)
		
		self.skyimage.scatter(xpoints_catalog,ypoints_catalog,\
			marker='+',s=20,c='yellow',alpha=0.75,label='Star in the catalog')
		
		self.skyimage.scatter(xpoints_phot,ypoints_phot,\
			marker='+',s=20,c='r',alpha=0.75,label='Stars found')
		
		self.skyimage.axis([0,self.ImageInfo.resolution[0],0,self.ImageInfo.resolution[1]])
		information=str(self.ImageInfo.date_string)+" UTC\n"+str(self.ImageInfo.latitude)+5*" "+\
			str(self.ImageInfo.longitude)+"\n"+self.ImageInfo.used_filter
		self.skyimage.text(0.005,0.005,information,fontsize='small',color='white',\
			transform = self.skyimage.transAxes,backgroundcolor=(0,0,0,0.75))
		self.skyimage.legend(('In catalog','Detected'),'upper right')
	
	def draw_polar_axes(self):
		''' Draws meridian and altitude isolines. '''
		
		zenith_xy = zenith_position(self.ImageInfo)
		
		for each_altitude in np.arange(0,90,15):
			coord_altitude_0 = horiz2xy(0,each_altitude,self.ImageInfo)
			radius = math.sqrt(\
				(coord_altitude_0[0]-zenith_xy[0])**2 +\
				(coord_altitude_0[1]-zenith_xy[1])**2)
			self.skyimage.add_patch(\
				mpp.Circle((zenith_xy[0],zenith_xy[1]),radius,
				facecolor='k',fill=False, alpha=0.2,label='_nolegend_'))
			self.skyimage.annotate(\
				str(each_altitude),
				xy=(radius+zenith_xy[0],zenith_xy[1]),
				alpha=0.2,
				fontsize=10)
		
		key_azimuths = {0: "N",90: "E", 180: "S", 270: "W"}
		
		for each_azimuth in np.arange(0,360,30):
			coord_azimuth_0 = horiz2xy(each_azimuth,0,self.ImageInfo)
			self.skyimage.plot(\
				[zenith_xy[0],coord_azimuth_0[0]],
				[zenith_xy[1],coord_azimuth_0[1]],
				color='k',
				alpha=0.2,)
			
			if each_azimuth in key_azimuths:
				azimuth_label = str(key_azimuths[each_azimuth])
			else:
				azimuth_label = str(each_azimuth)
			self.skyimage.annotate(\
				azimuth_label,
				xy=horiz2xy(each_azimuth,self.ImageInfo.min_altitude,self.ImageInfo),
				color='k',
				alpha=0.2,
				fontsize=10)
		
	def annotate_skymap(self,Star):
		# Draw identified stars and measuring circ<les.
		# Annotate HD catalog code and Magnitude for each star.
		self.skyimage.scatter(Star.Xcoord,Star.Ycoord,marker='x',c='r',alpha=0.2,label='Identified stars')
		self.skyimage.add_patch(mpp.Circle((Star.Xcoord,Star.Ycoord),Star.R1,facecolor='none',edgecolor=(0,0,0.8),\
			linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
		self.skyimage.add_patch(mpp.Circle((Star.Xcoord,Star.Ycoord),Star.R2,facecolor='none',edgecolor=(0,0.8,0),\
			linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
		self.skyimage.add_patch(mpp.Circle((Star.Xcoord,Star.Ycoord),Star.R3,facecolor='none',edgecolor=(0.8,0,0),\
			linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
		self.skyimage.annotate(Star.name,xy=(Star.Xcoord,Star.Ycoord), xycoords='data',xytext=(0, 3),\
			textcoords='offset points',fontsize=8)
		self.skyimage.annotate(Star.FilterMag,xy=(Star.Xcoord,Star.Ycoord), xycoords='data',xytext=(0,-10),\
			textcoords='offset points',fontsize=8)
			
	def show_figure(self):
		def skymap_filename():
			filename = self.ImageInfo.skymap_path +\
				"/Skymap_"+self.ImageInfo.obs_name+"_"+self.ImageInfo.fits_date+"_"+\
				self.ImageInfo.used_filter+".png"
			return(filename)
		
		if self.ImageInfo.skymap_path=="screen":
			plt.show()
		else:
			plt.savefig(skymap_filename(),bbox_inches='tight')
		
		#plt.clf()
		#plt.close('all')
