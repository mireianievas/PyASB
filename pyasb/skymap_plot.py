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
	#from scipy.ndimage import uniform_filter as denoise
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
	
	def __init__(self,ImageInfo,FitsImage):
		# Set ImageInfo as local sub-object, we will use
		# it a lot.
		self.ImageInfo = ImageInfo

		if self.ImageInfo.skymap_path==False:
			# Don't draw anything
			print('Skipping Skymap Graph ...')
		else:
			print('Star Map plot ...')
			self.stretch_data(\
				FitsImage.fits_data,\
				ImageInfo.perc_low,\
				ImageInfo.perc_high)
			#self.setup_skymap()
			#self.draw_catalog_stars()
			#self.draw_detected_stars()
			#self.astrometry_solver()
			#self.draw_polar_axes()
			#self.show_figure()
	
	def setup_skymap(self):
		'''
		To be executed at the beginning (no stars)
		'''
		if (self.ImageInfo.skymap_path!=False):
			self.define_skymap()
			self.draw_skymap_data()
			self.skyfigure.canvas.draw()
			self.skyfigure.canvas.flush_events()
			plt.show(block=False)
	
	def complete_skymap(self):
		''' 
		To be executed when an astrometric solution is found
		'''
		if (self.ImageInfo.skymap_path!=False):
			self.draw_catalog_stars()
			self.draw_detected_stars()
			self.draw_polar_axes()
			self.skyfigure.canvas.draw()
			self.skyfigure.canvas.flush_events()
			self.show_figure()
	
	def set_starcatalog(self,StarCatalog):
		self.StarCatalog = StarCatalog
	
	def draw_catalog_stars(self):
		for Star in self.StarCatalog.StarList_Tot:
			self.draw_annotate_star(Star, full=False)
	
	def draw_detected_stars(self):
		for Star in self.StarCatalog.StarList_Det:
			self.draw_annotate_star(Star, full=True)
	
	def stretch_data(self,fits_data,pmin,pmax):
		#log_fits_data = np.log(fits_data-np.min(fits_data)+1,dtype="float32")
		log_fits_data = np.arcsinh(fits_data-np.min(fits_data)+1,dtype="float32")
		valuemin = np.percentile(log_fits_data,pmin)
		valuemax = np.percentile(log_fits_data,pmax)
		self.stretched_fits_data = log_fits_data.clip(valuemin,valuemax)
	
	def define_skymap(self):
		''' Create figure and self.skyimage subplot. '''
		self.skyfigure = plt.figure(figsize=(10,10))
		self.skyimage  = self.skyfigure.add_subplot(111)
		self.skyfigure.canvas.draw()#(block=False)
	
	def mouse_press_callback(self,event):
		''' Coordinate input '''
		if event.button == 3:
			ix, iy = event.xdata, event.ydata
			print('x = %d, y = %d'%(ix, iy))
			self.identified_stars.append([self.name,self.azim,self.alti,ix,iy])
			self.star_index +=1
			self.astrometry_optimizer(full=(self.star_index>3))
			self.scatter_stars.append(\
				self.skyimage.scatter(ix,iy,marker='o',c='red',alpha=0.2))
			self.label_stars.append(\
				self.skyimage.annotate(\
					self.name,xy=(ix,iy), \
					xycoords='data',xytext=(0, 3),\
					textcoords='offset points',fontsize=8,alpha=0.8))
			
			self.name=self.StarCatalog.StarList_Tot[self.star_index].name
			self.azim=self.StarCatalog.StarList_Tot[self.star_index].azimuth
			self.alti=self.StarCatalog.StarList_Tot[self.star_index].altit_real
			px,py = horiz2xy(self.azim,self.alti,self.ImageInfo,derotate=True)
			
			try: self.preliminary_star.remove()
			except: pass
			
			self.preliminary_star = \
				self.skyimage.scatter(px,py,marker='o',c='yellow',alpha=0.5)
			print('Name: %s, Az: %s, Alt: %s' %(self.name,self.azim,self.alti))
			self.skyfigure.canvas.draw()
			self.skyfigure.canvas.flush_events()
		
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
			px,py = horiz2xy(self.azim,self.alti,self.ImageInfo,derotate=True)
			self.preliminary_star.remove()
			self.preliminary_star = \
				self.skyimage.scatter(px,py,marker='o',c='yellow',alpha=0.5)
			print('Name: %s, Az: %s, Alt: %s' %(self.name,self.azim,self.alti))
		elif event.key=='p':
			print('Previous star')
			self.preliminary_star.remove()
			self.scatter_stars[-1].remove()
			self.label_stars[-1].remove()
			self.scatter_stars.pop()
			self.label_stars.pop()
			self.identified_stars.pop()
			self.star_index -= 1
			self.name=self.StarCatalog.StarList_Tot[self.star_index].name
			self.azim=self.StarCatalog.StarList_Tot[self.star_index].azimuth
			self.alti=self.StarCatalog.StarList_Tot[self.star_index].altit_real
			px,py = horiz2xy(self.azim,self.alti,self.ImageInfo,derotate=True)
			self.preliminary_star = \
				self.skyimage.scatter(px,py,marker='o',c='yellow',alpha=0.5)
			self.skyfigure.canvas.draw()
			self.skyfigure.canvas.flush_events()
			print('Name: %s, Az: %s, Alt: %s' %(self.name,self.azim,self.alti))

		if event.key=='q':
			print('End')
			self.skyfigure.canvas.mpl_disconnect(self.cid_mouse)
			self.skyfigure.canvas.mpl_disconnect(self.cid_keyboard)
			print(self.identified_stars)
			plt.close()
		
		self.astrometry_optimizer(full=(self.star_index>3))

		return(None)
	
	def astrometry_optimizer(self,full=True):
		from scipy.optimize import minimize
		from astrometry import horiz2xy
		
		def horiz2xy_chi2(sol,az,alt,x,y):
			self.ImageInfo.radial_factor     = sol[0]
			self.ImageInfo.azimuth_zeropoint = sol[1]
			if (full==True):
				self.ImageInfo.delta_x           = sol[2]
				self.ImageInfo.delta_y           = sol[3]
				self.ImageInfo.latitude_offset   = sol[4]
				self.ImageInfo.longitude_offset  = sol[5]
			else:
				self.ImageInfo.delta_x           = 0
				self.ImageInfo.delta_y           = 0
				self.ImageInfo.latitude_offset   = 0
				self.ImageInfo.longitude_offset  = 0
			
			xf,yf = horiz2xy(az,alt,self.ImageInfo,derotate=True)
			return(np.sum((xf-x)**2 + (yf-y)**2))
		
		coords = np.array(self.identified_stars)[:,1:] # Remove star name
		coords = np.array(coords,dtype=float)          # Convert to float
		[_az,_alt,_x,_y] = np.transpose(coords)        # Transpose and split
		print('Solving equation system')
		
		if (full==True):
			initial=[10,0,0,0,0,0]
		else:
			initial=[0,0]
		
		res = minimize(horiz2xy_chi2,initial,args = (_az,_alt,_x,_y),tol=1e-3)
		print("Parameters (radial_factor, azimuth_zeropoint, delta_x, delta_y, lat_offset, lon_offset): ")
		print(res.x)
		print("Score [sum(dev^2)] = %.3f" %horiz2xy_chi2(res.x,_az,_alt,_x,_y))
		print("Success: %s" %res.success)

	def astrometry_solver(self):
		print(\
			'*** Star select tool. Press right-click to begin. *** \n'+\
			'Right-click: assign star coords. \n'+\
			'n:           next star (skip current). \n'+\
			'p:           previous star (remove last entry). \n'+\
			'q:           quit star select tool. \n')
		
		self.identified_stars = []
		self.scatter_stars = []
		self.label_stars = []
		self.star_index = 0
		self.completed=0

		# For the northern hemisphere, put Polaris as the first star
		if (self.ImageInfo.latitude>0):
			polaris_index=[Star.HDcode for Star in self.StarCatalog.StarList_Tot].index("HD8890")
			AuxStar = self.StarCatalog.StarList_Tot[polaris_index]
			self.StarCatalog.StarList_Tot[polaris_index] = \
				self.StarCatalog.StarList_Tot[0]
			self.StarCatalog.StarList_Tot[0] = AuxStar
		
		self.name=self.StarCatalog.StarList_Tot[0].name
		self.azim=self.StarCatalog.StarList_Tot[0].azimuth
		self.alti=self.StarCatalog.StarList_Tot[0].altit_real
		print('Name: %s, Az: %s, Alt: %s' %(self.name,self.azim,self.alti))
		
		self.cid_mouse = self.skyfigure.canvas.mpl_connect('button_press_event', self.mouse_press_callback)
		self.cid_keyboard = self.skyfigure.canvas.mpl_connect('key_press_event', self.key_press_callback)
		plt.show(block=True)
	
	def draw_skymap_data(self):
		''' Draw image '''
		self.skyimage.imshow(self.stretched_fits_data,cmap=mpl.cm.gray)

		self.skyimage.axis([0,self.ImageInfo.resolution[0],0,self.ImageInfo.resolution[1]])
		information=str(self.ImageInfo.date_string)+" UTC\n"+str(self.ImageInfo.latitude)+5*" "+\
			str(self.ImageInfo.longitude)+"\n"+self.ImageInfo.used_filter
		
		self.skyimage.text(0.005,0.005,information,fontsize='small',color='white',\
			transform = self.skyimage.transAxes,backgroundcolor=(0,0,0,0.75))
		
		plt.draw()
		
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
		
	def draw_annotate_star(self,Star,full=False):
		# Draw identified stars and measuring circles.
		# Annotate HD catalog code and Magnitude for each star.
		
		if (full==False):
			self.skyimage.scatter(Star.Xcoord,Star.Ycoord,marker='x',c='yellow',alpha=0.2,label='Identified stars')
			self.skyimage.annotate(\
				Star.name,xy=(Star.Xcoord,Star.Ycoord), \
				xycoords='data',xytext=(0, 3),\
				textcoords='offset points',fontsize=8,alpha=0.8)
		else:
			
			self.skyimage.scatter(Star.Xcoord,Star.Ycoord,marker='+',c='red',alpha=0.2,label='Identified stars')
			self.skyimage.add_patch(mpp.Circle(\
				(Star.Xcoord,Star.Ycoord),Star.R1,facecolor='none',edgecolor=(0,0,0.8),\
				linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
			self.skyimage.add_patch(mpp.Circle(\
				(Star.Xcoord,Star.Ycoord),Star.R2,facecolor='none',edgecolor=(0,0.8,0),\
				linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
			self.skyimage.add_patch(mpp.Circle(\
				(Star.Xcoord,Star.Ycoord),Star.R3,facecolor='none',edgecolor=(0.8,0,0),\
				linewidth=1, fill=False, alpha=0.5,label='_nolegend_'))
			self.skyimage.annotate(Star.FilterMag,xy=(Star.Xcoord,Star.Ycoord), xycoords='data',xytext=(0,-10),\
				textcoords='offset points',fontsize=8)
			
	def show_figure(self):
		self.skyimage.legend(('In catalog','Detected'),loc='upper right')
		
		def skymap_filename():
			filename = self.ImageInfo.skymap_path +\
				"/Skymap_"+self.ImageInfo.obs_name+"_"+self.ImageInfo.fits_date+"_"+\
				self.ImageInfo.used_filter+".png"
			return(filename)
		
		if self.ImageInfo.skymap_path=="screen":
			plt.show()
			#self.skyfigure.canvas.draw()
			#self.skyfigure.canvas.flush_events()
		else:
			plt.savefig(skymap_filename(),bbox_inches='tight')
		
		#plt.clf()
		#plt.close('all')
