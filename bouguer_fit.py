#!/usr/bin/env python

'''
Star detection module

Fit fluxes and star data to an extinction law to obtain 
extinction and instrument zeropoint.
____________________________

This module is part of the PyASB project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

__author__ = "Miguel Nievas"
__copyright__ = "Copyright 2012, PyASB project"
__credits__ = ["Miguel Nievas"]
__license__ = "GNU GPL v3"
__version__ = "1.99.0"
__maintainer__ = "Miguel Nievas"
__email__ = "miguelnr89[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


try:
	import matplotlib.pyplot as mpl
	import matplotlib.colors as mpc
	import matplotlib.patches as mpp
	import scipy.stats as stats
	import math
	import numpy as no
except:
	print 'One or more modules missing: numpy,matplotlib'
	raise SystemExit


def bouguer_fit(StarMeasured, ImageInfo, ObsPyephem):
	''' Fit measured fluxes to an extinction model
		Return regression parameters (ZeroPoint, Extinction) '''
	
	class Regression:
		x    = [Star.airmass for Star in StarMeasured]
		y    = [Star.m25logF for Star in StarMeasured]
		yerr = [Star.m25logF_unc for Star in StarMeasured]
	
	try:
		class fixed_zp:
			fixed_y     = ImageInfo.zeropoint
			fixed_y_unc = ImageInfo.zeropoint_unc
		Regression = theil_sen(Regression,fixed_zp)
		return Regression
	except:
		try:
			Regression = theil_sen(Regression)
			return Regression
		except:
			raise

def bouguer_plot(Regression,ImageInfo,ObsPyephem):
	''' Plot photometric data from the bouguer fit '''

	xfit = linspace(1,calculate_airmass(ImageInfo.min_altitude),10)
	yfit = polyval([Regression.slope,Regression.zp],xfit)

	xdata    = Regression.xdata
	ydata    = Regression.ydata
	yerrdata = Regression.yerrdata

	bouguerfigure = figure(figsize=(8,6),dpi=100)
	bouguerplot = figurabouguer.add_subplot(111)
	bouguerplot.set_title('Bouguer extinction law fit',size="xx-large")
	bouguerplot.errorbar(xdata, ydata, yerr=yerrdata, fmt='*', ecolor='g')
	bouguerplot.plot(x_fit,y_fit,'r-')
	
	try:
		plot_infotext = \
			ImageInfo.imagesdate+str(ObsPyephem.lat)+5*" "+str(ObsPyephem.lon)+"\n"+\
			ImageInfo.used_filter+4*" "+"Rcorr="+str("%.3f"%float(Regression.Kendall_tau))+"\n"+\
			"C="+str("%.3f"%float(Regression.zp))+"+/-"+str("%.3f"%float(Regression.zp_err))+"\n"+\
			"K="+str("%.3f"%float(Regression.slope))+"+/-"+str("%.3f"%float(Regression.slope_err))+"\n"+\
			str("%.0f"%(100.*Regression.Nrel))+"% of "+str(Norig)+" photometric measures shown"
		bouguerplot.text(0.1,0.1,plot_infotext,fontsize='x-small',transform = plot_infotext.transAxes)
	except:
		raise
	
	if ImageInfo.bouguerplot_file!=False:
		# Show or save the bouguer plot
		show_or_save_bouguerplot(bouguerfigure,ImageInfo,ObsPyephem)

class TheilSenRegression
	def __init__(self,Xpoints,Ypoints,y0=None, y0err=None, x0=None, x0err=None):
		assert len(Xpoints)==len(Ypoints)
		self.Xpoints = Xpoints
		self.Ypoints = Ypoints
		if y0!=None:
			self.fixed_zp = True
			if y0err!=None:
				self.y0err=y0err
			else: 
				self.y0err=0.0
			
			if x0!=None:
				self.x0 = x0
				if x0err!=None:
					self.x0err = 0.0
			else:
				self.x0 = 0.0
				self.x0err = 0.0
		else: self.fixed_zp = False
		self.Nstars_initial = len(Ypoints)
		self.Nstars_final = self.Nstars_initial
		self.pair_blacklist = []
		# Perform the regression
		self.perform_regression()
		
	def perform_regression(self):
		# Prepare data for regression
		self.build_matrix_values()
		self.build_complementary_matrix()
		self.build_slopes_matrix()
		self.upper_diagonal_slope_matrix_values()
		# Slope
		self.calculate_mean_slope()
		# Zero point
		self.build_zeropoint_array()
		self.calculate_mean_zeropoint()
		# Errors and fit quality
		self.calculate_residuals()
		self.calculate_kendall_tau()
		self.calculate_errors()
		if self.fixed_zp == True:
			self.mean_zeropoint = self.y0
			self.error_zeropoint = self.y0err
	
	def build_matrix_values(self):
		self.X_matrix_values = np.array([self.Xpoints for line in Xpoints])
		self.Y_matrix_values = np.array([self.Ypoints for line in Ypoints])
	
	def build_complementary_matrix(self):
		if self.fixed_zp == False:
			self.X_complementary_values = self.X_matrix_values.transpose()
			self.Y_complementary_values = self.Y_matrix_values.transpose()
		if self.fixed_zp == True:
			self.X_complementary_values = np.array([[self.x0\
				for column in self.Xpoints] for line in self.Xpoints])
			self.Y_complementary_values = np.array([[self.y0\
				for column in self.Ypoints] for line in self.Ypoints])
	
	def build_slopes_matrix(self):
		self.slopes_matrix = \
			((self.Y_matrix_values-self.Y_complementary_values)/ \
			(self.X_matrix_values-self.X_complementary_values))
	
	def upper_diagonal_slope_matrix_values(self):
		self.upper_diag_slopes = \
			np.array([self.slopes_matrix[l][c] \
			for l in xrange(len(self.slopes_matrix)) \
			for c in xrange(len(self.slopes_matrix)) if c>l])
	
	def calculate_mean_slope(self):
		self.mean_slope  = np.median(self.upper_diag_slopes)
		
	def build_zeropoint_array(self):
		self.zeropoint_array = self.Ypoints - self.Xpoints*self.median_slope
	
	def calculate_mean_zeropoint(self):
		self.mean_zeropoint  = np.median(self.zeropoint_array)
		
	def calculate_residuals(self):
		self.residuals = self.zeropoint_array-self.mean_zeropoint
		
	def calculate_errors(self):
		xmedcuad = np.median(self.Xpoints)**2
		xcuaddif = self.Xpoints**2 - xmedcuad
		xdensity = np.sum(xcuad)
		sigma2_res = (1./(self.Nstars_final-2))*self.residuals
		sigma2_slope = sigma2_res/abs(xdensity)
		sigma2_int = sigma2_res*(1./self.Nstars_final + xmedcuad/abs(xdensity))
		
		self.error_slope = stats.t.ppf(0.975,self.Nstars_final-2) * math.sqrt(sigma2_slope)
		self.error_zeropoint = stats.t.ppf(0.975,self.Nstars_final-2) * math.sqrt(sigma2_int)
	
	def calculate_kendall_tau(self):
		self.kendall_tau = \
			(len(self.upper_diag_slopes[self.upper_diag_slopes>0]) - \
			len(self.upper_diag_slopes[self.upper_diag_slopes<0])) / \
			len(self.upper_diag_slopes)
			

