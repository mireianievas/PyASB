#!/usr/bin/env python

'''
Star detection module

Fit fluxes and star data to an extinction law to obtain 
extinction and instrument zeropoint.
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
except:
	print 'One or more modules missing: pyfits,HeaderTest'
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
			fixed_y     = ImageInfo.Config.y
			fixed_y_unc = ImageInfo.Config.y_unc
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

	xfit = linspace(1,calculate_airmass(ImageInfo.Config.min_altitude),10)
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
	
	if ImageInfo.Config.bouguerplot_file!=False:
		# Show or save the bouguer plot
		show_or_save_bouguerplot(bouguerfigure,ImageInfo,ObsPyephem)



	


