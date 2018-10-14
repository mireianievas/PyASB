#!/usr/bin/env python

'''
Bouguer fitting module

Fit fluxes and star data to an extinction law to obtain
extinction and instrument zeropoint.
____________________________

This module is part of the PyASB project,
created and maintained by Mireia Nievas [UCM].
____________________________
'''

DEBUG=False

__author__ = "Mireia Nievas"
__copyright__ = "Copyright 2012, PyASB project"
__credits__ = ["Mireia Nievas"]
__license__ = "GNU GPL v3"
__shortname__ = "PyASB"
__longname__ = "Python All-Sky Brightness pipeline"
__version__ = "1.99.0"
__maintainer__ = "Mireia Nievas"
__email__ = "mirph4k[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


try:
    import sys,os,inspect
    import matplotlib.pyplot as plt
    import matplotlib.colors as mpc
    import matplotlib.patches as mpp
    import scipy.stats as stats
    import math
    import numpy as np
    import astrometry
except:
    print(str(inspect.stack()[0][2:4][::-1])+': One or more modules missing')
    raise SystemExit

class BouguerFit():
    def __init__(self,ImageInfo,PhotometricCatalog):
        print('Calculating Instrument zeropoint and extinction ...')
        self.can_continue = False
        # Get Zero Point from ImageInfo and Stars from the Catalog
        self.bouguer_fixedy(ImageInfo)
        self.bouguer_data(PhotometricCatalog)
        # Set the default values
        self.bouguer_setdefaults(ImageInfo)
        # Try to fit the data
        try:
            self.bouguer_fit(ImageInfo)
            self.bouguer_plot(ImageInfo)
        except:
            print(str(inspect.stack()[0][2:4][::-1])+\
             ' cannot fit the data. Npoints='+str(len(self.ydata)))
            assert(self.can_continue == True)
            raise
        else:
            self.can_continue = True
        
        assert(self.can_continue == True)
        print("Bouguer extinction fit results: \n"+\
         " -> C=%.3f+/-%.3f, K=%.3f+/-%.3f, r=%.3f" \
         %(self.Regression.mean_zeropoint,self.Regression.error_zeropoint,\
           self.Regression.extinction,self.Regression.error_extinction,\
           self.Regression.kendall_tau))

    def bouguer_data(self,StarCatalog):
        ''' Get Star data from the catalog '''
        self.xdata = np.array([Star.airmass \
         for Star in StarCatalog.StarList_Phot])
        self.ydata = np.array([Star.m25logF \
         for Star in StarCatalog.StarList_Phot])
        self.yerr  = np.array([Star.m25logF_unc \
         for Star in StarCatalog.StarList_Phot])

    def bouguer_fixedy(self,ImageInfo):
        ''' Try to get the fixed Y (zero point)'''
        try:
            self.fixed_y     = ImageInfo.used_zero_point[0]
            self.fixed_y_unc = ImageInfo.used_zero_point[1]
            self.yisfixed=True
        except:
            if DEBUG==True: print(str(inspect.stack()[0][2:4][::-1])+' dont fix the Zero Point')
            self.yisfixed=False

    def bouguer_setdefaults(self,ImageInfo):
        ''' Set default values (just in case that the bouguer fit fails '''
        if self.yisfixed == True:
            class Regression:
                mean_zeropoint  = self.fixed_y
                error_zeropoint = self.fixed_y_unc
                mean_slope = 10.0
                error_slope = 10.0
                extinction  = 10.0
                error_extinction = 10.0
                kendall_tau = 0.0
                Nstars_initial = 0
                Nstars_final = 0
                Nstars_rel = 0
            self.Regression = Regression()
            self.can_continue = True

    def bouguer_fit(self,ImageInfo):
        '''
        Fit measured fluxes to an extinction model
        Return regression parameters (ZeroPoint, Extinction)
        '''

        if self.yisfixed:
            self.Regression = TheilSenRegression(\
             Xpoints = self.xdata,\
             Ypoints = self.ydata,\
             ImageInfo = ImageInfo,\
             y0 = self.fixed_y,\
             y0err = self.fixed_y_unc)
        else:
            try:
                self.Regression = TheilSenRegression(\
                 Xpoints = self.xdata,\
                 Ypoints = self.ydata,\
                 ImageInfo = ImageInfo)
            except:
                print(inspect.stack()[0][2:4][::-1]);
                raise

        # Apply bad point filter to data
        self.xdata = self.xdata[self.Regression.badfilter]
        self.ydata = self.ydata[self.Regression.badfilter]
        self.yerr = self.yerr[self.Regression.badfilter]

    def bouguer_plot(self,ImageInfo):
        if ImageInfo.bouguerfit_path==False:
            # Don't draw anything
            print('Skipping BouguerFit Graph')
            return(None)

        ''' Plot photometric data from the bouguer fit '''

        xfit = np.linspace(1,astrometry.calculate_airmass(ImageInfo.min_altitude),10)
        yfit = np.polyval([self.Regression.mean_slope,self.Regression.mean_zeropoint],xfit)

        bouguerfigure = plt.figure(figsize=(8,6))
        bouguerplot = bouguerfigure.add_subplot(111)
        bouguerplot.set_title('Bouguer extinction law fit\n',size="xx-large")
        bouguerplot.set_xlabel('Airmass')
        bouguerplot.set_ylabel(r'$m_0+2.5\log_{10}(F)$',size="large")
        bouguerplot.errorbar(self.xdata, self.ydata, yerr=self.yerr, fmt='*', ecolor='g')
        bouguerplot.plot(xfit,yfit,'r-')

        try:
            plot_infotext = \
                ImageInfo.date_string+"\n"+str(ImageInfo.latitude)+5*" "+str(ImageInfo.longitude)+"\n"+\
                ImageInfo.used_filter+4*" "+"Rcorr="+str("%.3f"%float(self.Regression.kendall_tau))+"\n"+\
                "C="+str("%.3f"%float(self.Regression.mean_zeropoint))+\
                "+/-"+str("%.3f"%float(self.Regression.error_zeropoint))+"\n"+\
                "K="+str("%.3f"%float(self.Regression.extinction))+"+/-"\
                +str("%.3f"%float(self.Regression.error_slope))+"\n"+\
                str("%.0f"%(self.Regression.Nstars_rel))+"% of "+\
                str(self.Regression.Nstars_initial)+" photometric measures shown"
            bouguerplot.text(0.05,0.05,plot_infotext,fontsize='x-small',transform = bouguerplot.transAxes)
        except:
            print(inspect.stack()[0][2:4][::-1])
            raise

        # Show or save the bouguer plot
        if ImageInfo.bouguerfit_path=="screen":
            plt.show()
        else:
            bouguer_filename = str("%s/BouguerFit_%s_%s_%s.png" %(\
                ImageInfo.bouguerfit_path, ImageInfo.obs_name,\
                ImageInfo.fits_date, ImageInfo.used_filter))
            plt.tight_layout(pad=0)
            plt.savefig(bouguer_filename,bbox_inches='tight')

        #plt.clf()
        #plt.close('all')

class TheilSenRegression():
    # Robust Theil Sen estimator, instead of the classic least-squares.
    def __init__(self,Xpoints,Ypoints,ImageInfo,y0=None,y0err=None,x0=None,x0err=None):
        assert(len(Xpoints)==len(Ypoints) and len(Ypoints)>2)
        self.Xpoints = np.array(Xpoints)
        self.Ypoints = np.array(Ypoints)
        if y0!=None:
            self.fixed_zp = True
            self.y0 = y0
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
        self.Nstars_initial = len(self.Ypoints)
        self.Nstars_final = self.Nstars_initial
        self.pair_blacklist = []
        # Perform the regression
        self.perform_regression()
        # Delete bad points
        self.delete_bad_points(ImageInfo)
        # Try to improve the regression with filtered data
        self.perform_regression()

        self.Nstars_final = sum(self.badfilter)
        self.Nstars_rel = 100.*self.Nstars_final/self.Nstars_initial

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

        self.Nstars_final = len(self.Ypoints)

    def build_matrix_values(self):
        self.X_matrix_values = \
            np.array([[column for column in self.Xpoints] for line in self.Xpoints])
        self.Y_matrix_values = \
            np.array([[line for line in self.Ypoints] for line in self.Ypoints])

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
            ((self.Y_matrix_values-self.Y_complementary_values +1e-20)/ \
            (self.X_matrix_values-self.X_complementary_values +1e-20))
        # +1e-20 lets us hide Numpy warning with 0/0

    def upper_diagonal_slope_matrix_values(self):
        self.upper_diag_slopes = \
            np.array([self.slopes_matrix[l][c] \
            for l in xrange(len(self.slopes_matrix)) \
            for c in xrange(len(self.slopes_matrix[0])) if c>l])

    def calculate_mean_slope(self):
        self.mean_slope = np.median(self.upper_diag_slopes)
        self.extinction = -self.mean_slope

    def build_zeropoint_array(self):
        self.zeropoint_array = self.Ypoints - self.Xpoints*self.mean_slope

    def calculate_mean_zeropoint(self):
        self.mean_zeropoint  = np.median(self.zeropoint_array)

    def calculate_residuals(self):
        self.residuals = self.zeropoint_array-self.mean_zeropoint

    def delete_bad_points(self,ImageInfo):
        # 3*std_residuals threshold
        std_residual = np.std(self.residuals)
        self.badfilter = np.abs(self.residuals)<np.abs(ImageInfo.lim_Kendall_tau*std_residual)
        self.Xpoints = self.Xpoints[self.badfilter]
        self.Ypoints = self.Ypoints[self.badfilter]

    def calculate_errors(self):
        xmedcuad = np.median(self.Xpoints)**2
        xcuaddif = self.Xpoints**2 - xmedcuad
        xdensity = np.sum(xcuaddif)
        sigma2_res = (1./(self.Nstars_final-2))*abs(sum(self.residuals))
        sigma2_slope = sigma2_res/abs(xdensity)
        sigma2_int = sigma2_res*(1./self.Nstars_final + 1.*xmedcuad/abs(xdensity))

        self.error_slope = stats.t.ppf(0.975,self.Nstars_final-2) * math.sqrt(sigma2_slope)
        self.error_zeropoint = stats.t.ppf(0.975,self.Nstars_final-2) * math.sqrt(sigma2_int)
        self.error_extinction = self.error_slope

    def calculate_kendall_tau(self):
        self.kendall_tau = \
            (1.*np.sum(self.upper_diag_slopes>0)-1.*np.sum(self.upper_diag_slopes<0))\
            /(1.*np.size(self.upper_diag_slopes))
