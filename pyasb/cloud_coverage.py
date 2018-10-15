

'''
Cloud covering module

Determine if clouds are present in the image.
____________________________

This module is part of the PyASB project, 
created and maintained by Mireia Nievas [UCM].
____________________________
'''

DEBUG = False

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
    import copy
    import numpy as np
    import warnings
    import scipy.interpolate as sint
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.colors as mpc
    import matplotlib.cm as mpcm
    import matplotlib.patches as mpp
    from star_calibration import StarCatalog
except:
    print(str(inspect.stack()[0][2:4][::-1])+': One or more modules missing')
    raise SystemExit

warnings.simplefilter("ignore", category=RuntimeWarning)

class CloudCoverage():
    ''' This is a completely optional module.
        If user doesn't want neither the map nor the table, 
        just don't do anything'''
    def __init__(self,Image,ImageAnalysis,BouguerFit):
        ''' Calculate the mean cloud covering (at the whole image) '''
        #Old method, directly from the % stars. Gives higher values
        #TotalStars = len(ImageAnalysis.StarCatalog.StarList_woPhot)
        #TotalStarsWithPhot = len(ImageAnalysis.StarCatalog.StarList)

        #TotalStars = 0;
        #TotalStarsWithPhot = 0;
        ImageAnalysis.StarCatalog.look_for_nearby_stars(Image.FitsImage,Image.ImageInfo)
        
        '''        
        for Star in ImageAnalysis.StarCatalog.StarList_TotVisible:
            #if Star.masked==True: continue
            #elif Star.saturated == True: continue
            #elif Star.cold_pixels==True: continue
            #if Star.altit_real<30: continue
            # Star is not saturated and present in the sky, add it [weighted]
            TotalStars += 1
            if Star not in ImageAnalysis.StarCatalog.StarList_WithNearbyStar: continue
            # Star has photometric measures, add it [weighted]
            TotalStarsWithPhot += 1
        '''
        
        #TotalPercStars = (TotalStarsWithPhot+1e-6)*1./(TotalStars+1e-6)
        TotalStars     = len(ImageAnalysis.StarCatalog.StarList_TotVisible)
        TotalDetected  = len(ImageAnalysis.StarCatalog.StarList_WithNearbyStar) 
        TotalPercStars = TotalDetected/TotalStars
        
        # Calculate the mean Cloud Coverage value (AllSky)
        self.mean_cloudcover = self.cloud_coverage_value(\
          np.array([TotalPercStars]))[0]
        
        # Binomial error for the All Sky Cloud Coverage value.
        self.error_cloudcover = self.cloud_coverage_error(\
          self.mean_cloudcover,
          np.array([TotalPercStars]),
          np.array([TotalStars]))[0]
        
        print("Mean cloud coverage: %.3f+/-%.3f" \
         %(self.mean_cloudcover,self.error_cloudcover))
        
        ''' If asked, calculate the cloud covering by sectors in the image '''
        try:
            assert(\
                Image.ImageInfo.clouddata_path!=False or\
                Image.ImageInfo.cloudmap_path!=False)
        except Exception as e:
            #print(inspect.stack()[0][2:4][::-1])
            print('Skipping cloud coverage detection')
            #print(type(e))
            #print(e)
            return(None)
        else:
            print('Measuring Cloud Covering ...')
            self.create_bins()
            #self.star_detection(Image)
            self.StarCatalog = ImageAnalysis.StarCatalog
            self.cloud_coverage(
                self.StarCatalog.StarList_WithNearbyStar,
                self.StarCatalog.StarList_TotVisible,BouguerFit)
            self.cloud_map(BouguerFit,ImageInfo=Image.ImageInfo)
            self.clouddata_table(ImageInfo=Image.ImageInfo)
    
    def star_detection(self,Image):
        # relax requisites to get more stars
        #normal_stdout = sys.stdout 
        #sys.stdout = open(os.devnull,'w')
        II = copy.deepcopy(Image.ImageInfo)
        II.baseflux_detectable = II.baseflux_detectable/5.
        II.base_radius   *= 5
        #II.max_magnitude += 0.5
        II.min_altitude   = 0
        II.skymap_path    = False
        
        self.StarCatalog = StarCatalog(ImageInfo=II)
        self.StarCatalog.process_catalog_specific(Image.FitsImage,II)
        #sys.stdout = normal_stdout
        
    def create_bins(self):
        ''' Returns binned sky '''
        self.azseparation  = 45
        self.zdseparation = 20
        
        self.AZdirs  = np.arange(0,360+1,self.azseparation)
        self.ZDdirs = np.arange(0,90+1,self.zdseparation)
        self.ALTdirs = 90-self.ZDdirs
        
        self.AZgrid,self.ZDgrid = np.meshgrid(self.AZdirs,self.ZDdirs)
        self.ALTgrid = 90.-self.ZDgrid
        self.AZgrid = self.AZgrid*np.pi/180.0
    
    @staticmethod
    def cloud_coverage_value(PercentageStars):
        # Stretch the data and estimate the coverage
        TheCloudCoverage = 1.0-PercentageStars
        #self.CloudCoverage = self.CloudCoverage/0.2
        TheCloudCoverage[TheCloudCoverage<0.0]=0.0 # Cloudless
        TheCloudCoverage[TheCloudCoverage>1.0]=1.0 # Overcast
        return(TheCloudCoverage)
    
    @staticmethod
    def cloud_coverage_error(TheCloudCoverage,PercentageStars,TotalStars):
        '''
        Estimate the cloud coverage uncertainty as binomial distribution
        '''
        TheCloudCoverageErr = \
          TheCloudCoverage*(1-TheCloudCoverage)/np.sqrt(TotalStars)
        return(TheCloudCoverageErr)
    
    def cloud_coverage(self,StarList_Photom,StarList_woPhot,BouguerFit):
        Stars_in_field = [Star for Star in StarList_woPhot]
        
        # The number of predicted/observable stars on a region
        PredictedStars = np.zeros((len(self.ALTdirs),len(self.AZdirs)))
        # The number of detected stars on that region 
        ObservedStars  = np.zeros((len(self.ALTdirs),len(self.AZdirs)))
        
        # Sum of flux percentage (~ mean extinction in absolute units)
        PercentageFlux    = np.zeros((len(self.ALTdirs),len(self.AZdirs)))
        PercentageFluxErr = np.zeros((len(self.ALTdirs),len(self.AZdirs)))
        
        minimum_stars = 3;
        
        for column in xrange(len(self.AZdirs)):
            for line in xrange(len(self.ALTdirs)):
                for Star in Stars_in_field: # Here we should divide /2, /1 will smooth data
                    # See if the Star is in the given azimuth and altitude
                    if (\
                      (abs(Star.altit_real-self.ALTdirs[line])>self.zdseparation) or\
                      (abs(Star.azimuth-self.AZdirs[column])>self.azseparation and\
                       abs(360-abs(Star.azimuth-self.AZdirs[column]))>self.azseparation and\
                       90-Star.altit_real>self.zdseparation)
                       ):
                        continue
                    
                    if Star.saturated == True: continue
                    elif Star.cold_pixels == True: continue
                    elif Star.masked == True: continue
                    
                    # Star is not saturated and present in the sky, add it [weighted]
                    PredictedStars[line][column] += 1
                    
                    if Star not in StarList_Photom:
                        continue
                    
                    # Star has photometric measures, add it [weighted]
                    ObservedStars[line][column] += 1
                    
                    # Alternative estimation based on flux extinction. To be completed.
                    Predicted_Flux      = \
                     10**(0.4*(BouguerFit.Regression.mean_zeropoint-Star.FilterMag))
                    Predicted_FluxError =\
                     Predicted_Flux*np.log(10)*0.4*BouguerFit.Regression.error_zeropoint
                    Measured_Flux       = Star.starflux
                    Measured_FluxError  = Star.starflux_err
                    
                    PercentageFlux[line][column] += np.clip(Measured_Flux/Predicted_Flux,0,1)
                
                # If there are not enough stars in the field, truncate the measure
                if PredictedStars[line][column] < minimum_stars:
                    ObservedStars[line][column] = 0;
                    PredictedStars[line][column] = 0+1e-6;
                    
        
        # Normalization of flux percentages
        PercentageFlux = PercentageFlux*1./(1e-5 + ObservedStars)
        
        PercentageStars = \
            np.array((0+ObservedStars*1.0)/(1e-5+PredictedStars*1.0))
        
        self.CloudCoverage = self.cloud_coverage_value(PercentageStars)
        self.CloudCoverageErr = self.cloud_coverage_error(self.CloudCoverage,PercentageStars,PredictedStars)
        self.CloudCoverage[PredictedStars<2] = None # not enough stars
    
    def clouddata_table(self,ImageInfo):
        try:
            assert(ImageInfo.clouddata_path!=False)
        except:
            print(inspect.stack()[0][2:4][::-1])
            print('Skipping write clouddata table to file')
        else:
            print('Write clouddata table to file')
            header = '#Altitude\Azimuth'
            for az_ in  self.AZdirs:
                header += ', '+str(az_)
            header+='\n'
            
            content = [header]
            
            for k,alt_ in enumerate(self.ALTdirs):
                line = str(alt_)
                for j in xrange(len(self.AZdirs)):
                    line += ', '+str("%.3f +/- %.3f" %\
                     (float(self.CloudCoverage[k][j]), \
                     float(self.CloudCoverageErr[k][j])))
                line+='\n'
                content.append(line)
            
            if ImageInfo.clouddata_path == "screen":
                print(content)
            else:
                cloudtable_filename = str("%s/CloudTable_%s_%s_%s.txt" %(\
                    ImageInfo.clouddata_path, ImageInfo.obs_name,\
                    ImageInfo.fits_date,ImageInfo.used_filter))
                cloudfile = open(cloudtable_filename,'w+')
                cloudfile.writelines(content)
                cloudfile.close()
    
    
    def cloud_map(self,BouguerFit,ImageInfo):
        try:
            assert(ImageInfo.cloudmap_path!=False)
        except:
            print(inspect.stack()[0][2:4][::-1])
            print('Skipping write cloudmap to file')
            return(None)
        else:
            print('Output cloudmap')
        
        ''' Create the cloud map '''
        self.Cloudfigure = plt.figure(figsize=(8,7.5))
        self.Cloudgraph  = self.Cloudfigure.add_subplot(111,projection='polar')
        
        
        # Grid and interpolate data
        self.AZgridi,self.ZDgridi = np.mgrid[0:2*np.pi:1000j, 0:90:1000j]
        self.ALTgridi = 90. - self.ZDgridi
        coord_reshape = np.array([[self.AZgrid[j][k],self.ZDgrid[j][k]] \
            for k in xrange(len(self.AZgrid[0])) for j in xrange(len(self.AZgrid))])
        data_reshape = np.array([self.CloudCoverage[j][k] \
            for k in xrange(len(self.AZgrid[0])) for j in xrange(len(self.AZgrid))])
        self.CloudCoveragei = sint.griddata(coord_reshape,data_reshape, \
            (self.AZgridi,self.ZDgridi), method='nearest')
        
        # Colormap
        #cloud_cmap = mpcm.get_cmap('gray', 5)
        
        cdict = {
         'red': [(0.0,0.2,0.2),
                 (1.0,0.8,0.8)],
         'green': [(0.0,0.2,0.2),
                 (1.0,0.8,0.8)],
         'blue': [(0.0,0.2,0.2),
                 (1.0,0.8,0.8)]}
        
        cloud_cmap = mpc.LinearSegmentedColormap('gray_colormap',cdict,N=5)
        
        # Create the graph
        self.ColorMesh = self.Cloudgraph.pcolormesh(\
            self.AZgridi, 
            self.ZDgridi, 
            self.CloudCoveragei,
            vmin=0,
            vmax=1,
            cmap=cloud_cmap)
        
        # Ticks
        def ticks_and_locators():
            ''' Add ticks to the graph '''
            radial_locator = np.arange(10,90+1,10)
            radial_label = ["$80$","$70$","$60$","$50$","$40$","$30$","$20$","$10$","$0$"]
            theta_locator = np.arange(0,360,45)
            theta_label = ["$N$","$NE$","$E$","$SE$","$S$","$SW$","$W$","$NW$"]
            
            self.Cloudgraph.set_rgrids(radial_locator,radial_label,\
             size="large",color='k',alpha=0.75)
            self.Cloudgraph.set_thetagrids(theta_locator,theta_label,size="large")
            # rotate the graph (North up)
            self.Cloudgraph.set_theta_direction(-1)
            self.Cloudgraph.set_theta_offset(np.pi/2)
        
        self.Cloudgraph.grid(True)
        
        # Colorbar
        def color_bar():
            ''' Add the colorbar '''
            # Separation between colour bar and graph
            self.Cloudfigure.subplots_adjust(right=1)
            # Color bar 
            self.Cloudcolorbar = plt.colorbar(\
             self.ColorMesh,orientation='vertical',pad=0.07,shrink=0.75)
            self.Cloudfigure.subplots_adjust(right=0.80) # Restore separation
            #self.ColorMesh.set_clim(0.0,1.0)
            self.Cloudcolorbar.set_ticks(np.arange(0,1+1e-6,0.1))
            self.Cloudcolorbar.set_label("Cloud Coverage",rotation="vertical",size="large")
        
        
        # Ticks and colorbar
        ticks_and_locators()
        color_bar()
        
        # Information text on image
        self.Cloudgraph.text(0,np.max(self.ZDgridi)+15, unicode(ImageInfo.cloudmap_title, 'utf-8'),\
            horizontalalignment='center',size='xx-large')
        # Image information
        image_information = str(ImageInfo.date_string)+" UTC\n"+str(ImageInfo.latitude)+5*" "+\
            str(ImageInfo.longitude)+"\n"+ImageInfo.used_filter+4*" "+\
            "K="+str("%.3f" % float(BouguerFit.Regression.extinction))+"+-"+\
            str("%.3f" % float(BouguerFit.Regression.error_extinction))+"\n"+\
            str("Cloud coverage: %.3f +/- %.3f" %\
             (float(self.mean_cloudcover),float(self.error_cloudcover)))
        
        self.Cloudgraph.text(5*np.pi/4,145,unicode(image_information,'utf-8'),fontsize='x-small')
        
        if ImageInfo.cloudmap_path=="screen":
            plt.show()
        else:
            cloudmap_filename = str("%s/CloudMap_%s_%s_%s.png" %(\
                ImageInfo.cloudmap_path, ImageInfo.obs_name,\
                ImageInfo.fits_date, ImageInfo.used_filter))
            plt.tight_layout(pad=-1.5,rect=[0.1,0.05,1.,0.95])
            plt.savefig(cloudmap_filename)
        
        #plt.clf()
        #plt.close('all')
        
        
        
        
    
