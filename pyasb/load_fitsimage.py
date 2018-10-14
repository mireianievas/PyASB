#!/usr/bin/env python

'''
Load FITS image and header

This module loads the AllSky FITS image and returns both
the Image binary data and the plain-text header.
____________________________

This module is part of the PyASB project,
created and maintained by Mireia Nievas [UCM].
____________________________
'''

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
    import numpy as np
    import astropy.io.fits as pyfits
    from astrometry import ImageCoordinates
    from read_config import *
except:
    print(str(inspect.stack()[0][2:4][::-1])+': One or more modules missing')
    raise SystemExit

class ImageTest():
    '''Perform some test on the image header and extract information'''

    @staticmethod
    def correct_exposure(file_header):
        # Exposure
        try: texp = float(file_header['EXPOSURE'])
        except:
            print(inspect.stack()[0][2:4][::-1])
            raise
        else:
            #assert texp>0., '0s exposure time detected.'
            return texp

    @staticmethod
    def correct_date(file_header):
        # Date and time
        try: date = file_header['DATE']
        except:
            print(inspect.stack()[0][2:4][::-1])
            raise
        else:
            assert len(date)==15 and date[8]=="_", 'Date format not YYYYMMDD_HHMMSS'
            return date

    @staticmethod
    def correct_resolution(file_header):
        # Resolution
        try: resolution = [int(file_header['NAXIS1']),int(file_header['NAXIS2'])]
        except:
            print(inspect.stack()[0][2:4][::-1])
            raise
        else:
            assert resolution[0]>0 and resolution[1]>0, 'Matrix not 2 dimensional'
            return resolution

    @staticmethod
    def correct_filter(file_header):
        # Test if there's a known filter
        try: used_filter = file_header['FILTER']
        except:
            print(inspect.stack()[0][2:4][::-1])
            raise
        else:
            # due to an inconsistent format in AstMon,
            # we found 4 possible formats 'Jonhson_V','JohnsonV','Johnson_V','JonhsonV'
            used_filter = used_filter.replace('_','')
            assert used_filter[0:7] in ['Johnson','Jonhson'], 'Filter type not Johnson'
            assert used_filter[7:] in ['U','B','V','R','I','common'], 'Filter not U,B,V,R or I'
            return 'Johnson_'+used_filter[7:]


class FitsImage(ImageTest):
    def __init__(self,input_file):
        self.load_science(input_file)
        # Backup original data
        print('Backup original (non-calibrated) data')
        self.fits_data_notcalibrated = np.array(self.fits_data)

    def load_science(self,input_file):
        print('Loading ScienceFrame ['+str(input_file)+'] ...'),
        try:
            file_opened = pyfits.open(input_file)
            self.fits_data   = file_opened[0].data
            self.fits_Header = file_opened[0].header
            self.fits_Texp   = float(ImageTest.correct_exposure(self.fits_Header))
        except:
            print(inspect.stack()[0][2:4][::-1])
            raise
        else:
            print('OK')

    def load_dark(self,MasterDark):
        print('Loading MasterDark ...'),
        try:
            MasterDark_HDU    = pyfits.open(MasterDark)
            self.MasterDark_Data   = MasterDark_HDU[0].data
            self.MasterDark_Header = MasterDark_HDU[0].header
            self.MasterDark_Texp   = float(ImageTest.correct_exposure(self.MasterDark_Header))
        except:
            print(inspect.stack()[0][2:4][::-1])
            raise
        else: print('OK')

    def load_flat(self,MasterFlat):
        print('Loading MasterFlat ...'),
        try:
            MasterFlat_HDU    = pyfits.open(MasterFlat)
            self.MasterFlat_Data   = MasterFlat_HDU[0].data
            # Normalize MasterFlat
            self.MasterFlat_Data = self.MasterFlat_Data / np.mean(self.MasterFlat_Data)
            self.MasterFlat_Header = MasterFlat_HDU[0].header
            self.MasterFlat_Texp   = float(ImageTest.correct_exposure(self.MasterFlat_Header))
        except:
            print(inspect.stack()[0][2:4][::-1])
            raise
        else: print('OK')

    def load_bias(self,MasterBias):
        print('Loading MasterBias ...'),
        try:
            MasterBias_HDU    = pyfits.open(MasterBias)
            self.MasterBias_Data   = MasterBias_HDU[0].data
            self.MasterBias_Header = MasterBias_HDU[0].header
            self.MasterBias_Texp   = float(ImageTest.correct_exposure(self.MasterBias_Header))
        except:
            print(inspect.stack()[0][2:4][::-1])
            raise
        else: print('OK')

        

    def reduce_science_frame(self,MasterDark=None,MasterFlat=None,MasterBias=None,ImageInfo=None):
        '''
        Load MasterDark and MasterFlat. MasterBias is neccesary only if working
        with different exposures between Dark and Science frames
        '''

        skip_dark = False
        skip_flat = False


                ### Load FLAT Field
        try:
            self.load_flat(MasterFlat)
        except:
            print(str(inspect.stack()[0][2:4][::-1])+\
            ' WARNING: MasterFlat cannot be loaded, SKIP the flat calibration')
            skip_flat = True
        else:
            skip_flat = False

                ### Load DARK Frame
        try:
            self.load_dark(MasterDark)
        except:
            ''' Try to use MasterDark as a fixed offset value '''
            try:
                self.SyntDark_Data = float(MasterDark)
            except:
                #raise
                print(str(inspect.stack()[0][2:4][::-1])+\
                 ' WARNING: MasterDark cannot be loaded, SKIP the dark calibration')
                skip_dark = True
            else:
                print(str(inspect.stack()[0][2:4][::-1])+\
                 ' WARNING: MasterDark used as a fixed offset value.\n'+\
                 ' Its *STRONGLY* recommended to use a proper MasterDark')
                skip_dark = False
        else:
            if self.MasterDark_Texp == self.fits_Texp:
                self.SyntDark_Data   = self.MasterDark_Data
                self.SyntDark_Texp   = self.MasterDark_Texp
                self.SyntDark_Header = self.MasterDark_Header
            elif self.MasterDark_Texp != self.fits_Texp and MasterBias==None:
                if MasterBias==None:
                    print("WARNING: Science and Dark don't have the same exposure ! ")
                    print('Science_Texp='+str(self.fits_Texp)+'; Dark_Texp='+str(self.MasterDark_Texp))
                self.SyntDark_Data   = self.MasterDark_Data
                self.SyntDark_Texp   = self.MasterDark_Texp
                self.SyntDark_Header = self.MasterDark_Header
            elif self.MasterDark_Texp != self.fits_Texp and MasterBias!=None:
                self.load_bias(MasterBias)
                print('Creating synthetic Dark ...'),
                try:
                    self.SyntDark_Data = (self.MasterDark_Data-self.MasterBias_Data)/ \
                     (self.MasterDark_Texp-self.MasterBias_Data) *\
                     (self.ScienceFrame_Texp-self.MasterBias_Texp)+\
                     self.MasterBias_Data
                    self.SyntDark_Texp   = self.fits_Texp
                    self.SyntDark_Header = self.MasterDark_Header
                    self.SyntDark_Header['EXPOSURE'] = self.SyntDark_Texp
                except:
                    print(inspect.stack()[0][2:4][::-1])
                    raise
                else: print('OK')

            skip_dark = False

        print('Calibrating image with MasterFlat and MasterDark ...'),
    
        # Subtract dark frame
        if skip_dark == False:
            self.fits_data = self.fits_data - self.SyntDark_Data
        
        # Subtract background / bias (measure it in the non-illuminated corners of the image).
        try: assert(self.subtract_corners_background == True and ImageInfo!=None)
        except: pass
        else:
            ImageCoordinates_ = ImageCoordinates(ImageInfo)
            data_corners = self.fits_data[ImageCoordinates_.altitude_map<-20]
            self.bias_image_median = np.median(data_corners)
            self.bias_image_std    = np.std(data_corners)
            self.bias_image_err    = self.bias_image_std/np.sqrt(np.size(data_corners))
            self.fits_data = self.fits_data-self.bias_image_median
            print("Removed: %.2f +/- %.2f counts from measured background" \
             %(self.bias_image_median,self.bias_image_err))
            
            if ImageInfo.summary_path not in [ False, "False", "false", "F", "screen" ]:
                if not os.path.exists(ImageInfo.summary_path):
                    os.makedirs(ImageInfo.summary_path)
                measured_bias_log = open(ImageInfo.summary_path+'/measured_image_bias.txt','a+')
                text_to_log = str(ImageInfo.date_string)+','+str(ImageInfo.used_filter)+','+\
                 str(self.bias_image_median)+','+str(self.bias_image_err)+'\r\n'
                measured_bias_log.write(text_to_log)
                measured_bias_log.close()
            
        
        # Flat field correction
        if skip_flat == False:
            # Skip flat correction for points with <10% illumination?
            #self.MasterFlat_Data[self.MasterFlat_Data<np.mean(self.MasterFlat_Data)/10.]=1.
            self.fits_data = self.fits_data/self.MasterFlat_Data

        print('OK')

    def flip_image_if_needed(self,ImageInfo):
        if (ImageInfo.flip_image==True):
            try:
                self.fits_data_notcalibrated = np.fliplr(self.fits_data_notcalibrated)
            except:
                print('Warning. Cannot flip raw image as requested')

            try:
                self.fits_data = np.fliplr(self.fits_data)
            except:
                print('Warning. Cannot flip calibrated image as requested')



    def __clear__(self):
        backup_attributes = [\
            "fits_data","fits_Header","fits_data_notcalibrated"]

        for atribute in list(self.__dict__):
            #if atribute[0]!="_" and atribute not in backup_attributes:
            if atribute not in backup_attributes:
                del vars(self)[atribute]

