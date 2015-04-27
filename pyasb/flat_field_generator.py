#!/usr/bin/env python

'''
Synthetic FlatField generator

Create AllSky FlatFields from:
 - Low freq flatfields + High freq flatfields
 - Low freq flatfields only
 - Radial profiles + High freq flatfields
 - Radial profiles only
____________________________

This module is part of the PyASB project, 
created and maintained by Miguel Nievas [UCM].
____________________________
'''

DEBUG = False

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
    import signal
    import numpy as np
    import scipy
    import scipy.interpolate
    import scipy.ndimage
    import astropy.io.fits as pyfits
    import re
    # Aux functions from PyASB
    from fits_operator import *
    from read_config import *
    from image_info import *
    from astrometry import *
except:
    print(str(inspect.stack()[0][2:4][::-1])+': One or more modules missing')
    raise SystemExit

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~ Halt handler ~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

def handler(signum, frame):
    print 'Signal handler called with signal', signum
    print "CTRL-C pressed"
    sys.exit(0)

signal.signal(signal.SIGTERM, handler)
signal.signal(signal.SIGINT, handler)


'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~ Exec Function in verbose mode ~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

def verbose(function, *args):
    '''
    Run a function in verbose mode
    '''
    try:
        out = function(*args)
    except:
        # Something happened while runing function
        raise
        if DEBUG==True:
            print(str(inspect.stack()[0][2:4][::-1])+' Error')
            raise
    else:
        return(out)


'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~ Help message ~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

class PlatformHelp():
    def __init__(self):
        self.make_title()
        self.make_welcome()
        self.make_requisites()
        self.make_options()
    
    def make_title(self):
        nameversion = 10*'#'+3*' '+__shortname__+' v'+__version__+3*' '+10*'#'
        self.separator = len(nameversion)*'#'
        self.title = nameversion
    
    def make_welcome(self):
        self.welcome = 'Welcome to '+__shortname__+' ('+__longname__+')\n'
    
    def make_requisites(self):
        self.requisites = \
            'Requisites: Python 2.7; Scipy; Numpy;\n'+\
            '  Matplotlib; Pyfits; Pyephem\n\n'
    
    def make_options(self):
        self.options = \
            '-h: print this help message\n\n'+\
            '-b base_dir: \n'+\
            '  Use alternative base dir for the input/output\n'+\
            '-c config_file: \n'+\
            '  Use alternative config file\n'+\
            '-lowfreq file: Use an existing FlatField for the low freqs\n'+\
            '-highfreq file: Use an existing FlatField for the high freqs\n'+\
            '-radialprofile file: Use a CSV table with angular response as\n'+\
            '  the low freq component\n'+\
            '--reference file: Reference light image to determine the shape\n'+\
            '  of the synthetic FlatField if only radialprofile is given\n'+\
            '\n'
    
    def show_help(self):
        print(\
         self.separator+'\n'+self.title+'\n'+self.separator+'\n'+self.welcome+\
         self.requisites+self.options+self.separator)
        sys.exit(0)
        
    def incorrect_parameter(self,parameter):
        print('ERROR. Incorrect parameter: '+str(parameter))
        self.show_help()
    
    def date_or_file_input_error(self):
        print('ERROR. Date or file input')
        self.show_help()
    
    def no_parameters_error(self):
        print('ERROR. No input parameters especified')
        self.show_help()


'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~ Config file loading ~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

def get_config_filename(InputOptions):
    config_file = None
    try:
        assert(InputOptions.configfile!=False)
    except:
        print(str(inspect.stack()[0][2:4][::-1])+\
        ': config file not specified, use the default one:')
    else:
        config_file = InputOptions.configfile

    return(config_file)



class ReadOptions():
    def __init__(self,input_options):
        '''
        Read user input.
        Returns an Object with all user input parameters
        '''
        self.show_help = False
        self.process_input(input_options)
        self.test_flatfield_possible_builds()
        
        if self.show_help == True:
            Help = PlatformHelp()
            Help.show_help()
        
    
    def process_input(self,input_options):
        '''
        Process all input options
        '''
        
        # Lambda: http://docs.python.org/release/2.5.2/tut/node6.html [4.7.5]
        self.options = {\
         '-h': lambda: 'show_help',\
         '-c': lambda: 'use_configfile',\
         '-b': lambda: 'use_basedir',\
         '-lowfreq': lambda: 'lowfreq',\
         '-highfreq': lambda: 'highfreq',\
         '-radialprofile': lambda: 'radialprofile',\
         '--reference': lambda: 'file_reference'\
         }
        
        # By default, we wont show on screen nor save on disk.
        self.configfile = False;
        
        print('Input Options: '+str(input_options))
        
        self.input_options = input_options
        try: self.input_options[1]
        except Exception as e: 
            print(str(inspect.stack()[0][2:4][::-1])+'ERR. No imput options')
            self.no_parameters()
        else:
            while len(self.input_options)>1:
                input_option = self.options.get(self.input_options[1], lambda : None)()
                if input_option == 'show_help':
                    self.show_help = True
                    # Stop reading options. Program will halt
                    self.input_options = []
                elif input_option == 'use_configfile':
                    self.configfile = self.reference_file()
                elif input_option == 'use_basedir':
                    self.base_path = self.reference_file()
                elif input_option == 'lowfreq':
                    self.lowfreq_path = self.reference_file()
                elif input_option == 'highfreq':
                    self.highfreq_path = self.reference_file()
                elif input_option == 'file_reference':
                    self.file_reference = self.reference_file()
                elif input_option == 'radialprofile':
                    self.radialprofile = self.reference_file()
                else:
                    self.input_options.pop(1)
                
    def no_parameters(self):
        print('\nERR: Need more than one parameter')
        self.input_options = []
        self.show_help = True
    
    def reference_file(self):
        print('Path specified with '+self.input_options[1]+'. Extracting path')
        file_reference = None
        try: self.input_options[2]
        except:
            self.input_options.remove(self.input_options[1])
        else:
            if self.options.get(self.input_options[2], lambda : None)():
                self.input_options.remove(self.input_options[1])
            else:
                file_reference=self.input_options[2]
                self.input_options.remove(self.input_options[2]) 
                self.input_options.remove(self.input_options[1])
                return(file_reference)
    
    def not_enough_data(self):
        print(\
         '\n'+\
         'ERR: Not enough data given to build a synthetic FlatField\n'+\
         '     Need at least the low freq component or a file with\n'+\
         '     the angular response of the fisheye'\
         )
        
        self.input_options = []
        self.show_help = True
    
    def test_flatfield_possible_builds(self):
        '''
        Test what type of synthetic flatfield can be build from
        the user input data. We will try preferably to build
        a full Low+High Freq high quality synthetic FlatField.
        
        Returns the best synthetic FlatField that can be build (if any)
        as .build_type property.
        '''
        
        self.build_type=None
        
        def test(param):
            try: vars(self)[param]
            except: return(False)
            else: return(True)
        
        if test("lowfreq_path")== True and test("highfreq_path")==True:
            self.file_reference = self.lowfreq_path
            self.build_type = 'LowHighFreq'
        elif test("radialprofile")==True and test("highfreq_path")==True:
            self.file_reference = self.highfreq_path
            self.build_type = 'RadialHighFreq'
        elif test("lowfreq_path")== True:
            self.file_reference = self.lowfreq_path
            self.build_type = 'OnlyLowFreq'     
        elif test("radialprofile")==True and test("file_reference")==True:
            self.build_type = 'OnlyRadial'
        
        if self.build_type == None:
            self.not_enough_data()
        

def load_config_file(config_file):
    ''' 
    Open the config file
    This will set-up the observatory properties 
    such as Latitude and Longitude of the Observatory
    '''
    ConfigOptions_ = ConfigOptions(config_file)
    ImageInfoCommon = ImageInfo()
    
    class FakeInputOptions:
        void = None
    
    ImageInfoCommon.config_processing_common(ConfigOptions_,FakeInputOptions)
    return(ImageInfoCommon)


'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~ UTF-8 to latin1 conv. ~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

def encode_utf8_to_iso88591(utf8_text):
    '''
    Encode and return the given UTF-8 text as ISO-8859-1 (latin1) with
    unsupported characters replaced by '?', except for common special
    characters like smart quotes and symbols that we handle as well as we can.
    For example, the copyright symbol => '(c)' etc.

    If the given value is not a string it is returned unchanged.

    References:
    en.wikipedia.org/wiki/Quotation_mark_glyphs#Quotation_marks_in_Unicode
    en.wikipedia.org/wiki/Copyright_symbol
    en.wikipedia.org/wiki/Registered_trademark_symbol
    en.wikipedia.org/wiki/Sound_recording_copyright_symbol
    en.wikipedia.org/wiki/Service_mark_symbol
    en.wikipedia.org/wiki/Trademark_symbol
    '''
    if not isinstance(utf8_text, basestring):
        return utf8_text
    # Replace "smart" and other single-quote like things
    utf8_text = re.sub(
        u'[\u02bc\u2018\u2019\u201a\u201b\u2039\u203a\u300c\u300d]',
        "'", utf8_text)
    # Replace "smart" and other double-quote like things
    utf8_text = re.sub(
        u'[\u00ab\u00bb\u201c\u201d\u201e\u201f\u300e\u300f]',
        '"', utf8_text)
    # Replace copyright symbol
    utf8_text = re.sub(u'[\u00a9\u24b8\u24d2]', '(c)', utf8_text)
    # Replace registered trademark symbol
    utf8_text = re.sub(u'[\u00ae\u24c7]', '(r)', utf8_text)
    # Replace sound recording copyright symbol
    utf8_text = re.sub(u'[\u2117\u24c5\u24df]', '(p)', utf8_text)
    # Replace service mark symbol
    utf8_text = re.sub(u'[\u2120]', '(sm)', utf8_text)
    # Replace trademark symbol
    utf8_text = re.sub(u'[\u2122]', '(tm)', utf8_text)
    # Replace/clobber any remaining UTF-8 characters that aren't in ISO-8859-1
    return utf8_text.encode('ISO-8859-1', 'replace')


'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~ FlatField generation ~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

'''def load_image(image):
    print('Loading Data and Header for the given Image ...'),
    try:
        Image_HDU    = pyfits.open(image)
        Image_Data   = Image_HDU[0].data
        Image_Header = Image_HDU[0].header
        # Convert to string
        #Image_Header_text = Image_Header.tostring()
        #Image_Header_text = encode_utf8_to_iso88591(Image_Header_text)
        #Image_Header.fromstring(Image_Header_text)
        #print(Image_Header)
    except:
        print(inspect.stack()[0][2:4][::-1])
        raise
    else:
        print('OK')
        return(Image_Data,Image_Header) 
'''

def low_pass_filter(image,sigma=5):
    ''' Returns filtered image '''
    return(scipy.ndimage.filters.median_filter(image,sigma))

def high_pass_filter(image,sigma=5):
    ''' Returns filtered image '''
    return(image - low_pass_filter(image,sigma))

def angular_response(radialprofile,interporder=3,smoothing=None):
    ''' 
    Process a table (comma separated values) with angular
    response of the lens.
    Returns a tck tuple (knots,Bspline coef and degree). 
      See scipy.interpolate.splrep for further details.
    '''
    
    separation = ','
    endline = ''
    
    the_file = open(radialprofile)
    raw_data = the_file.readlines()
    
    data = raw_data
    
    # Remove trailing characters
    data = [item.replace(endline+'\n','') for item in data]
    # Remove spaces
    data = [item.replace(' ','') for item in data]
    # Remove empty lines
    data = [item for item in data if item!='']
    # Remove comments
    data = [item for item in data if item[0]!='#']
    
    # Separate angle and value
    data = np.array([item.split(separation) \
     for item in data],dtype='float')
    
    # Interpolation (spline)
    tck = scipy.interpolate.splrep(data[:,0],data[:,1],k=interporder,s=smoothing)
    return(tck)

def create_synthetic_image(radialprofile,ImageInfo):
    tck = angular_response(radialprofile,interporder=3,smoothing=None)
    IC = ImageCoordinates(ImageInfo)
    sx,sy = np.shape(IC.altitude_map)
    map_response = scipy.interpolate.splev(90.0-IC.altitude_map.flatten(),tck)
    map_response = np.array(map_response).reshape(sx,sy)
    map_response[np.isnan(map_response)] = 1
    return(map_response)

class FlatField():
    '''
    Synthetic FlatField generator
    '''

    def flat_from_low_and_highfreq(self,lowfreq_path,highfreq_path):
        # Load Low and High freq flats
        FlatLowFreq_Data,FlatLowFreq_Header = \
         verbose(load_image,lowfreq_path)
        FlatHighFreq_Data,FlatHighFreq_Header = \
         verbose(load_image,highfreq_path)
        
        # Check if the size match
        assert np.shape(FlatLowFreq_Data)==np.shape(FlatHighFreq_Data)
        
        # Normalization of Flats
        FlatHighFreq_Data = FlatHighFreq_Data/np.mean(FlatHighFreq_Data)
        FlatLowFreq_Data = FlatLowFreq_Data/np.mean(FlatLowFreq_Data)
        
        # Fourier filters
        print('Applying low and high pass filters ...')
        # Low pass
        FlatLow_filt = low_pass_filter(FlatLowFreq_Data)  
        # High pass
        FlatHigh_filt = high_pass_filter(FlatLowFreq_Data)
        
        # Build the synthetic flat field
        print('Creating synthetic Flatfield ...')
        self.FlatField = pyfits.HDUList(hdus=\
         [pyfits.PrimaryHDU(\
           data = FlatLow_filt+FlatHigh_filt,
           header = FlatLowFreq_Header)])
        self.FlatField[0].header.add_comment('Synthetic FF [PyASB]. Low + High')
        
        self.prefixname = 'Synthetic_FlatField_'+str(FlatLowFreq_Header['FILTER'])
    
    def flat_from_lowfreq(self,lowfreq_path):
        # Load Low and High freq flats
        FlatLowFreq_Data,FlatLowFreq_Header = \
         verbose(load_image,lowfreq_path)
        
        # Normalization of Flats
        FlatLowFreq_Data = FlatLowFreq_Data/np.mean(FlatLowFreq_Data)
        
        # Fourier filters
        print('Applying low and high pass filters ...')
        # Low pass
        FlatLow_filt = low_pass_filter(FlatLowFreq_Data)
        
        # Build the synthetic flat field
        print('Creating synthetic Flatfield ...')
        self.FlatField = pyfits.HDUList(hdus=\
         [pyfits.PrimaryHDU(\
           data = FlatLow_filt, header = FlatHighFreq_Header)])
        self.FlatField[0].header.add_comment('Synthetic FF [PyASB]. Low freqs')
        
        self.prefixname = 'Synthetic_FlatField_any'
    
    def flat_from_radial_and_highfreq(self,radialprofile,hightfreq_flat,ImageInfo):     
        FlatHighFreq_Data,FlatHighFreq_Header = \
         verbose(load_image,highfreq_path)
        
        print('Loading LowFreq from radial profile ...')
        ImageInfo.read_header(FlatHighFreq_Header)
        # Do not offset the center
        ImageInfo.delta_x = 0
        ImageInfo.delta_y = 0
        map_response = create_synthetic_image(radialprofile,ImageInfo)
        
        # Check if the size match
        assert np.shape(map_response)==np.shape(FlatHighFreq_Data)
        
        # Normalization of Flats
        FlatHighFreq_Data = FlatHighFreq_Data/np.mean(FlatHighFreq_Data)
        
        # Fourier filters
        print('Applying low and high pass filters ...')
        # High pass
        FlatHigh_filt = high_pass_filter(FlatLowFreq_Data)
        
        # Build the synthetic flat field
        print('Creating synthetic Flatfield ...')
        self.FlatField = pyfits.HDUList(hdus=\
         [pyfits.PrimaryHDU(\
           data = map_response+FlatHigh_filt,
           header = FlatHighFreq_Header)])
        self.FlatField[0].header.add_comment('Synthetic FF [PyASB]. Radial profile + High')
        
        self.prefixname = 'Synthetic_FlatField_'+str(FlatHighFreq_Header['FILTER'])
    
    def flat_from_radial(self,radialprofile,file_reference,ImageInfo):
        Header = verbose(load_image,file_reference)[1]
        
        # Tune a bit the Reference Header
        Header['FILTER'] = 'Johnson_common'
    
        print('Loading LowFreq from radial profile ...')
        ImageInfo.read_header(Header)
        # Do not offset the center
        ImageInfo.delta_x = 0
        ImageInfo.delta_y = 0
        map_response = create_synthetic_image(radialprofile,ImageInfo)
        
        # Build the synthetic flat field
        print('Creating synthetic Flatfield ...')
        self.FlatField = pyfits.HDUList(hdus=\
         [pyfits.PrimaryHDU(\
           data = map_response,
           header = Header)])
        self.FlatField[0].header.add_comment('Synthetic FF [PyASB]. Radial profile')
        
        self.prefixname = 'Synthetic_FlatField_any'
        
    
    def save_generated_flatfield(self,path_to_save):
        file_to_save = path_to_save+'/'+self.prefixname +'.fits'
        print('Save file to '+str(file_to_save))
        self.FlatField.writeto(file_to_save,clobber=True)
        
    
if __name__ == '__main__':
    # Read config and Observatory properties
    InputOptions = ReadOptions(sys.argv)
    config_file  = get_config_filename(InputOptions)
    ImageInfo    = load_config_file(config_file)

    print('Synthetic Flat Field generator ...')
    print('Generating ['+str(InputOptions.build_type)+'] flatfield')
    FlatField_ = FlatField()
    
    if InputOptions.build_type == 'LowHighFreq':
        FlatField_.flat_from_low_and_highfreq(\
         InputOptions.lowfreq_path,InputOptions.highfreq_path)
    elif InputOptions.build_type == 'RadialHighFreq':
        FlatField_.flat_from_radial_and_highfreq(\
         InputOptions.radialprofile,InputOptions.highfreq_path,ImageInfo)
    elif InputOptions.build_type == 'OnlyLowFreq':
        FlatField_.flat_from_lowfreq(InputOptions.lowfreq_path)
    elif InputOptions.build_type == 'OnlyRadial':
        FlatField_.flat_from_radial(\
         InputOptions.radialprofile,\
         InputOptions.file_reference,\
         ImageInfo)
    
    try: path_to_save = InputOptions.base_dir
    except:
        path_to_save = os.path.split(\
         os.path.abspath(InputOptions.file_reference))[0]
    
    FlatField_.save_generated_flatfield(path_to_save)
    
    

    
