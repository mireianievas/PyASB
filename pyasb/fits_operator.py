#!/usr/bin/env python

'''
FITS operations

Calculate the sum, difference, multiply or divide fits images.
The header is taken from the first image, with added comment.
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
    import signal
    import numpy as np
    import astropy.io.fits as pyfits
    import datetime
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


def load_image(image):
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

class FitsOperation():
    def __init__(self,fits1,fits2):
        self.loadfits(fits1,fits2)
    
    @staticmethod
    def get_datetime_filename(self):
        return(str(datetime.datetime.now()).replace(" ","_").replace(":","-").split(".")[0])
    
    def loadfits(self,fits1,fits2):
        self.Data1,self.Header1 = \
         verbose(load_image,fits1)
        self.Data2,self.Header2 = \
         verbose(load_image,fits2)

    def sumfits(self):
        self.DataResult   = self.Data1+self.Data2
        self.HeaderResult = self.Header1
    
    def subtractfits(self):
        self.DataResult   = self.Data1-self.Data2
        self.HeaderResult = self.Header1
    
    def multiplyfits(self):
        self.DataResult   = self.Data1*self.Data2
        self.HeaderResult = self.Header1
    
    def dividefits(self):
        self.DataResult   = self.Data1*1./self.Data2
        self.HeaderResult = self.Header1
    
    def normalizefits(self):
        self.DataResult   = self.DataResult*1./np.median(self.DataResult)
    
    def addnote(self,note="Fits edited with PyASB.fits_operator.py"):
        self.HeaderResult.add_comment(note)
    
    def savefits(self,filename=None):
        if filename==None: filename=self.get_datetime_filename
        pyfits.writeto(filename, self.DataResult, self.HeaderResult, clobber=True)
    
