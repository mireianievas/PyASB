#!/bin/env python2

import os,sys
import os.path
import astropy.io.fits as pyfits
import numpy as np
from PIL import Image
import datetime

def modification_date(filename):
    t  = os.path.getmtime(filename)
    dt =  datetime.datetime.fromtimestamp(t)
    return(dt.strftime("%Y%m%d_%H%M%S")

fullpath=sys.argv[1]
filepath=os.path.dirname(os.path.abspath(filename))
filename=fullpath.split("/")[-1]

img   = Image.open(fullpath)
array = np.array(img)

header = pyfits.Header()
header['SIMPLE']   = 'T'
header['BITPIX']   = '8'
header['NAXIS']    = '2'
header['NAXIS1']   = array.shape[1]
header['NAXIS2']   = array.shape[0]
header['BZERO']    = 0
header['BSCALE']   = 1
header['DATAMIN']  = 0
header['DATAMAX']  = 255
header['EXPOSURE'] = 1
header['FILTER']   = 'Johnson_V'
header['DATE']     = modification_date(fullpath)
header.add_historhy('This file was created with Python/PIL/numpy/astropy')
header.add_comment('Created by M. Nievas-Rosillo <mirph4k@gmail.com>') 

pyfits.writeto(fullpath+'.fits',data=array,header=header,overwrite=True)

