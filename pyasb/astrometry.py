#!/usr/bin/env python

'''
PyASB astrometry functions.

Convert from one coordinate system to another.
____________________________

This module is part of the PyASB project, 
created and maintained by Mireia Nievas [UCM].
____________________________
'''

try:
    import sys,os,inspect
    import numpy as np
    import math
    #from numpy.math import pi,sin,cos,sqrt,atan2,asin
    import ephem
except:
    print(str(inspect.stack()[0][2:4][::-1])+': One or more modules missing')
    raise SystemExit

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


''' Astrometry with Pyephem '''

# Setup Pyephem Observatory

def pyephem_setup_common(ImageInfo):
    ObsPyephem = ephem.Observer()
    ObsPyephem.pressure = 0 # Dont consider atmospheric effects
    ObsPyephem.date = ImageInfo.date_string
    return ObsPyephem

def pyephem_setup_image(ImageInfo):
    ObsPyephem = pyephem_setup_common(ImageInfo)
    ObsPyephem.lat = (ImageInfo.latitude-ImageInfo.latitude_offset)*np.pi/180
    ObsPyephem.lon = (ImageInfo.longitude-ImageInfo.longitude_offset)*np.pi/180
    return ObsPyephem

def pyephem_setup_real(ImageInfo):
    ObsPyephem = pyephem_setup_common(ImageInfo)
    ObsPyephem.lat = (ImageInfo.latitude)*np.pi/180
    ObsPyephem.lon = (ImageInfo.longitude)*np.pi/180
    return ObsPyephem

''' 
Standalone functions. 
To be used on single points
'''

def horiz2xy_old(azimuth,altitude,ImageInfo):
    '''
    Return X,Y position in the image from azimuth/altitude horizontal coord.
    azimuth and altitude must be in degrees.
    '''
    
    Rfactor = ImageInfo.radial_factor*(180.0/np.pi)*np.sqrt(2*(1-np.sin(altitude*np.pi/180.0)))
    X = ImageInfo.resolution[0]/2 + ImageInfo.delta_x -\
        Rfactor*np.cos(azimuth*np.pi/180.0-ImageInfo.azimuth_zeropoint*np.pi/180.0)
    Y = ImageInfo.resolution[1]/2 + ImageInfo.delta_y +\
        Rfactor*np.sin(azimuth*np.pi/180.0-ImageInfo.azimuth_zeropoint*np.pi/180.0)
    return(X,Y)

def horiz2xy(azimuth,altitude,ImageInfo,derotate=True):
    '''
    Return X,Y position in the image from azimuth/altitude horizontal coord.
    azimuth and altitude must be in degrees.
    '''
    
    if derotate==True and (ImageInfo.latitude_offset!=0 or ImageInfo.longitude_offset!=0):
        # We have a real azimuth and altitude coordinates. If the camera is not
        # pointing to the zenith, we need to derotate the image.
        ra_appa,dec_appa = horiz2eq(\
        azimuth,altitude,\
        ImageInfo,\
        lat = ImageInfo.latitude,\
        lon = ImageInfo.longitude)
        
        azimuth,altitude = eq2horiz(\
         ra_appa,dec_appa,\
         ImageInfo,\
         lat = ImageInfo.latitude-ImageInfo.latitude_offset,\
         lon = ImageInfo.longitude-ImageInfo.longitude_offset)
    
    ### allow for different projections
    if ImageInfo.projection == 'ZEA':
        Rfactor = ImageInfo.radial_factor*(180.0/np.pi)*np.sqrt(2*(1-np.sin(altitude*np.pi/180.0)))
    elif ImageInfo.projection == 'ARC':
        Rfactor = ImageInfo.radial_factor*(180.0/np.pi)*(1-altitude/90.0)
    
    #X = ImageInfo.resolution[0]/2 + ImageInfo.delta_x -\
    #    Rfactor*np.cos(azimuth*np.pi/180.0-ImageInfo.azimuth_zeropoint*np.pi/180.0)
    
    X = ImageInfo.resolution[0]/2 - ImageInfo.delta_x +\
        Rfactor*np.cos(azimuth*np.pi/180.0-ImageInfo.azimuth_zeropoint*np.pi/180.0)
    Y = ImageInfo.resolution[1]/2 + ImageInfo.delta_y +\
        Rfactor*np.sin(azimuth*np.pi/180.0-ImageInfo.azimuth_zeropoint*np.pi/180.0)
    return(X,Y)

def xy2horiz(X,Y,ImageInfo,derotate=True):
    '''
    Return horizontal coordinates from X,Y position in the image.
    azimuth and altitude are in degrees.
    '''
    
    X = X - ImageInfo.resolution[0]/2.-ImageInfo.delta_x
    Y = Y - ImageInfo.resolution[1]/2.-ImageInfo.delta_y
    X = ImageInfo.resolution[0] - X # flip the image horizontally
    Rfactor = np.sqrt(X**2 + Y**2)/ImageInfo.radial_factor
    
    if np.size(Rfactor)>1:
        Rfactor[Rfactor>360./np.pi]=360./np.pi
    
    ### allow for different projections
    if ImageInfo.projection == 'ZEA':
        alt_factor = np.array(1-0.5*(np.pi*Rfactor/180.0)**2)
        altitude = (180.0/np.pi)*np.arcsin(alt_factor)
    
    elif ImageInfo.projection == 'ARC':
        altitude = 90*(1-(Rfactor/ImageInfo.radial_factor)*(np.pi/180.0))
    
    azimuth  = (360+ImageInfo.azimuth_zeropoint + 180.0*np.arctan2(Y,-X)/np.pi)%360
    
    if derotate==True and (ImageInfo.latitude_offset!=0 or ImageInfo.longitude_offset!=0):
        # We have a real azimuth and altitude coordinates. If the camera is not
        # pointing to the zenith, we need to rotate the image.
        
        ra_real,dec_real = horiz2eq(\
         azimuth,altitude,\
         ImageInfo,\
         lat = ImageInfo.latitude-ImageInfo.latitude_offset,\
         lon = ImageInfo.longitude-ImageInfo.longitude_offset)          
        
        azimuth,altitude = eq2horiz(\
         ra_real,dec_real,\
         ImageInfo,\
         lat = ImageInfo.latitude,\
         lon = ImageInfo.longitude)
    
    return(azimuth,altitude)

def eq2horiz(ra,dec,ImageInfo=None,sidtime=None,lat=None,lon=None):
    '''
    Calculate horizontal coordinates for the given observation site
    and the given point in the sky.
    The coordinates must be given in degrees or hours
    '''
        
    if lat==None: lat=ImageInfo.latitude
    if lon==None: lon=ImageInfo.longitude
    if sidtime==None: sidtime = ImageInfo.sidereal_time
    
    # Sidereal Time to Local Sidereal Time
    sidtime = sidtime + lon/15.
    
    lat = lat*np.pi/180.
    lon = lon*np.pi/180.
    sidtime = sidtime*np.pi/12.
    ra = ra*np.pi/12.
    dec = dec*np.pi/180.
    
    H = sidtime - ra

    _sina = np.sin(dec)*np.sin(lat)+np.cos(dec)*np.cos(lat)*np.cos(H)
    alt = np.arcsin(_sina)
    _cosa = np.cos(alt)
    _sinA = -np.sin(H)*np.cos(dec)/_cosa
    _cosA = (np.sin(dec)-np.sin(lat)*_sina)/(_cosa*np.cos(lat))
    az = np.arctan2(_sinA,_cosA)
    
    az = (az*180./np.pi)%360
    alt = alt*180./np.pi
    
    return(az,alt)
    
def horiz2eq(az,alt,ImageInfo=None,sidtime=None,lat=None,lon=None):
    '''
    Calculate equatorial coordinates for the given observation site
    and the given point in the sky
    The coordinates must be given in degrees or hours
    '''
    
    if lat==None: lat=ImageInfo.latitude
    if lon==None: lon=ImageInfo.longitude
    if sidtime==None: sidtime = ImageInfo.sidereal_time
    
    # Sidereal Time to Local Sidereal Time
    sidtime = sidtime + lon/15.
    
    lat = lat*np.pi/180.
    lon = lon*np.pi/180.
    sidtime = sidtime*np.pi/12.
    az = az*np.pi/180.
    alt = alt*np.pi/180.
    
    _sindec = np.sin(alt)*np.sin(lat)+np.cos(alt)*np.cos(lat)*np.cos(az)
    dec = np.arcsin(_sindec)
    _cosdec = np.cos(dec)
    _sinH = -np.sin(az)*np.cos(alt)/_cosdec
    _cosH = (np.sin(alt)-_sindec*np.sin(lat))/(_cosdec*np.cos(lat))
    
    H = np.arctan2(_sinH,_cosH)
    ra = sidtime - H
    
    ra = (ra*12./np.pi)%24
    dec = dec*180./np.pi
    
    return(ra,dec)

def zenith_position(ImageInfo):
    # Return X,Y position of zenith in the image.
    return horiz2xy(0,90,ImageInfo)

def optical_axis(ImageInfo):
    # Return horizontal coordinates of the optical axis
    return xy2horiz(ImageInfo.resolution[0]/2,ImageInfo.resolution[1]/2,ImageInfo)

def atmospheric_refraction(altitude,mode):
    # Return apparent (non-corrected from refraction) or 
    # real (corrected from refraction) altitude.
    # Garfinkel (1967), http://en.wikipedia.org/wiki/Atmospheric_refraction
    def cot(x):
        # Return cotangent of the given value
        return np.cos(x)/np.sin(x)
    
    if mode=='dir':
        # Return apparent altitude from the real one. 
        return altitude + (1.02/60)*cot(altitude + 10.3/(altitude+5.11))
    elif mode=='inv':
        # Return real altitude from the apparent one
        return altitude - (1.00/60)*cot(altitude + 7.31/(altitude+4.4))
    else:
        print 'Unknow mode '+mode+'. Cannot correct from atmospheric refraction.'
        return altitude

def calculate_airmass(altitude):
    # Estimate airmass from apparent altitude using Pickering (2002) model.
    # zdist in degrees
    return 1/np.sin((altitude+244./(165+47*altitude**1.1))*np.pi/180.)

'''
Vectorial functions.
Generate a class that contains a map of coordinates
that match the Image pixels
'''

class ImageCoordinates():
    def __init__(self,ImageInfo):
        self.calculate_altaz(ImageInfo)
    
    def calculate_altaz(self,ImageInfo):
        ''' Reimplementation with numpy arrays (fast on large arrays). 
            We need it as we will use a very large array'''
        
        # Image coordinates
        x = np.arange(ImageInfo.resolution[0])
        y = np.arange(ImageInfo.resolution[1])
        X,Y = np.meshgrid(x,y)
        
        # Unreal/projected altitude and azimuth
        az,alt = xy2horiz(X,Y,ImageInfo,derotate=False)
        self.azimuth_map = np.array(az,dtype='float16')
        self.altitude_map = np.array(alt,dtype='float16')
