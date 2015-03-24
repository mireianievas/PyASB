#!/usr/bin/env python

'''
PyASB launcher module

Concatenate processes
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
__status__ = "Prototype"  # "Prototype", "Development", or "Production"


import gc
import sys
import time

from pyasb.input_options import *
from pyasb.program_help import *
from pyasb.astrometry import *
from pyasb.star_calibration import *
from pyasb.load_fitsimage import *
from pyasb.bouguer_fit import *
from pyasb.sky_brightness import *
from pyasb.skymap_plot import *
from pyasb.write_summary import *


config_file = 'pyasb_config.cfg'

#@profile


class LoadImage(object):

    def __init__(self, InputOptions, ImageInfo, ConfigOptions, input_file=None):
        # Load Image file list
        if input_file == None:
            input_file = InputOptions.fits_filename_list[0]

        ''' Load fits image '''
        self.FitsImage = FitsImage(input_file)
        # Local copy of ImageInfo. We will process it further.
        self.ImageInfo = ImageInfo
        self.ImageInfo.read_header(self.FitsImage.fits_Header)
        self.ImageInfo.config_processing_specificfilter(ConfigOptions)

        try:
            self.FitsImage.reduce_science_frame(
                self.ImageInfo.darkframe,
                self.ImageInfo.sel_flatfield,
                MasterBias=None)
        except:
            raise
            print('Cannot reduce science frame')

        self.FitsImage.__clear__()
        self.output_paths(InputOptions)

    def output_paths(self, InputOptions):
        # Output file paths (NOTE: should be moved to another file or at least separated function)
        # Photometric table
        try:
            self.ImageInfo.photometry_table_path = InputOptions.photometry_table_path
        except:
            try:
                self.ImageInfo.photometry_table_path
            except:
                raise
                SystemExit

        # Star Map
        try:
            self.ImageInfo.skymap_path = InputOptions.skymap_path
        except:
            try:
                self.ImageInfo.skymap_path
            except:
                raise
                SystemExit

        # Bouguer Fit
        try:
            self.ImageInfo.bouguerfit_path = InputOptions.bouguerfit_path
        except:
            try:
                self.ImageInfo.bouguerfit_path
            except:
                raise
                SystemExit

        # SkyBrightness
        try:
            self.ImageInfo.skybrightness_map_path = InputOptions.skybrightness_map_path
        except:
            try:
                self.ImageInfo.skybrightness_map_path
            except:
                raise
                SystemExit

        try:
            self.ImageInfo.skybrightness_table_path = InputOptions.skybrightness_table_path
        except:
            try:
                self.ImageInfo.skybrightness_table_path
            except:
                raise
                SystemExit

        # Summary
        try:
            self.ImageInfo.summary_path = InputOptions.summary_path
        except:
            try:
                self.ImageInfo.summary_path
            except:
                raise
                SystemExit

#@profile


class ImageAnalysis(object):

    def __init__(self, Image):
        ''' Analize image and perform star astrometry & photometry. 
            Returns ImageInfo and StarCatalog'''

        ObsPyephem_ = pyephem_setup(Image.ImageInfo)

        Image.ImageInfo.catalog_filename = 'ducati_catalog.tsv'
        self.StarCatalog = StarCatalog(
            Image.FitsImage, Image.ImageInfo, ObsPyephem_)

        SkyMap_ = SkyMap(self.StarCatalog, Image.ImageInfo, Image.FitsImage)

#@profile


class MultipleImageAnalysis(object):

    def __init__(self, InputOptions):
        class StarCatalog_(object):
            StarList = []
            StarList_woPhot = []

        InputFileList = InputOptions.fits_filename_list

        for EachFile in InputFileList:
            EachImage = LoadImage(EachFile)
            EachAnalysis = ImageAnalysis(EachImage)
            self.StarCatalog.StarList.append(EachAnalysis.StarCatalog.StarList)
            self.StarCatalog.StarList_woPhot.append(
                EachAnalysis.StarCatalog.StarList_woPhot)

#@profile


class InstrumentCalibration(object):

    def __init__(self, ImageInfo, StarCatalog):
        self.BouguerFit = BouguerFit(ImageInfo, StarCatalog)
        self.BouguerFit.bouguer_plot(ImageInfo)

#@profile


class MeasureSkyBrightness(object):

    def __init__(self, FitsImage, ImageInfo, BouguerFit):
        ImageCoordinates_ = ImageCoordinates(ImageInfo)
        SkyBrightness_ = SkyBrightness(
            FitsImage, ImageInfo, ImageCoordinates_, BouguerFit)
        SkyBrightnessGraph_ = SkyBrightnessGraph(
            SkyBrightness_, ImageInfo, BouguerFit)
        self.SBzenith = SkyBrightness_.SBzenith
        self.SBzenith_err = SkyBrightness_.SBzenith_err

#@profile


def perform_complete_analysis(InputOptions, ImageInfoCommon, ConfigOptions, input_file):
    # Load Image into memory & reduce it.
        # Clean (no leaks)
    Image_ = LoadImage(
        InputOptions, ImageInfoCommon, ConfigOptions, input_file)

    # Look for stars that appears in the catalog, measure their fluxes. Generate starmap.
    # Clean (no leaks)
    ImageAnalysis_ = ImageAnalysis(Image_)

    # Calibrate instrument with image. Generate fit plot.
    # Clean (no leaks)
    InstrumentCalibration_ = InstrumentCalibration(
        Image_.ImageInfo,
        ImageAnalysis_.StarCatalog)

    # Measure sky brightness / background. Generate map.
    ImageSkyBrightness = MeasureSkyBrightness(
        Image_.FitsImage,
        Image_.ImageInfo,
        InstrumentCalibration_.BouguerFit)

    Summary_ = Summary(Image_, InputOptions, ImageAnalysis_,
                       InstrumentCalibration_, ImageSkyBrightness)

    gc.collect()
    # print(gc.garbage)


if __name__ == '__main__':
    gc.set_debug(gc.DEBUG_STATS)
    PlatformHelp_ = PlatformHelp()
    InputOptions = ReadOptions(sys.argv)
    ConfigOptions_ = ConfigOptions(config_file)
    ImageInfoCommon = ImageInfo()
    ImageInfo.config_processing_common(ImageInfoCommon, ConfigOptions_)
    try:
        assert(InputOptions.show_help == False)
    except:
        # Show help and halt
        PlatformHelp_.show_help()
        raise SystemExit

    for input_file in InputOptions.fits_filename_list:
        perform_complete_analysis(
            InputOptions, ImageInfoCommon, ConfigOptions_, input_file)

    gc.collect()

    d = dict()
    for o in gc.get_objects():
        name = type(o).__name__
        if name not in d:
            d[name] = 1
        else:
            d[name] += 1

    items = d.items()
    items.sort(key=lambda x: x[1])
    debug_file = open("debug_objects.txt", 'w')
    debug_file.close()
    debug_file = open("debug_objects.txt", 'a+')
    for key, value in items:
        print key, value
        debug_file.write(str(key) + ",\t" + str(value) + "\n")

    debug_file.close()
