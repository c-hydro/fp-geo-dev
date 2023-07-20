#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
fp-geo - Make Land Maps

__date__ = '20230228'
__version__ = '2.0.1'
__author__ = 'Andrea Libertino (andrea.libertino@cimafoundation.org')
__library__ = 'fp-geo'

General command line:
python3 MakeLandMaps -settings_file "config/FILE.json"

Version(s):
20230228 (2.0.1)    --> Introduce editable grass setup
                    --> Add option to preserve original dem grid
20201203 (1.5.0)    --> Replaced alpha and beta routine with traditional one
                                     Changed mask routine (use shapefile mask)
                                     Added "-9999" external border - Hole DEM at endhoreic outlets
                                     Produce area km2 ancillary map
20201104 (1.1.5)    --> Add depressions support for forcing endhoreic basins - Produce macrobasins ancillary map
20201016 (1.1.1)    --> Changed alpha and beta routines - Bug fixes
20200922 (1.1.0)    --> Added json settings support
20200705 (1.0.2)    --> Fixed area in km calculation and latitude calculation
20200611 (1.0.1)    --> Fixed some export issues (null values of Int16 maps, non-masked layers LAT,LON,AREACELL) - Added logger
20200516 (1.0.0)    --> Changed the resample routine (approximated - does not require passage to UTM anymore)
			            Added the possibility to use an external mask for final clipping
			            Additionally export of the drainage directions in grass codification for r.water.outlet in the root folder
			            Various bux fixes
20200430 (0.0.5) -->    Beta release
"""

# -------------------------------------------------------------------------------------
# Complete library
import os

from grass_session import Session
import rasterio as rio
from copy import deepcopy
import numpy as np
import logging
from subprocess import PIPE
import math
import time
from datetime import date
from argparse import ArgumentParser
import json
import gdal

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information

alg_name = 'fp-geo tools - Make Land Maps'
alg_version = '2.0.1'
alg_release = '2023-02-28'
# Algorithm parameter(s)
time_format = '%Y%m%d%H%M'


# -------------------------------------------------------------------------------------
# Script main
def main():
    start_time = time.time()
    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings = get_args()

    # Set algorithm settings
    data_settings = read_file_json(alg_settings)

    # Check settings file
    if not "preserveGrid" in data_settings["algorithm"]["flags"].keys():
        logging.warning(' ---> WARNING! Preserve grid settings not set! Grid will not be preserved!')
        data_settings["algorithm"]["flags"]["preserveGrid"] = 0
    if not "grassbin" in data_settings["grass"].keys():
        logging.warning(' ---> WARNING! grassbin settings not set! Use standard value!')
        data_settings["grass"]["grassbin"] = "grass"

    # Configure grass
    os.environ['GISDBASE'] = data_settings["grass"]["grassdata"]
    os.environ['GRASSBIN'] = data_settings["grass"]["grassbin"]
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.pygrass.modules.shortcuts import vector as v
    from grass.pygrass.modules import Module
    from grass.script import parse_key_val, raster_info, pipe_command

    # Import useful settings
    InRes = data_settings['algorithm']['resolutions']['InResKm']
    OutRes = data_settings['algorithm']['resolutions']['OutResKm']
    dominio = data_settings['algorithm']['general']['domain']
    soglia_basins_originale = data_settings["algorithm"]["resolutions"]["soglia_basin_originale_ncell"]
    InResStr = str(InRes).replace(".", "_")
    OutResStr = str(OutRes).replace(".", "_")

    # Fill path
    path_settings = fillScriptSettings(data_settings, data_settings['algorithm']['general']['domain'])

    # Make useful paths
    os.makedirs(path_settings["outPath"], exist_ok=True)
    os.makedirs(path_settings["ancillaryPath"], exist_ok=True)
    os.makedirs(path_settings['logPath'], exist_ok=True)

    # Set algorithm logging
    set_logging(logger_file=os.path.join(path_settings['logPath'], data_settings['algorithm']['general']['logName']))

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Set algorithm logging and initialize algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME : ' + date.today().strftime("%d-%B-%Y %H:%m"))
    logging.info(' ==> START ... ')
    logging.info(' ==> ALGORITHM SETTINGS <== ')
    logging.info(' ----> Domain: ' + dominio)
    logging.info(' ----> Input resolution: ' + str(InRes) + ' km Output resolution: ' + str(OutRes) + ' km')
    logging.info(' ----> Stream threshold: ' + str(soglia_basins_originale * InRes * InRes) + ' km2')

    logging.info(' ')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Rescale to OutRes with cubic resample
    OutResDegree = kilometers2degrees(OutRes)

    if data_settings["algorithm"]["flags"]["changeDemRes"]:
        logging.info(' ----> Rescaling DEM to ' + str(OutRes) + ' km ... ')
        os.system('gdalwarp ' + os.path.join(path_settings["demPath"],
                                             data_settings["algorithm"]["general"]["demName"]) + ' ' + os.path.join(
            path_settings["ancillaryPath"], 'tempDem.tif') + ' -tr ' + str(OutResDegree) + ' ' + str(
            OutResDegree) + ' -r cubic -overwrite -q')
    else:
        OutRes = InRes

    # Calculate accumulation threshold on rescaled DEM
    soglia_basins_post = math.ceil((soglia_basins_originale * InRes * InRes) / (OutRes * OutRes))

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Depitting dem with saga GIS
    logging.info(' ---> Data loading to SAGA GIS')
    if data_settings["algorithm"]["flags"]["changeDemRes"]:
        os.system("gdal_translate -q -of SAGA " + os.path.join(path_settings["ancillaryPath"],
                                                               'tempDem.tif') + " " + os.path.join(
            path_settings['ancillaryPath'], dominio + "_dem.sdat"))
    else:
        os.system("gdal_translate -q -of SAGA " + os.path.join(path_settings["demPath"],
                                                               data_settings["algorithm"]["general"][
                                                                   "demName"]) + " " + os.path.join(
            path_settings['ancillaryPath'], data_settings['algorithm']['general']['domain'] + "_dem.sdat"))

    logging.info(' ----> SAGA: Identify sink routes')
    os.system("saga_cmd -f=s ta_preprocessor 1 -ELEVATION " + os.path.join(path_settings['ancillaryPath'],
                                                                           dominio + "_dem.sdat") + " -SINKROUTE " + os.path.join(
        path_settings['ancillaryPath'], dominio + "_SinkRoute.sdat"))
    logging.info(' ----> SAGA: Fill sinks')
    os.system("saga_cmd -f=s ta_preprocessor 2 -DEM " + os.path.join(path_settings['ancillaryPath'],
                                                                     dominio + "_dem.sdat") + " -SINKROUTE " + os.path.join(
        path_settings['ancillaryPath'], dominio + "_SinkRoute.sdat") + " -DEM_PREPROC " + os.path.join(
        path_settings['ancillaryPath'], dominio + "_DemProcessed.sdat") + " -METHOD 1")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Inputs preparation
    # Open an existing mapset or create a new one
    logging.info(' ----> Opening mapset: ' + dominio + '_land_data ... ')
    with Session(gisdb=data_settings["grass"]["grassdata"], location=data_settings["grass"]["location"],
                 mapset=dominio + "_land_data", create_opts=''):
        try:
            r.mask(flags='r')
        except:
            pass

        g.gisenv(set="GRASS_VERBOSE=-1")
        # Import grass modules
        r.in_gdal = Module('r.in.gdal')
        r_import = Module('r.import')
        v.in_ogr = Module('v.in.ogr')
        v.to_rast = Module('v.to.rast')
        r.to_vect = Module('r.to.vect')
        v.to_points = Module('v.to.points')
        r.stream_channel = Module('r.stream.channel')

        # Import DEM (original or rescaled one)
        logging.info(' ----> Importing DEM ... ')
        os.system('gdal_translate -of GTiff ' + os.path.join(path_settings['ancillaryPath'],
                                                             dominio + "_DemProcessed.sdat") + ' ' + os.path.join(
            path_settings['ancillaryPath'], dominio + "_DemProcessed.tif"))
        r_import(input=os.path.join(path_settings['ancillaryPath'], dominio + "_DemProcessed.sdat"),
                 output='dem.' + OutResStr + 'km', flags='o', overwrite=True, quiet=True)
        g.region(raster='dem.' + OutResStr + 'km', save='small_region', overwrite=True)
        logging.info(' ----> Setting null values ... ')
        r.null(map='dem.' + OutResStr + 'km', setnull=-1e+308 - -9999, quiet=True)

        # Start carving operations (if activated)
        if data_settings["algorithm"]["flags"]["carve"]:
            # If an external network is not provided for carving it is calculated with r.watershed on the original DEM
            if data_settings["algorithm"]["flags"]["existStream"]:
                try:
                    g.findfile(element='cell', file="stream." + InResStr + "km.th" + str(soglia_basins_originale))
                    logging.info(
                        ' ----> The stream at the resolution ' + InResStr + ' exist in the mapset. Loading ... ')
                except:
                    logging.info(
                        ' ----> The stream at the resolution ' + InResStr + ' does not exist in the mapset. Calculating ... ')
                    r.in_gdal(
                        input=os.path.join(path_settings["demPath"], data_settings["algorithm"]["general"]["demName"]),
                        output='dem.' + InResStr + 'km', overwrite=True, flags='o')
                    g.region(raster="dem." + InResStr + "km")
                    r.watershed(elevation="dem." + InResStr + "km", threshold=soglia_basins_originale,
                                basin="basin." + InResStr + "km.th" + str(soglia_basins_originale),
                                stream="stream." + InResStr + "km.th" + str(soglia_basins_originale), flags="sab",
                                overwrite=True)
                    r.to_vect = Module('r.to.vect')
                    r.to_vect(input="stream." + InResStr + "km.th" + str(soglia_basins_originale),
                              output="stream_" + InResStr + "km_th" + str(soglia_basins_originale), type='line',
                              overwrite=True, quiet=True)
                v.to_rast(input="stream_" + InResStr + "km_th" + str(soglia_basins_originale),
                          output="stream." + InResStr + "km.th" + str(soglia_basins_originale), type='line', use='cat',
                          flags='d', overwrite=True, quiet=True)
            # If an external stream is provided, it is imported
            else:
                logging.info(' ----> Stream for carving is provided. Loading ... ')
                v.in_ogr(input=os.path.join(path_settings["streamPath"],
                                            data_settings["algorithm"]["general"]["streamName"]),
                         output="stream_" + InResStr + "km_th" + str(soglia_basins_originale) + "_edit", flags='o',
                         overwrite=True, quiet=True)
                v.to_rast(input="stream_" + InResStr + "km_th" + str(soglia_basins_originale) + "_edit",
                          output="stream." + InResStr + "km.th" + str(soglia_basins_originale), type='line', use='cat',
                          flags='d', overwrite=True, quiet=True)

            g.region(raster='dem.' + OutResStr + 'km', save='small_region', overwrite=True)

            # Start carving routine
            logging.info(' ----> Carving operations ... ')
            # Neighbors MINIMUM near canals
            r.neighbors(input='dem.' + OutResStr + 'km',
                        selection='stream.' + InResStr + 'km.th' + str(soglia_basins_originale),
                        output='carver.nbmindem', method='minimum', size=3, overwrite=True, quiet=True)
            # Neighbour AVERAGE near canals
            r.neighbors(input='dem.' + OutResStr + 'km',
                        selection='stream.' + InResStr + 'km.th' + str(soglia_basins_originale),
                        output='carver.nbavedem', method='average', size=3, overwrite=True, quiet=True)
            # Higher dem value, nearby streams
            ret = Module('r.univar', flags='gr', map='carver.nbavedem', stdout_=PIPE)
            stats = parse_key_val(ret.outputs.stdout)
            maxDem = stats.max

            # Compute Dig
            r.mask(raster='stream.' + InResStr + 'km.th' + str(soglia_basins_originale), overwrite=True, quiet=True)
            # Convert elevation to degrees (proportion max vs 2pi)
            r.mapcalc('carver.nbmindem.zero = max(carver.nbmindem,0)', overwrite=True, quiet=True)
            # Convert zeros assigning min value (without negatives) or 1 if min is zero
            r.mapcalc('carver.nbavedem.zero = max(carver.nbavedem,abs(carver.nbmindem),1)', overwrite=True, quiet=True)
            r.mapcalc('carver.dig.proportion = ((2*3.14)*carver.nbmindem)/' + maxDem, overwrite=True, quiet=True)
            # Compute Numerator
            r.mapcalc('''carver.dig.numerator = (1-exp(-''' + str(
                data_settings["algorithm"]["carve_settings"]["LAMBDA"]) + '''*carver.dig.proportion))''',
                      overwrite=True, quiet=True)
            # Compute Denominator
            r.mapcalc('''carver.dig.denominator = (1-exp(-''' + str(
                data_settings["algorithm"]["carve_settings"]["LAMBDA"]) + '''*(2*3.14)))''', overwrite=True, quiet=True)
            # Compute normalised correction function
            r.mapcalc('carver.dig.normalised = carver.dig.numerator/carver.dig.denominator', overwrite=True, quiet=True)
            # Compute weighted function
            r.mapcalc('carver.dig.weighted = carver.dig.normalised*' + str(
                data_settings["algorithm"]["carve_settings"]["maxDig"]), overwrite=True, quiet=True)
            r.mask(flags='r')
            g.remove(type='raster', name='MASK', flags='f', quiet=True)

            ## General dig, alfa correction (dig more in flat areas)
            if not data_settings["algorithm"]["carve_settings"]["alfa"] == -9999:
                logging.info(' ----> Alpha is provided. Carving more flat areas ... ')
                r.mask(raster='stream.' + InResStr + 'km.th' + str(soglia_basins_originale), overwrite=True, quiet=True)
                r.mapcalc('''carver.alfa.corr = ''' + str(data_settings["algorithm"]["carve_settings"][
                                                              "alfa"]) + ''' + (carver.nbmindem.zero/abs(carver.nbavedem.zero))''',
                          overwrite=True, quiet=True)
                r.mapcalc('out.' + dominio + '.dig=carver.dig.weighted*carver.alfa.corr', overwrite=True, quiet=True)
                r.mask(flags='r')
                g.remove(type='raster', name='MASK', flags='f', quiet=True)
                g.remove(type='raster', name='carver.alfa.corr', flags='f', quiet=True)
            else:
                r.mapcalc('out.' + dominio + '.dig=carver.dig.weighted', overwrite=True, quiet=True)

            # Carve in meters
            r.mapcalc('''carver.canals = if(isnull(stream.''' + InResStr + '''km.th''' + str(
                soglia_basins_originale) + '''), null(), carver.nbmindem-out.''' + dominio + '''.dig)''',
                      overwrite=True, quiet=True)
            r.mapcalc(
                '''out.''' + dominio + '''.dem.''' + OutResStr + '''km = if(isnull(stream.''' + InResStr + '''km.th''' + str(
                    soglia_basins_originale) + '''), dem.''' + OutResStr + '''km, carver.canals)''', overwrite=True,
                quiet=True)

            # g.remove(type='raster', pattern='carver.*', flags='f', quiet=True)
        else:
            r.mapcalc('''out.''' + dominio + '''.dem.''' + OutResStr + '''km = dem.''' + OutResStr + '''km''',
                      overwrite=True, quiet=True)
        logging.info(' ----> Carving operations ... DONE')

        # Export unmasked non-walled DEM for rasterize the depression on the same grid
        r.out_gdal(input="out." + dominio + ".dem." + OutResStr + "km",
                   output=os.path.join(path_settings["ancillaryPath"], "dem_unwalled.tif"), format='GTiff',
                   type='Float32', overwrite=True, nodata=-9999, flags='fc', quiet=True)
        convertAIIGrid(os.path.join(path_settings["ancillaryPath"], "dem_unwalled.tif"),
                       os.path.join(path_settings["ancillaryPath"], dominio + ".dem_unwalled.txt"), 'Float32')

        # Condition the DEM with the walls shapefile
        logging.info(' ----> Conditioning the DEM with walls... ')
        if data_settings["algorithm"]["flags"]["existWall"]:
            v.in_ogr(input=os.path.join(path_settings["wallPath"], data_settings["algorithm"]["general"]["wallName"]),
                     output="wall_shp", flags='o', overwrite=True, quiet=True)
            v.to_rast(input='wall_shp', output='wall.' + OutResStr + 'km', type='line', use='val', value=1,
                      overwrite=True, quiet=True)
            r.null(map='wall.' + OutResStr + 'km', setnull=-9999, null=0)
            r.mapcalc('''out.''' + dominio + '''.dem.walled = if(wall.''' + OutResStr + '''km, ''' + str(
                data_settings["algorithm"]["wall_settings"][
                    "height"]) + '''*wall.''' + OutResStr + '''km, out.''' + dominio + '''.dem.''' + OutResStr + '''km)''',
                      overwrite=True, quiet=True)
        else:
            r.mapcalc('''out.''' + dominio + '''.dem.walled = out.''' + dominio + '''.dem.''' + OutResStr + '''km''',
                      overwrite=True, quiet=True)
            logging.info(' ----> DEM is ready ... ')

        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Derivatives from carved DEM
        if data_settings["algorithm"]["flags"]["existDepressions"]:
            logging.info(' ----> Importing real depressions ... ')
            logging.info(' ---> Rasterizing depressions ... ')

            # Extract raster grid for rasterizing the depression shapefile (v.to.rast actually not working with point data)
            raster = gdal.Open(os.path.join(path_settings["ancillaryPath"], dominio + ".dem_unwalled.txt"))
            gt = raster.GetGeoTransform()
            coordRst = [gt[0], gt[3] + raster.RasterYSize * gt[5], gt[0] + raster.RasterXSize * gt[1], gt[3]]
            os.system('gdal_rasterize -l ' + data_settings["algorithm"]["general"]["depressionsName"].replace('.shp',
                                                                                                              '') + ' -burn 2.0 -ts ' + str(
                abs(raster.RasterXSize)) + ' ' + str(abs(raster.RasterYSize)) + ' -a_nodata -9999.0 -te ' + str(
                coordRst[0]) + ' ' + str(coordRst[1]) + ' ' + str(coordRst[2]) + ' ' + str(
                coordRst[3]) + ' -ot Float32 -of GTiff ' + os.path.join(path_settings['depressionsPath'],
                                                                        data_settings["algorithm"]["general"][
                                                                            "depressionsName"]) + ' ' + os.path.join(
                path_settings['ancillaryPath'], 'depression_rst.tif'))
            convertAIIGrid(os.path.join(path_settings['ancillaryPath'], 'depression_rst.tif'),
                           os.path.join(path_settings["ancillaryPath"], dominio + ".depressions.txt"), 'Float32')
            # os.system('rm -r ' + os.path.join(path_settings['ancillaryPath'], 'depression_rst.tif'))

            r.in_gdal(input=os.path.join(path_settings["ancillaryPath"], dominio + ".depressions.txt"),
                      output='depressions.' + OutResStr + 'km', overwrite=True, quiet=True)

            logging.info(' ----> Computing the hydrological derivatives ... ')
            r.watershed(elevation="out." + dominio + ".dem.walled",
                        threshold=soglia_basins_post,
                        basin="basin." + OutResStr + "km.th" + str(soglia_basins_post),
                        stream="stream." + OutResStr + "km.th" + str(soglia_basins_post),
                        depression="depressions." + OutResStr + "km",
                        accumulation="TCA." + OutResStr + "km.th" + str(soglia_basins_post),
                        drainage="pnt_grass." + OutResStr + "km.th" + str(soglia_basins_post),
                        flags="sab", overwrite=True, quiet=True)

        else:

            logging.info(' ----> Computing the hydrological derivatives ... ')
            r.watershed(elevation="out." + dominio + ".dem.walled",
                        threshold=soglia_basins_post,
                        basin="basin." + OutResStr + "km.th" + str(soglia_basins_post),
                        stream="stream." + OutResStr + "km.th" + str(soglia_basins_post),
                        accumulation="TCA." + OutResStr + "km.th" + str(soglia_basins_post),
                        drainage="pnt_grass." + OutResStr + "km.th" + str(soglia_basins_post),
                        flags="sab", overwrite=True, quiet=True)

        logging.info(' ----> Exporting the macrobasins map ... ')
        r.stream_basins(direction="pnt_grass." + OutResStr + "km.th" + str(soglia_basins_post),
                        stream_rast="stream." + OutResStr + "km.th" + str(soglia_basins_post), flags='l',
                        basins='mbsn.' + OutResStr + 'km', overwrite=True, quiet=True)
        r.out_gdal(input='mbsn.' + OutResStr + 'km',
                   output=os.path.join(path_settings['ancillaryPath'], dominio + ".mbsn.tif"), format='GTiff',
                   overwrite=True, flags='fc', quiet=True)

        # If external basin mask exist set the mask
        logging.info(' ----> Setting the basin mask ... ')
        if data_settings["algorithm"]["flags"]["existBasinMask"]:
            v.in_ogr(input=os.path.join(path_settings["maskPath"], data_settings["algorithm"]["general"]["maskName"]),
                     output=dominio + '_mask_v', overwrite=True, quiet=True)
            v.to_rast(input=dominio + '_mask_v', output="basinM", use='val', type='area', value=1, overwrite=True,
                      quiet=True)
        else:
            r.mapcalc(
                'basinM=int((out.' + dominio + '.dem.' + OutResStr + 'km+1000)/(out.' + dominio + '.dem.' + OutResStr + 'km+1000))',
                overwrite=True, quiet=True)
            r.to_vect(input='basinM', output=dominio + '_mask_v', type='area', overwrite=True, quiet=True)

        r.mapcalc("basinMask=if(out." + dominio + ".dem.walled<-9000,0,basinM)", overwrite=True)
        r.mask(raster="basinMask", maskcats=1, overwrite=True, quiet=True)

        v.buffer(input=dominio + '_mask_v', output=dominio + '_mask_v_buffer', distance=OutResDegree * (3 / 4),
                 quiet=True, overwrite=True)
        if not data_settings["algorithm"]["flags"]["preserveGrid"] == 1:
            g.region(vector=dominio + '_mask_v_buffer', align="out." + dominio + ".dem." + OutResStr + "km",
                     res=OutResDegree, save='extended_region', overwrite=True)
        else:
            g.region(raster='dem.' + OutResStr + 'km', save='extended_region', overwrite=True)

        # Export non-carved DEM (add holes in depressions locations)
        logging.info(' ----> Exporting DEM ... ')
        if data_settings["algorithm"]["flags"]["existDepressions"]:
            r.mapcalc(
                "out." + dominio + ".dem." + OutResStr + "km_holes=if(isnull(depressions." + OutResStr + "km),out." + dominio + ".dem." + OutResStr + "km,null())",
                quiet=True, overwrite=True)
            r.out_gdal(input="out." + dominio + ".dem." + OutResStr + "km_holes",
                       output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff', type='Float32',
                       overwrite=True, nodata=-9999, flags='fc', quiet=True)
        else:
            r.out_gdal(input="out." + dominio + ".dem." + OutResStr + "km",
                       output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff', type='Float32',
                       overwrite=True, nodata=-9999, flags='fc', quiet=True)
        convertAIIGrid(os.path.join(path_settings["outPath"], "temp.tiff"),
                       os.path.join(path_settings["outPath"], dominio + ".dem.txt"), 'Float32')

        # Compute upstream area
        logging.info(' ----> Exporting upstream area ... ')
        r.out_gdal(input="TCA." + OutResStr + "km.th" + str(soglia_basins_post),
                   output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff', overwrite=True,
                   nodata=-9999, type='Int32', flags='fc', quiet=True)
        convertAIIGrid(os.path.join(path_settings["outPath"], "temp.tiff"),
                       os.path.join(path_settings["outPath"], dominio + ".area.txt"), 'Int32')

        # Compute choice (channel network) fixing null values
        logging.info(' ----> Computing and exporting choice ... ')
        r.mapcalc('''choice.''' + OutResStr + '''km.th''' + str(
            soglia_basins_post) + '''= if(isnull(stream.''' + OutResStr + '''km.th''' + str(
            soglia_basins_post) + '''), 0, 1)''', overwrite=True)
        r.out_gdal(input='choice.' + OutResStr + 'km.th' + str(soglia_basins_post),
                   output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff', overwrite=True,
                   type='Int16', flags='fc', quiet=True)
        convertAIIGrid(os.path.join(path_settings["outPath"], "temp.tiff"),
                       os.path.join(path_settings["outPath"], dominio + ".choice.txt"), 'Int16')

        # Compute hydrological pointers in HMC codification (save also the Drainage Direction map in Grass codification in the path root for r.water.outlet) fixing null values
        logging.info(' ----> Computing and exporting hydrological pointers ... ')
        writeGrass2HMC(path_settings['ancillaryPath'])
        r.reclass(input="pnt_grass." + OutResStr + "km.th" + str(soglia_basins_post),
                  output="pnt_HMC." + OutResStr + "km.th" + str(soglia_basins_post),
                  rules=os.path.join(path_settings['ancillaryPath'], "grass2cont.txt"), overwrite=True, quiet=True)
        r.out_gdal(input="pnt_grass." + OutResStr + "km.th" + str(soglia_basins_post),
                   output=os.path.join(path_settings['ancillaryPath'], "DrainageDirection.Grass.tif"), format='GTiff',
                   overwrite=True, type='Int16', flags='fc', quiet=True)
        r.out_gdal(input="pnt_HMC." + OutResStr + "km.th" + str(soglia_basins_post),
                   output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff', overwrite=True,
                   type='Int16', flags='fc', quiet=True)
        convertAIIGrid(os.path.join(path_settings["outPath"], "temp.tiff"),
                       os.path.join(path_settings["outPath"], dominio + ".pnt.txt"), 'Int16')

        # Compute and export mask fixing null values
        logging.info(' ----> Exporting mask ... ')
        r.mapcalc("mask = TCA." + OutResStr + "km.th" + str(soglia_basins_post) + "/TCA." + OutResStr + "km.th" + str(
            soglia_basins_post), overwrite=True, quiet=True)
        r.out_gdal(input="mask", output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff', flags='fc',
                   overwrite=True, type='Int16', quiet=True)
        convertAIIGrid(os.path.join(path_settings["outPath"], "temp.tiff"),
                       os.path.join(path_settings["outPath"], dominio + ".mask.txt"), 'Int16')

        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Geographical references layers
        # Remove mask for exporting non-masked maps
        try:
            r.mask(flags='r')
        except:
            pass

        # Compute Lat and Lon
        logging.info(' ----> Computing and exporting latitude and longitude ... ')
        r.mapcalc("latitude = y()", overwrite=True, quiet=True)
        r.mapcalc("longitude= x()", overwrite=True, quiet=True)
        r.out_gdal(input="latitude", output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff',
                   overwrite=True, type='Float32', flags='fc', quiet=True)
        convertAIIGrid(os.path.join(path_settings["outPath"], "temp.tiff"),
                       os.path.join(path_settings["outPath"], dominio + ".lat.txt"), 'Float32', precision=4)
        r.out_gdal(input="longitude", output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff',
                   overwrite=True, type='Float32', flags='fc', quiet=True)
        convertAIIGrid(os.path.join(path_settings["outPath"], "temp.tiff"),
                       os.path.join(path_settings["outPath"], dominio + ".lon.txt"), 'Float32', precision=4)

        # Compute areacell
        logging.info(' ----> Computing areacell ... ')
        r.mapcalc("areacell = area()", overwrite=True, quiet=True)
        r.out_gdal(input="areacell", output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff',
                   overwrite=True, type='Float32', flags='fc', quiet=True)
        convertAIIGrid(os.path.join(path_settings["outPath"], "temp.tiff"),
                       os.path.join(path_settings["outPath"], dominio + ".areacell.txt"), 'Float32')

        g.region(region='small_region')
        # Calculate area km2
        logging.info(' ----> Writing area km2 ancillary map ... ')
        r.mapcalc("areaCellkm = areacell/(10^6)", overwrite=True, quiet=True)
        r.accumulate(direction="pnt_grass." + OutResStr + "km.th" + str(soglia_basins_post), format='auto',
                     accumulation='temp.TCA.km2', weight='areaCellkm', overwrite=True, quiet=True)
        r.out_gdal(input='temp.TCA.km2', output=os.path.join(path_settings['ancillaryPath'], 'upstream_area_km2.tif'),
                   format='GTiff', type='Float32', overwrite=True, nodata=-9999, flags='fc', quiet=True)

        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Alpha and beta
        # Compute alpha
        logging.info(' ----> Computing alpha and beta... ')
        makeAlphaBeta(path_settings["outPath"], path_settings["ancillaryPath"], dominio)

        convertAIIGrid(os.path.join(path_settings["ancillaryPath"], "temp_alpha.tif"),
                       os.path.join(path_settings["outPath"], dominio + ".alpha.txt"), 'Float32', precision=6)
        convertAIIGrid(os.path.join(path_settings["ancillaryPath"], "temp_beta.tif"),
                       os.path.join(path_settings["outPath"], dominio + ".beta.txt"), 'Float32', precision=6)

        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Width and depth
        # Compute width&depth
        if data_settings["algorithm"]["flags"]["compute_WeD"]:

            logging.info(' ----> With and depth module is active ...')
            logging.info(' ----> Calculate width and depth with provided parameters ...')

            if data_settings["algorithm"]["flags"]["existBasinMask"]:
                r.mask(raster="basinMask", maskcats=1, overwrite=True, quiet=True)

            r.mapcalc("areaCellkm = areacell/(10^6)", overwrite=True, quiet=True)
            r.accumulate(direction="pnt_grass." + OutResStr + "km.th" + str(soglia_basins_post), format='auto',
                         accumulation='temp.TCA.km2', weight='areaCellkm', overwrite=True)
            r.mapcalc("width." + OutResStr + "km.th" + str(soglia_basins_post) + " = " + str(
                data_settings["algorithm"]["width_depth_settings"]["aAcc2Width"]) + " * (temp.TCA.km2 ^ " + str(
                data_settings["algorithm"]["width_depth_settings"]["bAcc2Width"]) + ")", overwrite=True, quiet=True)
            r.mapcalc("temp.depth = " + str(data_settings["algorithm"]["width_depth_settings"][
                                                "aWidth2Depth"]) + " * (width." + OutResStr + "km.th" + str(
                soglia_basins_post) + " ^ " + str(
                data_settings["algorithm"]["width_depth_settings"]["bWidth2Depth"]) + ")", overwrite=True, quiet=True)
            r.mapcalc(
                "depth." + OutResStr + "km.th" + str(soglia_basins_post) + " = if(choice." + OutResStr + "km.th" + str(
                    soglia_basins_post) + " ==1, temp.depth, null())", overwrite=True, quiet=True)

            g.remove(type='raster', pattern='temp.*', flags='f', quiet=True)

            g.region(region='extended_region')
            logging.info(' ----> Exporting width and depth ...')
            r.out_gdal(input="width." + OutResStr + "km.th" + str(soglia_basins_post),
                       output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff', nodata=-9999,
                       overwrite=True, type='Float32', flags='f', quiet=True)
            convertAIIGrid(os.path.join(path_settings["outPath"], "temp.tiff"),
                           os.path.join(path_settings["outPath"], dominio + ".width.txt"), 'Float32')

            r.out_gdal(input="depth." + OutResStr + "km.th" + str(soglia_basins_post),
                       output=os.path.join(path_settings["outPath"], "temp.tiff"), format='GTiff', nodata=-9999,
                       overwrite=True, type='Float32', flags='f', quiet=True)
            convertAIIGrid(os.path.join(path_settings["outPath"], "temp.tiff"),
                           os.path.join(path_settings["outPath"], dominio + ".lfl.txt"), 'Float32')
            os.system("cp " + os.path.join(path_settings["outPath"], dominio + ".lfl.txt") + " " + os.path.join(
                path_settings["outPath"], dominio + ".rfl.txt"))

        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Cleaning and closing

        logging.info(' ----> Cleaning temporary data ...')

        # Clean geographic reference data (useless for Continuum)
        os.system("rm " + os.path.join(path_settings["outPath"], "*.prj"))
        os.system("rm " + os.path.join(path_settings["outPath"], "*.xml"))

        try:
            r.mask(flags='r')
        except:
            pass

    if data_settings["algorithm"]["flags"]["temporaryMapset"]:
        logging.info(' ----> Cleaning temporary mapset ...')
        os.system('rm -r ' + os.path.join(data_settings["grass"]["grassdata"], data_settings["grass"]["location"],
                                          dominio + "_land_data", ""))

    # Info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    logging.info(' ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_values = parser_handle.parse_args()

    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings
    else:
        alg_settings = 'configuration.json'

    return alg_settings


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to read file json
def read_file_json(file_name):
    env_ws = {}
    for env_item, env_value in os.environ.items():
        env_ws[env_item] = env_value

    with open(file_name, "r") as file_handle:
        json_block = []
        for file_row in file_handle:

            for env_key, env_value in env_ws.items():
                env_tag = '$' + env_key
                if env_tag in file_row:
                    env_value = env_value.strip("'\\'")
                    file_row = file_row.replace(env_tag, env_value)
                    file_row = file_row.replace('//', '/')

            # Add the line to our JSON block
            json_block.append(file_row)

            # Check whether we closed our JSON block
            if file_row.startswith('}'):
                # Do something with the JSON dictionary
                json_dict = json.loads(''.join(json_block))
                # Start a new block
                json_block = []

    return json_dict


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Convenience function to convert kilometers to degrees assuming a perfectly spherical Earth.
def kilometers2degrees(kilometer, radius=6371):
    return kilometer / (2.0 * radius * math.pi / 360.0)


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to set logging information
def set_logging(logger_file='log.txt', logger_format=None):
    if logger_format is None:
        logger_format = '%(asctime)s %(name)-12s %(levelname)-8s ' \
                        '%(filename)s:[%(lineno)-6s - %(funcName)20s()] %(message)s'

    # Remove old logging file
    if os.path.exists(logger_file):
        os.remove(logger_file)

    # Set level of root debugger
    logging.root.setLevel(logging.INFO)

    # Open logging basic configuration
    logging.basicConfig(level=logging.INFO, format=logger_format, filename=logger_file, filemode='w')

    # Set logger handle
    logger_handle_1 = logging.FileHandler(logger_file, 'w')
    logger_handle_2 = logging.StreamHandler()
    # Set logger level
    logger_handle_1.setLevel(logging.INFO)
    logger_handle_2.setLevel(logging.INFO)
    # Set logger formatter
    logger_formatter = logging.Formatter(logger_format)
    logger_handle_1.setFormatter(logger_formatter)
    logger_handle_2.setFormatter(logger_formatter)

    # Add handle to logging
    logging.getLogger('').addHandler(logger_handle_1)
    logging.getLogger('').addHandler(logger_handle_2)


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to set logging information
def convertAIIGrid(inFile, outFile, outType, precision=2):
    os.system(
        "gdal_translate -co FORCE_CELLSIZE=YES -q -co DECIMAL_PRECISION=" + str(
            precision) + " -ot " + outType + " -a_nodata -9999.0 -of AAIGrid " + inFile + " " + outFile)
    os.system("rm " + inFile)

    if outType is 'Int16':
        os.system("sed -i s/32768/9999/g " + outFile)


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to fill path names
def fillScriptSettings(data_settings, domain):
    path_settings = {}

    for k, d in data_settings['path'].items():
        for k1, strValue in d.items():
            if isinstance(strValue, str):
                if '{' in strValue:
                    strValue = strValue.replace('{domain}', domain)
            path_settings[k1] = strValue

    return path_settings


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to write HMC pointer conversion
def writeGrass2HMC(ancillary_path):
    file1 = open(os.path.join(ancillary_path, "grass2cont.txt"), "w")
    L = ["1 -1=9\n", "2 -2=8\n", "3 -3=7\n", "4 -4=4\n", "5 -5=1\n", "6 -6=2\n", "7 -7=3\n", "8 -8=6\n", "*=-9999"]
    file1.writelines(L)


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to clean DEM
def cleanlakes(a2dDem, iRows, iCols):
    i = 1
    j = 1
    kk = 1
    kkk = 0

    while kk != 0:
        kk = 0
        while i != iCols:
            while j != iRows:
                exitLoop = 0

                if a2dDem[(j, i)] <= 0:
                    j = j + 1
                    continue
                zmin = 1e20
                zx = a2dDem[(j, i)]

                for iii in np.arange(-1, 2, 1):
                    if exitLoop == 1:
                        break
                    for jjj in np.arange(-1, 2, 1):
                        if iii == 0 and jjj == 0:
                            continue
                        ii = iii + i
                        jj = jjj + j

                        if zx > a2dDem[(jj, ii)]:
                            exitLoop = 1
                            break
                        else:
                            if zmin > a2dDem[(jj, ii)]:
                                zmin = a2dDem[(jj, ii)]

                if exitLoop == 1:
                    j = j + 1
                    continue
                else:
                    if zx <= zmin:
                        a2dDem[(j, i)] = zmin + 0.4
                        kkk = kkk + 1
                        kk = 1
                        if i != 1:
                            i = i - 1
                        if j != 1:
                            j = j - 1

            j = 1
            i = i + 1
            # print('Re-iter ' + str(kkk))

    return a2dDem


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to compute alpha and beta

def makeAlphaBeta(basePath, ancillaryPath, domain):
    logging.info(' ---> Importing data ...')

    Dem_in = rio.open(os.path.join(ancillaryPath, domain + ".dem_unwalled.txt"))
    Dem_in = rio.open(basePath + domain + '.dem.txt')
    iPun = rio.open(basePath + domain + '.pnt.txt')
    iChoice = rio.open(basePath + domain + '.choice.txt')
    AreaCell = rio.open(basePath + domain + '.areacell.txt')
    DD = 50
    a2dAlphaMap = ancillaryPath + 'temp_alpha.tif'
    a2dBetaMap = ancillaryPath + 'temp_beta.tif'

    a2dDem = np.flipud(Dem_in.read(1))
    a2iPun = np.flipud(iPun.read(1))
    a2dDem[a2iPun == -9999] = -9999
    a2iChoice = np.flipud(iChoice.read(1))
    a2dAreaCell = np.flipud(AreaCell.read(1))

    # Variable initialization
    diff_DD = np.ones(a2dDem.shape) * -9999
    LDD = np.zeros(a2dDem.shape)
    mask_perc_tot = np.ones(a2dDem.shape) * -9999
    a2dAlpha = np.ones(a2dDem.shape) * -9999
    a2dBeta = np.ones(a2dDem.shape) * -9999
    pend = np.zeros(a2dDem.shape)
    pend2 = np.zeros(a2dDem.shape)
    pend3 = np.zeros(a2dDem.shape)

    # GEt 2d variable dimensions and cell area mean value
    dDxM = math.sqrt(np.nanmean(a2dAreaCell))
    dDyM = dDxM

    iRows, iCols = Dem_in.shape

    # Checking distance t
    dDistanceT = 500
    if dDxM >= 100 and dDxM < 1000:
        dDistanceT = 2000
    if dDxM >= 5000 and dDxM < 20000:
        dDistanceT = 30000

    # Dem corrections
    logging.info(' ---> Correcting elevation ...')

    a2dDem[(a2dDem <= 0) & (a2dDem > -1000)] = 0.2
    a2dDemOrig = deepcopy(a2dDem)
    # a2dDem=cleanlakes(a2dDemOrig, iRows, iCols)

    # plt.imshow(a2dDemOrig-a2dDemFix)
    # plt.colorbar()
    # plt.show(block=True)

    # cleanlakes(a2dDem)
    # Define alpha matrix angle

    logging.info(' ---> Computing alpha ...')

    for i in np.arange(0, iRows):
        for j in np.arange(0, iCols):

            a = deepcopy(i)
            b = deepcopy(j)

            if a2dDem[(i, j)] > 0:
                fNumPen = 0

                while a2dDem[(a, b)] > 0 and diff_DD[(a, b)] == -9999:
                    if a >= 0 and a < iRows and b >= 0 and b < iCols:  # a cosa serve? ci sono dentro per forza...
                        iii = a + (int((a2iPun[(a, b)] - 1) / 3) - 1)
                        jjj = b + a2iPun[(a, b)] - 5 - 3 * (int((a2iPun[(a, b)] - 1) / 3) - 1)
                        LDD[(a, b)] = math.sqrt(((a - iii) * dDyM) ** 2 + ((b - jjj) * dDxM) ** 2)

                        if iii < 0 or jjj < 0 or iii >= iRows or jjj >= iCols:
                            break

                        diff_DD[(a, b)] = a2dDem[(a, b)] - a2dDem[(iii, jjj)]

                        # Pendenza media sui canali
                        if math.atan2(diff_DD[(a, b)], LDD[(a, b)]) > 0 and diff_DD[(a, b)] < 9000:
                            fNumPen = fNumPen + 1
                            pend[(a, b)] = pend[(a, b)] + math.atan2(diff_DD[(a, b)], LDD[(a, b)])

                        while a2dDem[(a, b)] - a2dDem[
                            (iii, jjj)] <= DD and iii >= 0 and iii < iRows - 1 and jjj >= 0 and jjj < iCols - 1 and \
                                a2dDem[(iii, jjj)] > 0 and LDD[(a, b)] < dDistanceT:
                            diff_DD[(a, b)] = a2dDem[(a, b)] - a2dDem[(iii,
                                                                       jjj)]  # Quesra non andrebbe calcolata dopo essere scesi di cella? (vedi commento sotto) MOD2023
                            ii = iii + (int((a2iPun[(iii, jjj)] - 1) / 3) - 1)
                            jj = jjj + a2iPun[(iii, jjj)] - 5 - 3 * (int((a2iPun[(iii, jjj)] - 1) / 3) - 1)

                            if a2dDem[(a, b)] - a2dDem[
                                (ii, jj)] <= DD and ii >= 0 and ii < iRows and jj >= 0 and jj < iCols:
                                LDD[(a, b)] = LDD[(a, b)] + math.sqrt(
                                    (((ii - iii) * dDyM) ** 2 + ((jj - jjj) * dDxM) ** 2))

                                # Pendenza media sui canali
                                if math.atan2(diff_DD[(a, b)], LDD[(a, b)]) > 0:
                                    if a2iChoice[(a, b)] == 1 and diff_DD[(a, b)] < 9000:
                                        fNumPen = fNumPen + 1
                                        pend[(a, b)] = pend[(a, b)] + math.atan2(diff_DD[(a, b)], LDD[(a,
                                                                                                       b)])  # ma è giusta? diff_DD è riferita a (a,b)(iii,jjj) mentre LDD a (a,b)(ii,jj)

                                    if a2iChoice[(a, b)] == 0 and LDD[(a, b)] < 500 and diff_DD[
                                        (a, b)] < 9000:  # valore fissato o è residuo di dDistanceT?
                                        fNumPen = fNumPen + 1
                                        pend[(a, b)] = pend[(a, b)] + math.atan2(diff_DD[(a, b)], LDD[(a, b)])

                            iii = deepcopy(ii)
                            jjj = deepcopy(jj)

                            if diff_DD[(iii, jjj)] != -9999:  # finale dei percorsi
                                while a2dDem[(a, b)] - a2dDem[
                                    (iii, jjj)] <= DD and iii >= 0 and iii < iRows and jjj >= 0 and jjj < iCols and \
                                        a2dDem[(iii, jjj)] > 0 and LDD[(a, b)] < dDistanceT:
                                    diff_DD[(a, b)] = a2dDem[(a, b)] - a2dDem[(iii, jjj)]
                                    ii = iii + (int((a2iPun[(iii, jjj)] - 1) / 3) - 1)
                                    jj = jjj + a2iPun[(iii, jjj)] - 5 - 3 * (int((a2iPun[(iii, jjj)] - 1) / 3) - 1)

                                    if (a2dDem[(a, b)] - a2dDem[
                                        (ii, jj)]) <= DD and ii >= 0 and ii < iRows and jj >= 0 and jj < iCols:
                                        LDD[(a, b)] = LDD[(a, b)] + math.sqrt(
                                            ((ii - iii) * dDyM) ** 2 + ((jj - jjj) * dDxM) ** 2)

                                        # Pendenza media sui canali
                                        if math.atan2(diff_DD[(a, b)], LDD[(a, b)]) > 0 and diff_DD[(a, b)] < 9000:
                                            if a2iChoice[(a, b)] == 1:
                                                fNumPen = fNumPen + 1
                                                pend[(a, b)] = pend[(a, b)] + math.atan2(diff_DD[(a, b)], LDD[(a, b)])

                                            if a2iChoice[(a, b)] == 0 and LDD[(a, b)] < 500 and diff_DD[(a, b)] < 9000:
                                                fNumPen = fNumPen + 1
                                                pend[(a, b)] = pend[(a, b)] + math.atan2(diff_DD[(a, b)], LDD[(a, b)])

                                    iii = deepcopy(ii)
                                    jjj = deepcopy(jj)

                        if fNumPen > 0:
                            pend[(a, b)] = pend[(a, b)] / fNumPen

                        a2dAlpha[(a, b)] = math.atan2(DD, LDD[(a, b)])  # Angolo in radianti

                        if diff_DD[(a, b)] == 0.9 or diff_DD[(a, b)] > 500:
                            diff_DD[(a, b)] = 0.9

                        if diff_DD[(a, b)] < 1 and LDD[(a, b)] < 4 * dDxM:
                            LDD[(a, b)] = 4 * dDxM

                        a2dAlpha[(a, b)] = math.atan2(diff_DD[(a, b)], LDD[(a, b)])

                        ii = a + (int((a2iPun[(a, b)] - 1) / 3) - 1)
                        jj = b + a2iPun[(a, b)] - 5 - 3 * (int((a2iPun[(a, b)] - 1) / 3) - 1)

                        # if a2dAlpha[(a, b)]<0.1:
                        #    print('ciao')

                        if a2dDem[(ii, jj)] >= 0:
                            a = deepcopy(ii)
                            b = deepcopy(jj)
                            fNumPen = 0
                        else:
                            continue  # esce ma conserva gli indici della fine percorso svolto

                # Fine di un percorso completo seguendo i puntatori
                ii = a + (int((a2iPun[(a, b)] - 1) / 3) - 1)
                jj = b + a2iPun[(a, b)] - 5 - 3 * (int((a2iPun[(a, b)] - 1) / 3) - 1)

    # a2dAlpha[a2dAlpha==-9999]=np.nan
    # plt.imshow(a2dAlpha)
    # plt.colorbar()
    # plt.show(block=True)

    logging.info(' ---> Computing beta ...')
    # Define beta matrix angle
    a2dBeta = deepcopy(pend)
    pend[a2iChoice < 1] = 0

    pend2 = deepcopy(pend)
    pend = pend * 0

    # Smoothing della pendenza sui canali
    for i in np.arange(0, iRows - 1):
        for j in np.arange(0, iCols - 1):
            if a2iChoice[(i, j)] == 1:
                fn = 0
                for ii in np.arange(i - 1, i + 2):
                    for jj in np.arange(j - 1, j + 2):
                        if pend2[(ii, jj)] > 0:
                            fn = fn + 1
                            pend[(i, j)] = pend[(i, j)] + pend2[(ii, jj)]

                if fn < 1:
                    fn = 1
                pend[(i, j)] = pend[(i, j)] / fn

                if LDD[(i, j)] <= 4 * dDxM and diff_DD[(i, j)] < 2:
                    pend[(i, j)] = a2dAlpha[(i, j)]

                if pend[(i, j)] > 0:
                    a2dAlpha[(i, j)] = pend[(i, j)]
                    a2dBeta[(i, j)] = pend[(i, j)]

    dBmin = max(0.0001, np.min(pend[pend > 0]))

    a2dBeta = np.where((a2dDem > 0) & (a2dBeta == 0), a2dAlpha, a2dBeta)
    a2dBeta = np.where((a2dDem > 0) & (a2dBeta < dBmin), dBmin, a2dBeta)
    a2dAlpha = np.where((a2dDem > 0) & (a2dAlpha < 0.0001), 0.0001, a2dAlpha)
    a2dBeta[a2dAlpha == -9999] = -9999

    # plt.imshow(a2dAlpha)
    # plt.colorbar()
    # plt.show(block=True)

    logging.info(' ----> Alpha max: ' + str(np.nanmax(a2dAlpha)))
    logging.info(' ----> Alpha min: ' + str(np.nanmin(a2dAlpha[a2dAlpha != -9999])))

    logging.info(' ----> Beta max: ' + str(np.nanmax(a2dBeta)))
    logging.info(' ----> Beta min: ' + str(np.nanmin(a2dBeta[a2dBeta != -9999])))

    profile = Dem_in.profile
    profile.update(driver='GTiff', count=1)

    logging.info(' ----> Writing alpha and beta ...')

    with rio.open(a2dBetaMap, 'w', **profile) as dst:
        dst.write(np.flipud(a2dBeta.astype('Float32')), 1)

    with rio.open(a2dAlphaMap, 'w', **profile) as dst:
        dst.write(np.flipud(a2dAlpha.astype('Float32')), 1)


# ----------------------------------------------------------------------------
# Call script from external library

if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------
