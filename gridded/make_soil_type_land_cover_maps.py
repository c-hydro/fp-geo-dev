#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fp-geo - Make Soil Type Land Cover Maps

__date__ = '20220430'
__version__ = '2.0.2'
__author__ =
        'Andrea Libertino (andrea.libertino@cimafoundation.org')
__library__ = 'fp-geo'

General command line:
python3 MakeSoilMaps

Version(s):
20220430 (2.0.2) -->    Added support to longitudes crossing the 180 degree for Fiji
                        Added first attempt Field Capacity estimation with Saxton et al., 1986
20211203 (2.0.1) -->    Add vegetation layers procedures
20210820 (2.0.0) -->	Unified soil layer procedures and add download capabilities
20201016 (1.1.0) -->	Add support to setting file and to external LULC maps
20200430 (1.0.0) --> 	Beta release on the Apollo server
"""

# -------------------------------------------------------------------------------------
# Complete library

import pandas as pd
from owslib.wcs import WebCoverageService
from pyproj.transformer import Transformer
import os
from requests.exceptions import HTTPError
import numpy as np
import rasterio as rio
from copy import deepcopy
import time
from osgeo import gdal, gdalconst
from scipy import ndimage as nd
import logging
import json
from argparse import ArgumentParser
from datetime import date
import warnings
from copy import deepcopy

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'fp-geo tools - Make Soil Type Land Cover Maps'
alg_version = '2.0.2'
alg_release = '2022-04-30'
# Algorithm parameter(s)
time_format = '%Y%m%d%H%M'


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Script main
def main():
    start_time = time.time()
    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings = get_args()

    # Set algorithm settings
    data_settings = read_file_json(alg_settings)

    # Set algorithm logging
    log_file = os.path.join(data_settings["log"]["folder"], data_settings["log"]["file_name"]).format(
        domain=data_settings["algorithm"]["general"]["domain"])
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    set_logging(logger_file=log_file)

    # Initialize folders
    logging.info(" --> Initialize folders...")
    ancillary_folder = data_settings["path"]["ancillary"].format(domain=data_settings["algorithm"]["general"]["domain"])
    out_folder = data_settings["path"]["output"].format(domain=data_settings["algorithm"]["general"]["domain"])
    os.makedirs(out_folder, exist_ok=True)
    os.makedirs(ancillary_folder, exist_ok=True)
    logging.info(" --> Initialize folders...DONE")

    gdal.SetConfigOption("GDAL_HTTP_UNSAFESSL", "YES")
    gdal.VSICurlClearCache()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Initialize algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME : ' + date.today().strftime("%d-%B-%Y %H:%m"))
    logging.info(' ==> START ... ')
    logging.info(' ==> ALGORITHM SETTINGS <== ')
    logging.info(' ----> Domain: ' + data_settings["algorithm"]["general"]["domain"])

    logging.info(' ')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Check on settings validity
    logging.info(" --> Check flags validity...")
    for input in ["soil_maps", "land_cover"]:
        flagged_settings = len(
            [data_settings["algorithm"]["flags"][input][x] for x in data_settings["algorithm"]["flags"][input] if
             data_settings["algorithm"]["flags"][input][x]])
        if flagged_settings > 1:
            logging.error(' ----> ERROR! Please choose if use local ' + input + ' maps or download it')
            raise ValueError(input + " sources flags are mutually exclusive!")
        elif flagged_settings < 1:
            logging.error(' ----> ERROR! Please choose if use local ' + input + ' maps or download it')
            raise ValueError("At least one " + input + " sources flags has to be choosen!")
    logging.info(" --> Check flags validity...DONE")

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get grid model info
    logging.info(" --> Derive grid from input dem...")
    proj_settings = data_settings["advanced"]["epsg_settings"]
    proj_dem = proj_settings[data_settings["data"]["dem"]["srs"]]

    model_srs = {}
    match_ds = gdal.Open(
        data_settings["data"]["dem"]["file_path"].format(domain=data_settings["algorithm"]["general"]["domain"]),
        gdalconst.GA_ReadOnly)
    out_profile = rio.open(data_settings["data"]["dem"]["file_path"].format(
        domain=data_settings["algorithm"]["general"]["domain"])).profile
    out_profile["driver"] = "GTiff"
    dem = np.array(match_ds.GetRasterBand(1).ReadAsArray())
    model_srs['proj'] = proj_dem
    model_srs['geotrans'] = match_ds.GetGeoTransform()
    model_srs['wide'] = match_ds.RasterXSize
    model_srs['high'] = match_ds.RasterYSize

    lon_min = min(model_srs['geotrans'][0], model_srs['geotrans'][0] + model_srs['wide'] * model_srs['geotrans'][1])
    lon_max = max(model_srs['geotrans'][0], model_srs['geotrans'][0] + model_srs['wide'] * model_srs['geotrans'][1])
    lat_min = min(model_srs['geotrans'][3], model_srs['geotrans'][3] + model_srs['high'] * model_srs['geotrans'][5])
    lat_max = max(model_srs['geotrans'][3], model_srs['geotrans'][3] + model_srs['high'] * model_srs['geotrans'][5])
    model_srs['bbox'] = [lon_min, lon_max, lat_min, lat_max]
    logging.info(" --> Derive grid from input dem...DONE")

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Compute soil features (for CN and soil layers)
    if data_settings["algorithm"]["flags"]["make_soil_maps"] or data_settings["algorithm"]["flags"]["make_cn_map"]:
        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Get soil data
        Depths = np.array([0, 5, 15, 30, 60, 100])

        # Download from SoilGrids wcs
        if data_settings["algorithm"]["flags"]["soil_maps"]["download_soil_grids"]:
            logging.info(" --> Download soil maps...")
            logging.info(" ---> Identify grid...")
            proj_soil = proj_settings["Homolosine"]
            transformer = Transformer.from_crs("epsg:4326", "+proj=igh +datum=WGS84 +no_defs +towgs84=0,0,0")
            xll, yll = transformer.transform(lat_min, lon_min)
            xur, yur = transformer.transform(lat_max, lon_max)
            subsets = [('X', xll, xur), ('Y', yll, yur)]
            logging.info(" ---> Identify grid...DONE!")

            depths = [str(Depths[i_d]) + "-" + str(Depths[i_d + 1]) + "cm" for i_d in np.arange(0, len(Depths) - 1)]

            for var in ["clay", "sand"]:
                logging.info(" ---> Compute " + var + " maps...")
                wcs_address = data_settings["data"]["soil_texture"]["download_soil_grids"]["wcs_server"].format(var=var)
                logging.info(" ----> Connect to: " + wcs_address)
                wcs = WebCoverageService(wcs_address, version='2.0.1')

                for depth, h in zip(depths, np.diff(Depths)):
                    logging.info(" ----> Depth: " + depth)
                    logging.info(" ----> Download map...")
                    var_depth = data_settings["data"]["soil_texture"]["download_soil_grids"]["map_id"].format(
                        depth=depth, var=var)
                    map_ref = wcs.contents[var_depth]

                    # Get Remote Coverage managing server connection problems that sometimes arise
                    connected = False
                    n_tries = 0
                    while not connected:
                        try:
                            response = wcs.getCoverage(
                                identifier=[var_depth],
                                crs=proj_soil,
                                subsets=subsets,
                                resx=250, resy=250,
                                format=map_ref.supportedFormats[0])
                            connected = True
                        except HTTPError:
                            n_tries = n_tries + 1
                            logging.info("----> Re-connection attempt #" + str(n_tries).zfill(2) + " in 3 seconds...")
                            time.sleep(3)
                            if n_tries > 10:
                                logging.info("ERROR! Can not get remote coverage in 10 attempts!")
                                raise HTTPError
                            pass

                    logging.info(" ----> Writing output map...")
                    with open(os.path.join(ancillary_folder, var_depth + ".tif"), 'wb') as file:
                        file.write(response.read())

                    if depth == depths[0]:
                        file = rio.open(os.path.join(ancillary_folder, var_depth + ".tif"))
                        map = file.read(1).astype(np.int32) * h
                        profile = file.profile
                    else:
                        map = map + deepcopy(
                            rio.open(os.path.join(ancillary_folder, var_depth + ".tif")).read(1).astype(np.int32)) * h
                    logging.info(" ----> Depth: " + depth + "...DONE")

                map = map / (np.sum(np.diff(Depths)))
                mask = np.where(map == 0, True, False)
                map = fill(map, invalid=mask)

                with rio.open(os.path.join(ancillary_folder, "average_" + var + "_to_reproj.tif"), 'w',
                              **profile) as dst:
                    dst.write(map.astype(np.int16), 1)
                write_tif_from_grid(proj_soil, model_srs,
                                    os.path.join(ancillary_folder, "average_" + var + "_to_reproj.tif"),
                                    os.path.join(ancillary_folder, "average_regrid_" + var + ".tif"),
                                    gdalconst.GRA_Bilinear)

                logging.info(" ---> Compute " + var + " maps...DONE")
            logging.info(" --> Download soil maps...DONE")

        # Import existing sources
        elif data_settings["algorithm"]["flags"]["soil_maps"]["use_local_data"]:
            logging.info(" --> Import local soil maps...")
            proj_soil = proj_settings[data_settings["data"]["soil_texture"]["local_source_settings"]["srs"]]
            for var in ["sand", "clay"]:
                if os.path.isfile(data_settings["data"]["soil_texture"]["local_source_settings"]["file_path_" + var]):
                    write_tif_from_grid(proj_soil, model_srs,
                                        data_settings["data"]["soil_texture"]["local_source_settings"][
                                            "file_path_" + var],
                                        os.path.join(ancillary_folder, "temp_average_regrid_" + var + ".tif"),
                                        gdalconst.GRA_Bilinear)

                    map = rio.open(os.path.join(ancillary_folder, "temp_average_regrid_" + var + ".tif"))
                    profile = map.profile
                    mask = np.where(map.read(1) == 0, True, False)
                    map_filled = fill(map.read(1), invalid=mask)
                    with rio.open(os.path.join(ancillary_folder, "average_regrid_" + var + ".tif"), 'w',
                                  **profile) as dst:
                        dst.write(map_filled, 1)
                else:
                    logging.error("ERROR! " + var + " map is mandatory when 'use_local_soil_texture' is flagged!")
                    raise FileNotFoundError(
                        "Cannot find file " + data_settings["data"]["soil_texture"]["local_source_settings"][
                            "file_path_" + var])
            logging.info(" --> Import local soil maps...DONE")

        logging.info(" --> Longitude max larger tha 180, correcting...")
        if lon_max > 180:
            for var in ["clay", "sand"]:
                soil_map = rio.open(os.path.join(ancillary_folder, "average_regrid_" + var + ".tif"))
                last_col = soil_map.index(180, 0)[1]
                n_cols = soil_map.shape[1] - last_col - 1
                model_srs_edit = deepcopy(model_srs)
                temp = list(model_srs_edit['geotrans'])
                temp[0] = -180
                model_srs_edit['geotrans'] = tuple(temp)
                model_srs_edit['wide'] = n_cols
                model_srs_edit['bbox'][0] = -180
                model_srs_edit['bbox'][1] = -180 + (model_srs_edit['bbox'][1] - 180)
                if data_settings["algorithm"]["flags"]["soil_maps"]["download_soil_grids"]:
                    logging.error(
                        " --> If longitude exceeds 180° only local data can be used for soil characterization")
                    raise NotImplementedError
                else:
                    write_tif_from_grid(proj_soil, model_srs_edit,
                                        data_settings["data"]["soil_texture"]["local_source_settings"][
                                            "file_path_" + var],
                                        os.path.join(ancillary_folder, "temp2_average_regrid_" + var + ".tif"),
                                        gdalconst.GRA_Bilinear)
                    map = rio.open(os.path.join(ancillary_folder, "temp2_average_regrid_" + var + ".tif"))
                    profile = map.profile
                    mask = np.where(map.read(1) == 0, True, False)
                    map_filled = fill(map.read(1), invalid=mask)
                    with rio.open(os.path.join(ancillary_folder, "temp3_average_regrid_" + var + ".tif"), 'w',
                                  **profile) as dst:
                        dst.write(map_filled, 1)

                soil_map_temp = rio.open(os.path.join(ancillary_folder, "temp3_average_regrid_" + var + ".tif"))
                soil_west = soil_map.read(1).astype(np.float32)
                soil_east = soil_map_temp.read(1).astype(np.float32)

                soil_west[:, last_col + 1:] = soil_east
                out_meta = deepcopy(soil_map.meta)
                with rio.open(os.path.join(ancillary_folder, "average_regrid_" + var + ".tif"), 'w', **out_meta) as dst:
                    dst.write_band(1, soil_west.astype(rio.float32))
                os.remove(os.path.join(ancillary_folder, "temp3_average_regrid_" + var + ".tif"))
                os.remove(os.path.join(ancillary_folder, "temp2_average_regrid_" + var + ".tif"))
        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Derive soil maps
        logging.info(" --> Compute soil features...")

        # Import sand and clay average for calculation and export them in Continuum format
        logging.info(" ---> Export average clay and sand fraction...")
        sand_map = rio.open(os.path.join(ancillary_folder, "average_regrid_sand.tif"))
        sand = sand_map.read(1) / 10
        clay_map = rio.open(os.path.join(ancillary_folder, "average_regrid_clay.tif"))
        clay = clay_map.read(1) / 10
        logging.info(" ---> Compute average clay and sand fraction...DONE")

        # Compute first attempt field capacity with Saxtoon et al., 1986
        a = np.exp(-4.396-0.0715*clay-0.000488*(sand**2)-0.00004285*(sand**2)*clay)
        b = -3.14-0.00222*(clay**2)-0.00003484*(sand**2)*clay
        fc = (0.33333/a)**(1/b)

        # Compute hydrological soil groups by assigning a soil class 1-12 according to USDA classification, a map is generated in the "ancillary" folder for CN
        if data_settings["algorithm"]["flags"]["make_cn_map"]:
            logging.info(" ---> Compute hydrological soil type...")
            silt = 1 - sand/100 - clay/100

            soil_class = -9999 * np.ones_like(sand)
            soil_class = classify_usda(soil_class, sand/100, clay/100, silt)
            soil_class[dem < -9000] = -9999

            # profile = sand_map.profile
            with rio.open(os.path.join(ancillary_folder, "soil_texture.tif"), 'w', **out_profile) as dst:
                dst.write(soil_class.astype(np.float32), 1)

            # Assign to each soil class an hydrological soil type, a map is generated in the "ancillary" folder
            hsg = -9999 * np.ones_like(dem)
            hsg = classify_hsg(hsg, soil_class)
            with rio.open(os.path.join(ancillary_folder, "hydrological_soil_group.tif"), 'w', **out_profile) as dst:
                dst.write(hsg.astype(np.float32), 1)

            logging.info(" ---> Compute hydrological soil type...DONE")

        # Write average clay and average sand maps for soil layers calibration
        if data_settings["algorithm"]["flags"]["make_soil_maps"]:
            logging.info(" ---> Write average clay and sand fraction...")
            # profile = sand_map.profile
            with rio.open(os.path.join(ancillary_folder, "sand_perc.tif"), 'w', **out_profile) as dst:
                dst.write((sand_map.read(1) / 10).astype(np.float32), 1)
            convertAIIGrid(os.path.join(ancillary_folder, "sand_perc.tif"),
                           os.path.join(out_folder, data_settings["algorithm"]["general"]["domain"] + ".AVG_sand.txt"),
                           'Int16')
            with rio.open(os.path.join(ancillary_folder, "clay_perc.tif"), 'w', **out_profile) as dst:
                dst.write((clay_map.read(1) / 10).astype(np.float32), 1)
            convertAIIGrid(os.path.join(ancillary_folder, "clay_perc.tif"),
                           os.path.join(out_folder, data_settings["algorithm"]["general"]["domain"] + ".AVG_clay.txt"),
                           'Int16')
            with rio.open(os.path.join(ancillary_folder, "ct_first_attempt.tif"), 'w', **out_profile) as dst:
                dst.write(fc.astype(np.float32), 1)
            convertAIIGrid(os.path.join(ancillary_folder, "ct_first_attempt.tif"),
                           os.path.join(out_folder, data_settings["algorithm"]["general"]["domain"] + ".ct.txt"),
                           'Float32')
        logging.info(" ---> Write average clay and sand fraction...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Compute land cover features (for CN and vegetation layers)
    if data_settings["algorithm"]["flags"]["make_vegetation_maps"] or data_settings["algorithm"]["flags"][
        "make_cn_map"]:
        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Get land cover data

        # Download from esa-cci website
        if data_settings["algorithm"]["flags"]["land_cover"]["download_esa_cci_lc"]:
            proj_lulc = proj_settings["EPSG:4326"]
            logging.info(" --> Download land cover map...")
            write_tif_from_grid(proj_lulc, model_srs,
                                '/vsicurl/' + data_settings["data"]["land_cover"]["download_esa_cci"]["download_path"],
                                os.path.join(ancillary_folder, 'LULC_regrid.tif'), gdalconst.GRA_Mode,
                                file_type=gdalconst.GDT_Float32)
            logging.info(" --> Download land cover map...DONE")

        # Import existing source
        elif data_settings["algorithm"]["flags"]["land_cover"]["use_local_data"]:
            logging.info(" --> Use local land cover map...")
            proj_lulc = proj_settings[data_settings["data"]["land_cover"]["local_source_settings"]["srs"]]
            if os.path.isfile(data_settings["data"]["land_cover"]["local_source_settings"]["file_path"]):
                write_tif_from_grid(proj_lulc, model_srs,
                                    data_settings["data"]["land_cover"]["local_source_settings"]["file_path"],
                                    os.path.join(ancillary_folder, "LULC_regrid.tif"), gdalconst.GRA_Mode,
                                    file_type=gdalconst.GDT_Float32)
            else:
                logging.error("ERROR! Land cover map is mandatory when 'use_local_land_cover' is flagged!")
                raise FileNotFoundError(
                    "Cannot find file " + data_settings["data"]["land_cover"]["local_source_settings"]["file_path"])
            logging.info(" --> Import local soil maps...DONE")

        logging.info(" --> Longitude max larger tha 180, correcting...")
        if lon_max > 180:
            lulc_map = rio.open(os.path.join(ancillary_folder, 'LULC_regrid.tif'))
            last_col = lulc_map.index(180, 0)[1]
            n_cols = lulc_map.shape[1] - last_col - 1
            model_srs_edit = deepcopy(model_srs)
            temp = list(model_srs_edit['geotrans'])
            temp[0] = -180
            model_srs_edit['geotrans'] = tuple(temp)
            model_srs_edit['wide'] = n_cols
            model_srs_edit['bbox'][0] = -180
            model_srs_edit['bbox'][1] = -180 + (model_srs_edit['bbox'][1] - 180)
            if data_settings["algorithm"]["flags"]["land_cover"]["download_esa_cci_lc"]:
                write_tif_from_grid(proj_lulc, model_srs_edit,
                                    '/vsicurl/' + data_settings["data"]["land_cover"]["download_esa_cci"][
                                        "download_path"], os.path.join(ancillary_folder, 'temp_LULC_regrid.tif'),
                                    gdalconst.GRA_Mode, file_type=gdalconst.GDT_Float32)
            else:
                write_tif_from_grid(proj_lulc, model_srs_edit,
                                    data_settings["data"]["land_cover"]["local_source_settings"]["file_path"],
                                    os.path.join(ancillary_folder, "temp_LULC_regrid.tif"), gdalconst.GRA_Mode,
                                    file_type=gdalconst.GDT_Float32)
            lulc_map_temp = rio.open(os.path.join(ancillary_folder, 'temp_LULC_regrid.tif'))
            lulc_west = lulc_map.read(1).astype(np.float32)
            lulc_east = lulc_map_temp.read(1).astype(np.float32)

            lulc_west[:, last_col + 1:] = lulc_east
            out_meta = deepcopy(lulc_map.meta)
            with rio.open(os.path.join(ancillary_folder, 'LULC_regrid.tif'), 'w', **out_meta) as dst:
                dst.write_band(1, lulc_west.astype(rio.float32))
            os.remove(os.path.join(ancillary_folder, 'temp_LULC_regrid.tif'))

        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Derive land cover maps
        lulc_map = rio.open(os.path.join(ancillary_folder, 'LULC_regrid.tif'))
        lulc = lulc_map.read(1).astype(np.float32)

        # Calculate cn with a lookup table assigning the land cover class to the cn according to the hsg.
        # The standard table is derived from Jaafar et al., 2019
        if data_settings["algorithm"]["flags"]["make_cn_map"]:

            logging.info(" ---> Calculate CN...")
            lulc_to_cn = pd.read_csv(data_settings["data"]["land_cover"]["conversion_cn_table"], header=0,
                                     usecols=["ID", "A", "B", "C", "D"], index_col="ID", sep=',')

            cn_map = -9999 * np.ones_like(dem)
            classified_maps = lulc_to_cn.loc[lulc.flatten()]

            for num, hsg_class in enumerate(["A", "B", "C", "D"], start=1):
                cn_map = np.where(hsg == num, classified_maps[hsg_class].values.reshape(cn_map.shape), cn_map)

            with rio.open(os.path.join(ancillary_folder, "temp_cn.tif"), 'w', **out_profile) as dst:
                dst.write(cn_map.astype(np.float32), 1)

            convertAIIGrid(os.path.join(ancillary_folder, "temp_cn.tif"),
                           os.path.join(out_folder, data_settings["algorithm"]["general"]["domain"] + ".cn.txt"),
                           'Int16')
            logging.info(" ---> Calculate CN...DONE")

        # Calculate vegetation layers with a lookup table assigning the land cover class to the different veg variables
        if data_settings["algorithm"]["flags"]["make_vegetation_maps"]:
            logging.info(" ---> Calculate vegetation maps...")
            lulc_to_veg = pd.read_csv(data_settings["data"]["land_cover"]["conversion_veg_table"], header=0,
                                      usecols=["ID", "BareSoil", "RSmin", "Hveg", "Gd"], index_col="ID", sep=',')

            var_type = {"BareSoil": "Int16", "RSmin": "Float32", "Gd": "Float32", "Hveg": "Float32"}

            for var_veg in ["BareSoil", "RSmin", "Hveg", "Gd"]:
                var_map = lulc_to_veg.loc[lulc.flatten()][var_veg].values.reshape(dem.shape)
                var_map = np.where(dem < -9000, -9999, var_map)
                with rio.open(os.path.join(ancillary_folder, "temp_" + var_veg + ".tif"), 'w', **out_profile) as dst:
                    dst.write(var_map.astype(np.float32), 1)
                convertAIIGrid(os.path.join(ancillary_folder, "temp_" + var_veg + ".tif"), os.path.join(out_folder,
                                                                                                        data_settings[
                                                                                                            "algorithm"][
                                                                                                            "general"][
                                                                                                            "domain"] + "." + var_veg + ".txt"),
                               var_type[var_veg])
                logging.info(" ---> Calculate " + var_veg + "...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Clean system
    logging.info(" --> Cleaning system...")
    os.system('rm ' + os.path.join(out_folder, "*.xml") + " || True")

    if data_settings["algorithm"]["flags"]["clear_ancillary"]:
        os.system("rm -r " + ancillary_folder)

    logging.info(" --> Cleaning system...DONE")

    # Info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    logging.info(' ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Function to regrid a raster dataset over the grid of a model dataset
def write_tif_from_grid(src_proj, model_srs, src_filename, dst_filename, resample_type,
                        file_type=gdalconst.GDT_Float32):
    """
        src_proj : projection of the map to convert, provided in OGG WKT format
        model_srs : dictinary including the information on the model grid: proj, geotransf, wide, high
        src_filename : filename of the map to convert
        dst_filename : filename of the converted map
        resample_type: type of resample to apply
    """
    src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
    dst = gdal.GetDriverByName('GTiff').Create(dst_filename, model_srs['wide'], model_srs['high'], 1, file_type)
    dst.SetGeoTransform(model_srs['geotrans'])
    dst.SetProjection(model_srs['proj'])
    gdal.ReprojectImage(src, dst, src_proj, model_srs['proj'], resample_type)


# ----------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Function to fill values of a numpy array with nearest neighhbours
def fill(data, invalid=None):
    if invalid is None: invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid,
                                    return_distances=False,
                                    return_indices=True)
    return data[tuple(ind)]


# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Function to generate the soil classes map
def classify_usda(out_map, sand, clay, silt):
    # clay
    out_map = np.where((clay >= 0.4) & (sand <= 0.45) & (silt < 0.4), 1, out_map)
    # silty-clay
    out_map = np.where((clay >= 0.4) & (silt >= 0.4), 2, out_map)
    # silty-clay-loam
    out_map = np.where((clay >= 0.27) & (clay < 0.4) & (sand <= 0.2), 3, out_map)
    # sandy-clay
    out_map = np.where((clay >= 0.35) & (sand >= 0.45), 4, out_map)
    # sandy-clay-loam
    out_map = np.where((clay >= 0.2) & (clay < 0.35) & (silt < 0.28) & (sand > 0.45), 5, out_map)
    # clay-loam
    out_map = np.where((clay >= 0.27) & (clay < 0.4) & (sand > 0.2) & (sand <= 0.45), 6, out_map)
    # silt
    out_map = np.where((silt >= 0.8) & (clay < 0.12), 7, out_map)
    # silt-loam
    out_map = np.where(
        ((silt >= 0.5) & (clay >= 0.12) & (clay < 0.27)) | ((silt >= 0.5) & (silt < 0.8) & (clay < 0.12)), 8, out_map)
    # loam
    out_map = np.where((clay >= 0.07) & (clay <= 0.27) & (silt >= 0.28) & (silt < 0.5) & (sand <= 0.52), 9, out_map)
    # sand
    out_map = np.where((silt + 1.5 * clay) < 0.15, 10, out_map)
    # loamy-sand
    out_map = np.where(((silt + 1.5 * clay) >= 0.15) & ((silt + 2 * clay) < 0.3), 11, out_map)
    # sandy-loam
    out_map = np.where(((clay >= 0.07) & (clay <= 0.2) & (sand > 0.52) & ((silt + 2 * clay) >= 0.3)) | (
                (clay < 0.07) & (silt < 0.5) & ((silt + 2 * clay) >= 0.3)), 12, out_map)
    # fill nan
    out_map = np.where(sand == 0, -9999, out_map)
    return out_map


# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Function to generate the hydrologicla soil type map
def classify_hsg(out_map, soil_class):
    # class A:1
    out_map = np.where(soil_class == 10, 1, out_map)
    # class B:2
    out_map = np.where((soil_class == 11) | (soil_class == 12), 2, out_map)
    # class C:3
    out_map = np.where(((soil_class > 4) & (soil_class < 10)) | (soil_class == 3), 3, out_map)
    # class D:4
    out_map = np.where((soil_class == 1) | (soil_class == 2) | (soil_class == 4), 4, out_map)
    return out_map


# ----------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to set logging information
def convertAIIGrid(inFile, outFile, outType, precision=2):
    os.system(
        "gdal_translate -co FORCE_CELLSIZE=YES -q -co DECIMAL_PRECISION=" + str(
            precision) + " -ot " + outType + " -a_nodata -9999.0 -of AAIGrid " + inFile + " " + outFile)
    os.system("rm " + inFile)

    if outType is 'Int16':
        os.system("sed -i s/32768/9999/g " + outFile)

    if outType is 'Int32':
        os.system("sed -i s/-2147483648/-9999/g " + outFile)


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

# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------