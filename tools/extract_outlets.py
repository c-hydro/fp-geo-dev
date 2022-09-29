#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fp-geo tools - Extract outlets

__date__ = '20220607'
__version__ = '1.0.0'
__author__ =
        'Andrea Libertino (andrea.libertino@cimafoundation.org')
__library__ = 'fp-geo'

General command line:
python3 extract_outlets.py -settings_file settings.json (-domain domain_name)

Version(s):
20220207 (1.0.0) --> 	Beta release
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
import numpy as np
from xrspatial import zonal_stats
from pysheds.grid import Grid
import time, json, logging, os
from argparse import ArgumentParser
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'fp-geo tools - Extract outlets'
alg_version = '1.0.0'
alg_release = '2022-06-07'
# Algorithm parameter(s)
time_format = '%Y%m%d%H%M'
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Script main
def main():
    start_time = time.time()
    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings, domain = get_args()

    # Set algorithm settings
    data_settings = read_file_json(alg_settings)

    if domain is not None:
        data_settings["algorithm"]["domain"] = domain

    # Set algorithm logging
    log_file = os.path.join(data_settings["data"]["log"]["folder"], data_settings["data"]["log"]["filename"]).replace("{domain}",data_settings["algorithm"]["domain"])
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    set_logging(logger_file=log_file)

    # Initialize folders
    logging.info(" --> Initialize folders...")
    ancillary_folder = data_settings["data"]["ancillary"]["folder"].replace("{domain}",data_settings["algorithm"]["domain"])
    out_folder = data_settings["data"]["output"]["folder"].replace("{domain}",data_settings["algorithm"]["domain"])
    os.makedirs(out_folder, exist_ok=True)
    os.makedirs(ancillary_folder, exist_ok=True)

    ancillary_file = os.path.join(ancillary_folder, "stream_ids.tif")
    output_file = os.path.join(out_folder, data_settings["data"]["output"]["filename"].replace("{domain}",data_settings["algorithm"]["domain"]))
    logging.info(" --> Initialize folders...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Initialize algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> START ... ')
    logging.info(' ==> ALGORITHM SETTINGS <== ')
    logging.info(' ----> Domain: ' + data_settings["algorithm"]["domain"])

    logging.info(" --> Read input data...")
    # Read pysheds input data
    grid = Grid.from_ascii(data_settings["data"]["input"]["pnt"].replace("{domain}",data_settings["algorithm"]["domain"]))
    pnt = grid.read_ascii(data_settings["data"]["input"]["pnt"].replace("{domain}",data_settings["algorithm"]["domain"]))
    pnt.nodata = -9999
    choice = grid.read_ascii(data_settings["data"]["input"]["choice"].replace("{domain}",data_settings["algorithm"]["domain"]))
    logging.info(" --> Read input data... DONE")

    logging.info(" --> Classify choice...")
    # Classify choice raster with unique ids for streams
    dirmap_HMC = (8, 9, 6, 3, 2, 1, 4, 7)
    streams = grid.extract_river_network(fdir=pnt, mask=choice==1, dirmap=dirmap_HMC, routing='d8')
    streams_coll = zip([i["geometry"] for i in streams["features"]], [i["id"] for i in streams["features"]])
    river_raster = grid.rasterize(streams_coll)
    grid.to_raster(river_raster,ancillary_file)
    logging.info(" --> Classify choice...DONE")

    logging.info(" --> Read data for zonal stats...")
    # Read rioxarray files for zonal stats
    acc = rxr.open_rasterio(data_settings["data"]["input"]["area"].replace("{domain}",data_settings["algorithm"]["domain"]))
    choice = rxr.open_rasterio(data_settings["data"]["input"]["choice"].replace("{domain}", data_settings["algorithm"]["domain"]))
    acc.values = np.where(choice.values<0,-9999,acc.values)
    streams = rxr.open_rasterio(ancillary_file)
    logging.info(" --> Read data for zonal stats...DONE")

    logging.info(" --> Find outlets...")
    logging.info(" ---> Calculate maximum accumulation per stream...")
    # Calculate zonal stats for finding maximum accumulation per stream
    streams.values = np.nan_to_num(streams.values, copy=False).astype(int)
    max_zones = zonal_stats(zones=streams.squeeze(), values=acc.squeeze(), stats_funcs=['max'])
    out_array = streams.copy()
    out_array = out_array * np.nan

    # Generate a raster with stream_id in correspondance of the maximum accumulatio nof the stream
    for stream_id in np.unique(streams.values):
        if stream_id == 0:
            continue
        out_array = xr.where((streams == stream_id) & (acc == max_zones.loc[max_zones.zone==stream_id, "max"].values), stream_id, out_array)

    logging.info(" ---> Find coordinates...")
    # Extract coordinate of non-null cells
    da_stacked = out_array.squeeze().stack(p=['y','x'])
    da_stacked = da_stacked[da_stacked.notnull()]

    logging.info(" ---> Write output...")

    gdf = gpd.GeoDataFrame({"id" : [data_settings["algorithm"]["domain"] + "P" + str(int(i)).zfill(5) for i in da_stacked.values], "domain" : data_settings["algorithm"]["domain"]}, index= da_stacked.values.astype(int), geometry=gpd.points_from_xy(da_stacked.x.values, da_stacked.y.values))

    for val, x, y in zip(gdf.index.values, da_stacked.x.values, da_stacked.y.values):
        gdf.loc[val,["x_HMC", "y_HMC"]] = [int(grid.nearest_cell(x,y,snap='center')[1] + 1) , int(grid.nearest_cell(x,y,snap='center')[0] + 1)]

    gdf = gdf.astype({"id": str, "domain": str, "x_HMC" : int, "y_HMC" : int})
    gdf.set_crs(epsg=4326, inplace=True).to_file(output_file)
    gdf.to_csv(output_file.replace(".shp",".txt"), columns=["x_HMC", "y_HMC", "domain", "id"], sep=' ', header=False, index=False)

    logging.info(" --> Find outlets...DONE")

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Clean system
    logging.info(" --> Clean system...")
    if data_settings["algorithm"]["flags"]["clear_ancillary"]:
        os.system("rm -r " + ancillary_file)
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

# -------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_handle.add_argument('-domain', action="store", dest="domain")
    parser_values = parser_handle.parse_args()

    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings
    else:
        alg_settings = 'configuration.json'

    if parser_values.domain:
        domain = parser_values.domain
    else:
        domain = None


    return alg_settings, domain

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





