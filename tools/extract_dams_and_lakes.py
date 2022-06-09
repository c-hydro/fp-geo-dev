#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fp-geo tools - Extract dams and lakes

__date__ = '20220607'
__version__ = '1.0.0'
__author__ =
        'Andrea Libertino (andrea.libertino@cimafoundation.org')
__library__ = 'fp-geo'

General command line:
python3 extract_dams_and_lakes.py -settings_file settings.json (-domain domain_name)

Version(s):
20220207 (1.0.0) --> 	Beta release
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
import geopandas as gpd
import rioxarray as rxr
import time, json, logging, os
from argparse import ArgumentParser
import shapely.geometry as sg
import rasterio
import numpy as np

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'fp-geo tools - Extract dams and lakes'
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

    ancillary_file = None
    output_file_dams = os.path.join(out_folder, data_settings["data"]["output"]["filename_dams"].replace("{domain}",data_settings["algorithm"]["domain"]))
    output_file_lakes = os.path.join(out_folder, data_settings["data"]["output"]["filename_lakes"].replace("{domain}",data_settings["algorithm"]["domain"]))

    logging.info(" --> Initialize folders...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Initialize algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> START ... ')
    logging.info(' ==> ALGORITHM SETTINGS <== ')
    logging.info(' ----> Domain: ' + data_settings["algorithm"]["domain"])

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Setting the domain
    # Import choice for domain mask
    choice = rxr.open_rasterio(data_settings["data"]["input"]["choice"].replace("{domain}", data_settings["algorithm"]["domain"]))

    logging.info(" --> Rasterize domain...")
    domain_shape = polygonize(choice.squeeze()).set_crs(epsg=4326, inplace=True)
    domain_shape['dummy'] = 'dummy'  # add dummy column to dissolve all geometries into one
    geom = domain_shape.dissolve(by='dummy').geometry[0]
    logging.info(" --> Rasterize domain...DONE")

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Dams selection
    logging.info(" --> Compute dams...")

    logging.info(" ---> Read and filter dams database...")
    dams_db = gpd.read_file(data_settings["data"]["database"]["dams"])
    dams_in = dams_db[dams_db.within(geom)]
    dams_in = dams_in[dams_in[data_settings["algorithm"]["tresholds_vol_MCM"]["dams_field"]] >= data_settings["algorithm"]["tresholds_vol_MCM"]["dams_value"]]

    logging.info(" ---> Verify if dams lie on the network...")
    coord_list = [(x, y) for x, y in zip(dams_in['geometry'].x, dams_in['geometry'].y)]
    with rasterio.open(data_settings["data"]["input"]["choice"].replace("{domain}", data_settings["algorithm"]["domain"])) as src:
        dams_in["choice"] = [x[0] for x in src.sample(coord_list)]

    logging.info(" ---> Write dam output...")
    if not len(dams_in)==0:
        dams_in.to_file(output_file_dams)
    else:
        logging.warning(" --> WARNING! No dams in the domain with the current settings!")

    logging.info(" --> Compute dams...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Dams selection
    logging.info(" --> Compute lakes...")

    logging.info(" ---> Read and filter lakes database...")
    choice_extent = str(min(choice.x.values)) + " " + str(min(choice.y.values)) + " " + str(max(choice.x.values)) + " " + str(max(choice.y.values))
    temp_lakes = os.path.join(ancillary_folder,"lakes_" + data_settings["algorithm"]["domain"] + "_temp.shp")
    os.system("ogr2ogr " + temp_lakes + " " + data_settings["data"]["database"]["lakes"] + " -clipsrc " + choice_extent)
    lakes_db = gpd.read_file(temp_lakes)
    lakes_in = lakes_db[lakes_db.within(geom)]
    lakes_in = lakes_in[lakes_in[data_settings["algorithm"]["tresholds_vol_MCM"]["lakes_field"]] >= data_settings["algorithm"]["tresholds_vol_MCM"]["lakes_value"]]

    if len(dams_in) > 0:
        logging.info(" ---> Remove lakes near to dams...")
        buff_dist = km2deg(data_settings["algorithm"]["buffer_dam_lakes_km"])
        dams_buffer = dams_in.copy()
        dams_buffer.geometry = dams_buffer.geometry.buffer(buff_dist)
        dams_buffer['dummy'] = 'dummy'  # add dummy column to dissolve all geometries into one
        geom = dams_buffer.dissolve(by='dummy').geometry[0]
        lakes_in = lakes_in[~lakes_in.within(geom)]

    if not len(lakes_in) == 0:
        logging.info(" ---> Verify if dams lie on the network...")
        coord_list = [(x, y) for x, y in zip(lakes_in['geometry'].x, lakes_in['geometry'].y)]
        with rasterio.open(data_settings["data"]["input"]["choice"].replace("{domain}", data_settings["algorithm"]["domain"])) as src:
            lakes_in["choice"] = [x[0] for x in src.sample(coord_list)]

        logging.info(" ---> Write lakes output...")
        lakes_in.to_file(output_file_lakes)
    else:
        logging.warning(" --> WARNING! No dams in the domain with the current settings!")

    os.system("rm " + temp_lakes.replace(".shp",".*"))

    logging.info(" --> Compute lakes...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Clean system
    logging.info(" --> Clean system...")
    if data_settings["algorithm"]["flags"]["clear_ancillary"] and ancillary_file is not None:
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

def polygonize(da):
    if da.dims != ("y", "x"):
        raise ValueError('Dimensions must be ("y", "x")')

    values = da.values
    values[values>0] = 1
    transform = da.rio.transform()
    shapes = rasterio.features.shapes(values, transform=transform)

    geometries = []
    colvalues = []
    for (geom, colval) in shapes:
        if colval < 0:
            continue
        else:
            geometries.append(sg.Polygon(geom["coordinates"][0]))
            colvalues.append(colval)

    gdf = gpd.GeoDataFrame({"value": colvalues, "geometry": geometries})
    gdf.crs = da.attrs.get("crs")
    return gdf

# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# Method to convert km to decimal degrees
def km2deg(km):
    # Earth radius
    dRE = 6378.1370
    deg = 180 * km / (np.pi * dRE)
    return deg
# --------------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------





