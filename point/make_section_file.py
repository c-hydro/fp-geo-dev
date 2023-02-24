#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fp-geo - Make Section file

__date__ = '20220607'
__version__ = '1.0.0'
__author__ =
        'Andrea Libertino (andrea.libertino@cimafoundation.org')
__library__ = 'fp-geo'

General command line:
python3 MakeSoilMaps

Version(s):
20220607 (1.0.0) -->    Beta release
"""

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Complete library
import geopandas as gpd
import rioxarray as rxr
import numpy as np
import time, json, logging, os
from argparse import ArgumentParser

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'fp-geo - Make dams and lakes file'
alg_version = '1.0.1'
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

    # Initialize folders and files
    logging.info(" --> Initialize folders and files...")
    out_folder = data_settings["data"]["output"]["folder"].replace("{domain}",data_settings["algorithm"]["domain"])
    os.makedirs(out_folder, exist_ok=True)

    output_file = os.path.join(out_folder,data_settings["data"]["output"]["filename"]).format(domain=data_settings["algorithm"]["domain"])

    input_file_section = data_settings["data"]["input"]["section"]["filename"].replace("{domain}",data_settings["algorithm"]["domain"])
    input_file_choice = data_settings["data"]["input"]["choice"].replace("{domain}",data_settings["algorithm"]["domain"])

    logging.info(" --> Initialize folders and files...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Initialize algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> START ... ')
    logging.info(' ==> ALGORITHM SETTINGS <== ')
    logging.info(' ----> Domain: ' + data_settings["algorithm"]["domain"])
    # -------------------------------------------------------------------------------------

    choice_map = rxr.open_rasterio(input_file_choice)
    choice = choice_map.squeeze().values

    res = np.abs(choice_map.rio.resolution()[0])
    x_ll = np.min(choice_map.x.values) - (res/2)
    y_ul = np.max(choice_map.y.values) + (res/2)

    # -------------------------------------------------------------------------------------
    # Compute
    logging.info(" --> Compute sections...")
    logging.info(" ---> Load input section file...")

    if os.path.isfile(input_file_section):
        import fiona
        try:
            gdf = gpd.read_file(input_file_section)
        except RecursionError:
            gdf = gpd.read_file(input_file_section, ignore_fields=["geometry"])
        logging.info(" ---> Load input dams file...DONE")

        gdf_col_name = data_settings["data"]["input"]["section"]["name_col"]
        gdf_col_river = data_settings["data"]["input"]["section"]["river_col"]

        unnamed_dam_code = 1
        index = []
        for name in gdf[gdf_col_name]:
            if name is not None:
                index = index + [name.replace(" ","_")]
            else:
                index = index + ["Section_" + domain + str(unnamed_dam_code).zfill(3)]
                unnamed_dam_code = unnamed_dam_code + 1

        gdf[gdf_col_name] = index
        gdf = gdf.set_index(gdf_col_name)

        logging.info(" ---> Check input sections file...")
        logging.info(" ----> " + str(len(gdf)) + " sections found!")

        with open(output_file ,"w") as section_file:
            for sect_name in gdf.index.values:
                # Sect name and river
                logging.info('----> Characterize section: ' + sect_name + '...')
                if gdf_col_river is None:
                    river = domain
                else:
                    try:
                        river = gdf.loc[sect_name, gdf_col_river].replace(" ","_")
                    except:
                        river = "-9999"
                # Sect coords
                Y_HMC = np.ceil(np.abs((gdf.loc[sect_name].geometry.x-x_ll))/res)
                X_HMC = np.ceil(np.abs((gdf.loc[sect_name].geometry.y-y_ul))/res)
                if choice[int(X_HMC) - 1, int(Y_HMC) - 1]==0:
                    raise ValueError("DOMAIN: " + domain + ".The section " + sect_name + " is not on the network")
                section_file.write(str(int(X_HMC)) + " " + str(int(Y_HMC)) + " " + river + " " + sect_name + "\n")

                logging.info('----> Characterize section: ' + sect_name + '...DONE')
    else:
        logging.warning(" --> WARNING! No sections input file found for the domain!")
    logging.info(" --> Compute sections...DONE!")

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
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