#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fp-geo - Make Dams and Lakes file

__date__ = '20220701'
__version__ = '1.0.2'
__author__ =
        'Andrea Libertino (andrea.libertino@cimafoundation.org')
__library__ = 'fp-geo'

General command line:
python3 MakeSoilMaps

Version(s):
20220701 (1.0.2) -->    Correct lake routine bug
20220607 (1.0.1) -->    Include json config file
20220113 (1.0.0) --> 	Beta release
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
alg_release = '2022-07-01'
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

    output_file_dams = os.path.join(out_folder,data_settings["data"]["output"]["filename_dams"]).format(domain=data_settings["algorithm"]["domain"])
    output_file_lakes = os.path.join(out_folder,data_settings["data"]["output"]["filename_lakes"]).format(domain=data_settings["algorithm"]["domain"])

    input_file_dams = data_settings["data"]["input"]["dams"].replace("{domain}",data_settings["algorithm"]["domain"])
    input_file_lakes = data_settings["data"]["input"]["lakes"].replace("{domain}",data_settings["algorithm"]["domain"])
    input_file_choice = data_settings["data"]["input"]["choice"].replace("{domain}",data_settings["algorithm"]["domain"])
    input_file_pnt = data_settings["data"]["input"]["pnt"].replace("{domain}",data_settings["algorithm"]["domain"])

    # Other initialization variables
    dam_sep = '##################################################################################################\n'
    lake_sep = '##################################################################################################\n'
    lakes_not_valid = []

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
    pnt = rxr.open_rasterio(input_file_pnt).squeeze().values

    res = np.abs(choice_map.rio.resolution()[0])
    x_ll = np.min(choice_map.x.values) - (res/2)
    y_ul = np.max(choice_map.y.values) + (res/2)

    # -------------------------------------------------------------------------------------
    # Dams
    logging.info(" --> Compute dams...")
    logging.info(" ---> Load input dams file...")
    if os.path.isfile(input_file_dams):
        gdf = gpd.read_file(input_file_dams)
        logging.info(" ---> Load input dams file...DONE")

        unnamed_dam_code = 1
        index = []
        for name in gdf["DAM_NAME"]:
            if name is not None:
                index = index + [name]
            else:
                index = index + ["Dam_" + data_settings["algorithm"]["domain"] + str(unnamed_dam_code).zfill(3)]
                unnamed_dam_code = unnamed_dam_code + 1

        gdf["DAM_NAME"] = index
        gdf = gdf.set_index('DAM_NAME')

        logging.info(" ---> Check input dams file...")
        logging.info(" ----> " + str(len(gdf)) + " dams found!")

        with open(output_file_dams ,"w") as dam_file:
            dam_file.write(str(len(gdf)) + "\t" + "#Number of dams\n")
            dam_file.write(str(len(gdf)) + "\t" + "#Number of plants\n")
            dam_file.write(dam_sep)

            for dam_name in gdf.index.values:
                # Dam name
                logging.info('----> Characterize dam: ' + dam_name + '...')
                dam_file.write(dam_name + "\t\t\t#Dam name " + str(gdf.loc[dam_name,"YEAR"]) + "\n")
                # Dam coords
                Y_HMC = np.ceil(np.abs((gdf.loc[dam_name].geometry.x-x_ll))/res)
                X_HMC = np.ceil(np.abs((gdf.loc[dam_name].geometry.y-y_ul))/res)
                if choice[int(X_HMC) - 1, int(Y_HMC) - 1]==0:
                    raise ValueError("DOMAIN: " + data_settings["algorithm"]["domain"] + ".The dam " + dam_name + " is not on the network")
                dam_file.write(str(int(X_HMC)) + " " + str(int(Y_HMC)) + "\t\t\t#Row and column dam coordinates\n")
                # Number of plants
                dam_file.write("1\t\t\t#Number of plants downstream the dam\n")
                # Code for distributed reservoirs
                dam_file.write("-9999\t\t\t#Code of the reservoirs cells of the dam (if point dam set to -9999)\n")
                # Storage features
                max_storage_m3 = gdf.loc[dam_name,"CAP_MCM"]*(10**6)
                initial_storage_m3 = max_storage_m3 * 0.7
                dam_file.write(str(max_storage_m3) + "\t\t\t#Max storage (m3)\n")
                dam_file.write(str(initial_storage_m3) + "\t\t\t#Initial storage (m3)\n")
                # Spillway features
                dam_file.write("99999\t\t\t#Critical discharge for surface spillway\n")
                dam_len = gdf.loc[dam_name, "DAM_LEN_M"]
                if dam_len > 0:
                    spillway_len = 0.15 * dam_len
                else:
                    spillway_len = dam_len
                dam_file.write(str(spillway_len) + "\t\t\t#Equivalent length of surface spillway\n")
                # Maximum reservoir depth
                max_depth = gdf.loc[dam_name, "DAM_HGT_M"]
                if max_depth < 0:
                    max_depth = np.round(2 * gdf.loc[dam_name, "CAP_MCM"] / gdf.loc[dam_name, "AREA_SKM"],0)
                dam_file.write(str(max_depth) + "\t\t\t#Maximum reservoir depth\n")
                # Linear tank coefficient
                dam_file.write("1e-006\t\t\t#Linear tank coefficient\n")
                # Ancillary files name
                dam_file.write("\t\t\t\t#Depth-volume curve file name\n")
                dam_file.write("\t\t\t\t#Turbines discharge file name\n")
                # Downstream cell
                X_HMC_OUT = int(X_HMC) - (int((pnt[int(X_HMC)-1, int(Y_HMC)-1] - 1) / 3) - 1)
                Y_HMC_OUT = int(Y_HMC) + pnt[int(X_HMC)-1, int(Y_HMC)-1] - 5 - 3 * (int((pnt[int(X_HMC)-1, int(Y_HMC)-1] - 1) / 3) - 1)
                if choice[int(X_HMC_OUT) - 1, int(Y_HMC_OUT) - 1]==0:
                    raise ValueError("DOMAIN: " + data_settings["algorithm"]["domain"] + ".The dam " + dam_name + " output is not on the network")
                dam_file.write(str(int(X_HMC_OUT)) + " " + str(int(Y_HMC_OUT)) + "\t\t\t#Row and column outlet dam coordinates\n")
                # Plant features
                max_discharge_m3_s =  0.001 * 2.5 * gdf.loc[dam_name, "DIS_AVG_LS"]
                dam_file.write("-9999\t\t\t#Plant corrivation time (minutes)\n")
                dam_file.write(str(max_discharge_m3_s) + "\t\t\t#Maximum plant discharge (m3/s)\n")
                dam_file.write("1\t\t\t#flag=1 if the plant discharge water\n")
                dam_file.write(dam_sep)
                logging.info('----> Characterize dam: ' + dam_name + '...DONE')
    else:
        logging.warning(" --> WARNING! No dams input file found for the domain!")
    logging.info(" --> Compute dams...DONE!")

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Lakes
    logging.info(" --> Compute lakes...")
    logging.info(" ---> Load input lakes file...")

    if os.path.isfile(input_file_lakes):
        gdf = gpd.read_file(input_file_lakes)
        logging.info(" ---> Load input lakes file...DONE")

        unnamed_lake_code = 1
        index = []
        for name in gdf["Lake_name"]:
            if name is not None:
                index = index + [name]
            else:
                index = index + ["Lake_" + data_settings["algorithm"]["domain"] + str(unnamed_lake_code).zfill(3)]
                unnamed_lake_code = unnamed_lake_code + 1

        gdf["Lake_name"] = index
        gdf = gdf.set_index('Lake_name')

        logging.info(" ---> Check input lakes file...")
        logging.info(" ---> Drop lakes with invalid residence time...")
        # Drops lake if residence time is not valid
        index_to_drops = gdf.index[gdf['Res_time'] < 0]
        if len(index_to_drops) > 0:
            lakes_not_valid = lakes_not_valid + [data_settings["algorithm"]["domain"] + " : " + i for i in index_to_drops]
            gdf.drop(gdf.index[gdf['Res_time'] < 0], inplace=True)

        logging.info(" ----> " + str(len(gdf)) + " valid lakes found!")

        with open(output_file_lakes, "w") as dam_file:
            dam_file.write(str(len(gdf)) + "\t" + "#Number of lakes\n")
            dam_file.write(lake_sep)

            for lake_name in gdf.index.values:
                res_time_day = gdf.loc[lake_name, "Res_time"]
                # Lake name
                logging.info('----> Characterize lake: ' + lake_name)
                dam_file.write(lake_name + "\t\t\t#Lake name\n")
                # Lake coords
                Y_HMC = np.ceil(np.abs((gdf.loc[lake_name].geometry.x-x_ll))/res)
                X_HMC = np.ceil(np.abs((gdf.loc[lake_name].geometry.y-y_ul))/res)
                if choice[int(X_HMC) - 1, int(Y_HMC) - 1] == 0:
                    raise ValueError("DOMAIN: " + data_settings["algorithm"]["domain"] + ".The lake " + lake_name + " " + str(X_HMC) + "-" + str(Y_HMC) + " is not on the network")
                dam_file.write(str(int(X_HMC)) + " " + str(int(Y_HMC)) + "\t\t\t#Row and column dam coordinates\n")
                # Code for distributed lakes
                dam_file.write("-9999\t\t\t#Code of the lakes cells of the lake (if point dam set to -9999)\n")
                # Volume features
                vol_tot = gdf.loc[lake_name, "Vol_total"] * (10 ** 6)
                dis_avg = gdf.loc[lake_name, "Dis_avg"]
                vol_min = 0
                vol_init = vol_tot
                # vol_min = vol_avg - (dis_avg * res_time_day * 24 * 3600)
                # dam_file.write(str(vol_min) + "\t\t\t#Minimum storage non-null discharge (m3)\n")
                # dam_file.write(str(vol_init) + "\t\t\t#Initial storage (m3)\n")
                # Discharge features
                res_time_hr = res_time_day * 24
                lake_const = 1 / res_time_hr
                # discharge_exp = gdf.loc[lake_name,"Vol_total"]*(10**6)/(res_time_hr*3600)
                dam_file.write(str(vol_min) + "\t\t\t#Minimum storage non-null discharge (m3)\n")
                dam_file.write(str(vol_init) + "\t\t\t#Initial storage (m3)\n")
                dam_file.write(str(lake_const) + "\t\t\t#Lake constant (1/h) \n")
                dam_file.write(lake_sep)
                logging.info('----> Characterize lake: ' + lake_name + '...DONE')

        if len(lakes_not_valid) > 0:
            logging.warning(" --> WARNING! Lakes skipped because residence time is missing: " + ";\n".join(lakes_not_valid))
    else:
        logging.warning(" --> WARNING! No lakes input file found for the domain!")
    logging.info(" --> Compute lakes...DONE!")

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