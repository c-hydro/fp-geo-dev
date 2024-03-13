"""
Q-HyDE : FP - Hydrological Conditioning
__date__ = '20240220'
__version__ = '2.0.0'
__author__ =
        'Andrea Libertino (andrea.libertino@cimafoundation.org',
__library__ = 'fp_geo_hmc'

Version(s):
20240220 (2.0.0) --> Optimized for Continuum 3.1.6 / Restablished original alpha and beta routine
20201023 (1.1.0) --> Alpha and beta integration
20200507 (1.0.0) --> Beta release
"""

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Complete library
from qgis.core import QgsProcessing
from qgis.core import QgsProcessingAlgorithm
from qgis.core import QgsProcessingMultiStepFeedback
from qgis.core import QgsProcessingParameterRasterLayer
from qgis.core import QgsProcessingParameterNumber
from qgis.core import QgsProcessingParameterBoolean
from qgis.core import QgsProcessingParameterVectorLayer
from qgis.core import QgsProcessingParameterRasterDestination
from qgis.core import QgsProcessingParameterDefinition
from qgis.core import QgsCoordinateReferenceSystem
from qgis.core import QgsProcessingParameterPoint
from qgis.core import QgsProcessingParameterFile
from qgis.core import QgsProcessingParameterString
from qgis.core import QgsRasterLayer
from qgis.core import QgsMessageLog
from qgis.core import Qgis
from qgis.core import QgsProject
from osgeo import gdal
import numpy as np
import linecache
import copy
import os
import processing
import math
import tempfile
import shutil
from copy import deepcopy
#import rasterio as rio
import subprocess
import sys
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'fp_geo_hmc - MAKE HYDRO DERIVATIVES'
alg_version = '2.0.0'
alg_release = '2024-02-20'
# Algorithm parameter(s)
time_format = '%Y%m%d%H%M'
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
class FpMakeHydrologicalDerivatives(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
                
        # Input / output--------------------------------------------------------
        # Mandatory
        # Input DEM
        self.addParameter(QgsProcessingParameterRasterLayer('CondDEM', 'Conditioned DEM', defaultValue=None))
        # Input DEM
        self.addParameter(QgsProcessingParameterRasterLayer('RawDEM', 'Raw DEM', defaultValue=None))
        # Stream definition threshold
        self.addParameter(QgsProcessingParameterNumber('streamdefinitionkm2', 'Stream definition threshold (km2)', type=QgsProcessingParameterNumber.Double, defaultValue=0))
        # Output directory
        self.addParameter(QgsProcessingParameterFile('OutputDirectory', 'Output Directory', behavior=QgsProcessingParameterFile.Folder, fileFilter='All files (*.*)', defaultValue=None))
        # Domain name
        self.addParameter(QgsProcessingParameterString('domain', 'Domain Name', defaultValue=None))
        
        # Advanced
        # Flag to mask output
        param = QgsProcessingParameterBoolean('UseMask', 'Mask outputs at outlet section', optional=True, defaultValue=False)
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        # Outlet section for masking
        param = QgsProcessingParameterPoint('OutletSection', 'Outlet section (must lie on the channel network)', optional=True, defaultValue='0,0')
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        # ----------------------------------------------------------------------
    
    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        tempPath = tempfile.mkdtemp()
        outPath = parameters['OutputDirectory']
        domain = parameters['domain']
        feedback = QgsProcessingMultiStepFeedback(13, model_feedback)
        results = {}
        outputs = {}

        # Additional functions -------------------------------------------------
        # Harvesine formula for approximate kms to degrees   
        def kilometers2degrees(kilometer, radius=6371):
            return kilometer / (2.0 * radius * math.pi / 360.0)
        
        # Harvesine formula for approximate degrees to kms 
        def degrees2kilometers(gradi, radius=6371):
            return gradi * (2.0 * radius * math.pi / 360.0)

        # Function for mask outputs
        def maskOutput(map,mask,cellsize,region,tempPath): 
            alg_params = {
                'GRASS_RASTER_FORMAT_META': '',
                'GRASS_RASTER_FORMAT_OPT': '',
                'GRASS_REGION_CELLSIZE_PARAMETER': cellsize,
                'GRASS_REGION_PARAMETER': region,
                'a': map,
                'b': mask,
                'c': None,
                'd': None,
                'e': None,
                'expression': 'if(B==1,A,null())',
                'f': None,
                'output': os.path.join(tempPath,'temp.tif')
                }
            outputs = processing.run('grass7:r.mapcalc.simple', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
            return outputs
            
        # Function for write Continuum ASCII grids
        def convertAAIGrid(varName, map, precision, outType, outPath, domain):  
            # Write 
            alg_params = {
                'COPY_SUBDATASETS': False,
                'DATA_TYPE': 6,
                'EXTRA': '-of AAIGrid',
                'INPUT': map,
                'NODATA': -9999,
                'OPTIONS': '-co FORCE_CELLSIZE=YES -co DECIMAL_PRECISION=' + str(precision),
                'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:4326'),
                'OUTPUT': os.path.join(outPath, domain + '.' + varName + '.txt')
            }
            outputs = processing.run('gdal:translate', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
            return outputs
        
        # Function for replace nodata values
        def replaceNoData(outPath,domain,varName,outType):
            fin = open(os.path.join(outPath, domain + '.' + varName + '.txt'),"rt")
            data = fin.read()
            if outType == 'Float32':
                data = data.replace('-nan','-9999')
            if outType == 'Int16':
                data = data.replace('255','-9999')
            fin.close()
            fin = open(os.path.join(outPath, domain + '.' + varName + '.txt'),"wt")
            fin.write(data)
            fin.close()
        
        
        def read_ascii_grid(file_path):
            """
            Read an ASCII grid file and extract header information and data matrix.

            Parameters:
            - file_path (str): Path to the ASCII grid file.

            Returns:
            - header (dict): Dictionary containing header information.
            - data_matrix (numpy.ndarray): NumPy array containing the data matrix.
            """

            header = {}
            data_matrix = None
            header_lines = 0

            with open(file_path, 'r') as file:
                # Read header information
                for line in file:
                    if line.lower().startswith("ncols"):
                        header["ncols"] = int(line.split()[1])
                    elif line.lower().startswith("nrows"):
                        header["nrows"] = int(line.split()[1])
                    elif line.lower().startswith("xllcorner"):
                        header["xllcorner"] = float(line.split()[1])
                    elif line.lower().startswith("yllcorner"):
                        header["yllcorner"] = float(line.split()[1])
                    elif line.lower().startswith("cellsize"):
                        header["cellsize"] = float(line.split()[1])
                    elif line.lower().startswith("nodata_value"):
                        header["NODATA_value"] = float(line.split()[1])
                    else:
                        break
                    header_lines += 1

                # Read data matrix
                data_matrix = np.loadtxt(file_path, skiprows=header_lines)
                # set data matrix values equal to NODATA_value to np.nan
                # data_matrix[data_matrix == header["NODATA_value"]] = np.nan

            return header, data_matrix

        def write_ascii_grid(file_path, header, data_matrix):
            """
            Write header information and data matrix to an ASCII grid file.

            Parameters:
            - file_path (str): Path to the ASCII grid file.
            - header (dict): Dictionary containing header information.
            - data_matrix (numpy.ndarray): NumPy array containing the data matrix.
            """
            
            data_matrix[np.isnan(data_matrix)] = header["NODATA_value"]

            with open(file_path, 'w') as file:
                # Write header information
                for key, value in header.items():
                    file.write(f"{key} {value}\n")

                # Write data matrix
                np.savetxt(file, data_matrix, fmt="%g")

    
        def makeAlphaBeta(basePath, ancillaryPath, domain):
            feedback.setProgressText(' ---> Importing data ...')

            header, a2dDem = read_ascii_grid(os.path.join(basePath, domain + '.dem.txt'))
            _, a2iPun = read_ascii_grid(os.path.join(basePath, domain + '.pnt.txt'))
            _, a2iChoice = read_ascii_grid(os.path.join(basePath, domain + '.choice.txt'))
            _, a2dAreaCell = read_ascii_grid(os.path.join(basePath, domain + '.areacell.txt'))
            DD = 50
            a2dAlphaMap = os.path.join(basePath, domain + '.alpha.txt')
            a2dBetaMap = os.path.join(basePath, domain + '.beta.txt')

            # a2dDem = np.flipud(a2dDem)
            # a2iPun = np.flipud(a2iPun)
            a2dDem[a2iPun == -9999] = -9999
            # a2iChoice = np.flipud(a2iChoice)
            # a2dAreaCell = np.flipud(a2dAreaCell)

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

            iRows, iCols = a2dDem.shape

            # Checking distance t
            dDistanceT = 500
            if dDxM >= 100 and dDxM < 1000:
                dDistanceT = 2000
            if dDxM >= 5000 and dDxM < 20000:
                dDistanceT = 30000

            # Dem corrections
            feedback.setProgressText(' ---> Correcting elevation ...')

            a2dDem[(a2dDem <= 0) & (a2dDem > -1000)] = 0.2
            a2dDemOrig = copy.deepcopy(a2dDem)
            # a2dDem=cleanlakes(a2dDemOrig, iRows, iCols)

            # plt.imshow(a2dDemOrig-a2dDemFix)
            # plt.colorbar()
            # plt.show(block=True)

            # cleanlakes(a2dDem)
            # Define alpha matrix angle

            feedback.setProgressText(' ---> Computing alpha ...')

            for i in np.arange(0, iRows):
                for j in np.arange(0, iCols):

                    a = deepcopy(i)
                    b = deepcopy(j)

                    if a2dDem[(i, j)] > 0:
                        fNumPen = 0

                        while a2dDem[(int(a), int(b))] > 0 and diff_DD[(int(a), int(b))] == -9999:
                            if a >= 0 and a < iRows and b >= 0 and b < iCols:  # a cosa serve? ci sono dentro per forza...
                                iii = a + (int((a2iPun[(int(a), int(b))] - 1) / 3) - 1)
                                jjj = b + a2iPun[(int(a), int(b))] - 5 - 3 * (int((a2iPun[(int(a), int(b))] - 1) / 3) - 1)
                                LDD[(int(a), int(b))] = math.sqrt(((a - iii) * dDyM) ** 2 + ((b - jjj) * dDxM) ** 2)

                                if iii < 0 or jjj < 0 or iii >= iRows or jjj >= iCols:
                                    break

                                diff_DD[(int(a), int(b))] = a2dDem[(int(a), int(b))] - a2dDem[(int(iii), int(jjj))]

                                # Pendenza media sui canali
                                if math.atan2(diff_DD[(int(a), int(b))], LDD[(int(a), int(b))]) > 0 and diff_DD[(int(a), int(b))] < 9000:
                                    fNumPen = fNumPen + 1
                                    pend[(int(a), int(b))] = pend[(int(a), int(b))] + math.atan2(diff_DD[(int(a), int(b))], LDD[(int(a), int(b))])

                                while a2dDem[(int(a), int(b))] - a2dDem[
                                    (int(iii), int(jjj))] <= DD and iii >= 0 and iii < iRows - 1 and jjj >= 0 and jjj < iCols - 1 and \
                                        a2dDem[(int(iii), int(jjj))] > 0 and LDD[(int(a), int(b))] < dDistanceT:
                                    diff_DD[(int(a), int(b))] = a2dDem[(int(a), int(b))] - a2dDem[(int(iii), int(jjj))]  # Quesra non andrebbe calcolata dopo essere scesi di cella? (vedi commento sotto) MOD2023
                                    ii = iii + (int((a2iPun[(int(iii), int(jjj))] - 1) / 3) - 1)
                                    jj = jjj + a2iPun[(int(iii), int(jjj))] - 5 - 3 * (int((a2iPun[(int(iii), int(jjj))] - 1) / 3) - 1)

                                    if a2dDem[(int(a), int(b))] - a2dDem[
                                        (int(ii), int(jj))] <= DD and ii >= 0 and ii < iRows and jj >= 0 and jj < iCols:
                                        LDD[(int(a), int(b))] = LDD[(int(a), int(b))] + math.sqrt(
                                            (((ii - iii) * dDyM) ** 2 + ((jj - jjj) * dDxM) ** 2))

                                        # Pendenza media sui canali
                                        if math.atan2(diff_DD[(int(a), int(b))], LDD[(int(a), int(b))]) > 0:
                                            if a2iChoice[(int(a), int(b))] == 1 and diff_DD[(int(a), int(b))] < 9000:
                                                fNumPen = fNumPen + 1
                                                pend[(int(a), int(b))] = pend[(int(a), int(b))] + math.atan2(diff_DD[(int(a), int(b))], LDD[(int(a),
                                                                                                               int(b))])  # ma è giusta? diff_DD è riferita a (a,b)(iii,jjj) mentre LDD a (a,b)(ii,jj)

                                            if a2iChoice[(int(a), int(b))] == 0 and LDD[(int(a), int(b))] < 500 and diff_DD[
                                                (int(a), int(b))] < 9000:  # valore fissato o è residuo di dDistanceT?
                                                fNumPen = fNumPen + 1
                                                pend[(int(a), int(b))] = pend[(int(a), int(b))] + math.atan2(diff_DD[(int(a), int(b))], LDD[(int(a), int(b))])

                                    iii = deepcopy(ii)
                                    jjj = deepcopy(jj)

                                    if diff_DD[(int(iii), int(jjj))] != -9999:  # finale dei percorsi
                                        while a2dDem[(int(a), int(b))] - a2dDem[(
                                        int(iii), int(jjj))] <= DD and iii >= 0 and iii < iRows - 1 and jjj >= 0 and jjj < iCols - 1 and \
                                                a2dDem[(int(iii), int(jjj))] > 0 and LDD[(int(a), int(b))] < dDistanceT:
                                            diff_DD[(int(a), int(b))] = a2dDem[(int(a), int(b))] - a2dDem[(int(iii), int(jjj))]
                                            ii = iii + (int((a2iPun[(int(iii), int(jjj))] - 1) / 3) - 1)
                                            jj = jjj + a2iPun[(int(iii), int(jjj))] - 5 - 3 * (int((a2iPun[(int(iii), int(jjj))] - 1) / 3) - 1)

                                            if (a2dDem[(int(a), int(b))] - a2dDem[
                                                (int(ii), int(jj))]) <= DD and ii >= 0 and ii < iRows and jj >= 0 and jj < iCols:
                                                LDD[(int(a), int(b))] = LDD[(int(a), int(b))] + math.sqrt(
                                                    ((ii - iii) * dDyM) ** 2 + ((jj - jjj) * dDxM) ** 2)

                                                # Pendenza media sui canali
                                                if math.atan2(diff_DD[(int(a), int(b))], LDD[(int(a), int(b))]) > 0 and diff_DD[(int(a), int(b))] < 9000:
                                                    if a2iChoice[(int(a), int(b))] == 1:
                                                        fNumPen = fNumPen + 1
                                                        pend[(int(a), int(b))] = pend[(int(a), int(b))] + math.atan2(diff_DD[(int(a), int(b))], LDD[(int(a), int(b))])

                                                    if a2iChoice[(int(a), int(b))] == 0 and LDD[(int(a), int(b))] < 500 and diff_DD[(int(a), int(b))] < 9000:
                                                        fNumPen = fNumPen + 1
                                                        pend[(int(a), int(b))] = pend[(int(a), int(b))] + math.atan2(diff_DD[(int(a), int(b))], LDD[(int(a), int(b))])

                                            iii = deepcopy(ii)
                                            jjj = deepcopy(jj)

                                if fNumPen > 0:
                                    pend[(int(a), int(b))] = pend[(int(a), int(b))] / fNumPen

                                a2dAlpha[(int(a), int(b))] = math.atan2(DD, LDD[(int(a), int(b))])  # Angolo in radianti

                                if diff_DD[(int(a), int(b))] == 0.9 or diff_DD[(int(a), int(b))] > 500:
                                    diff_DD[(int(a), int(b))] = 0.9

                                if diff_DD[(int(a), int(b))] < 1 and LDD[(int(a), int(b))] < 4 * dDxM:
                                    LDD[(int(a), int(b))] = 4 * dDxM

                                a2dAlpha[(int(a), int(b))] = math.atan2(diff_DD[(int(a), int(b))], LDD[(int(a), int(b))])

                                ii = a + (int((a2iPun[(int(a), int(b))] - 1) / 3) - 1)
                                jj = b + a2iPun[(int(a), int(b))] - 5 - 3 * (int((a2iPun[(int(a), int(b))] - 1) / 3) - 1)

                                # if a2dAlpha[(int(a), int(b))]<0.1:
                                #    print('ciao')

                                if a2dDem[(int(ii), int(jj))] >= 0:
                                    a = deepcopy(ii)
                                    b = deepcopy(jj)
                                    fNumPen = 0
                                else:
                                    continue  # esce ma conserva gli indici della fine percorso svolto

                        # Fine di un percorso completo seguendo i puntatori
                        ii = a + (int((a2iPun[(int(a), int(b))] - 1) / 3) - 1)
                        jj = b + a2iPun[(int(a), int(b))] - 5 - 3 * (int((a2iPun[(int(a), int(b))] - 1) / 3) - 1)

            # a2dAlpha[a2dAlpha==-9999]=np.nan
            # plt.imshow(a2dAlpha)
            # plt.colorbar()
            # plt.show(block=True)

            feedback.setProgressText(' ---> Computing beta ...')
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
                                if pend2[(int(ii), int(jj))] > 0:
                                    fn = fn + 1
                                    pend[(i, j)] = pend[(i, j)] + pend2[(int(ii), int(jj))]

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

            feedback.setProgressText(' ----> Alpha max: ' + str(np.nanmax(a2dAlpha)))
            feedback.setProgressText(' ----> Alpha min: ' + str(np.nanmin(a2dAlpha[a2dAlpha != -9999])))

            feedback.setProgressText(' ----> Beta max: ' + str(np.nanmax(a2dBeta)))
            feedback.setProgressText(' ----> Beta min: ' + str(np.nanmin(a2dBeta[a2dBeta != -9999])))

            # profile = Dem_in.profile
            # profile.update(driver='GTiff', count=1)

            feedback.setProgressText(' ----> Writing alpha and beta ...')

            write_ascii_grid(a2dBetaMap, header, a2dBeta.astype(np.float32))
            write_ascii_grid(a2dAlphaMap, header, a2dAlpha.astype(np.float32))
            # with rio.open(a2dBetaMap, 'w', **profile) as dst:
            #    dst.write(np.flipud(a2dBeta.astype('Float32')), 1)

            # with rio.open(a2dAlphaMap, 'w', **profile) as dst:
            #    dst.write(np.flipud(a2dAlpha.astype('Float32')), 1)
        # ----------------------------------------------------------------------
        
        # Preliminary operations------------------------------------------------
        # Create working paths and initialize variables

        feedback.setProgressText(' ============================================================================ ')
        feedback.setProgressText(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
        feedback.setProgressText(' ==> START ... ')
        feedback.setProgressText(' ')
        
        # Export inputDEM for gdal
        alg_params = {
            'COPY_SUBDATASETS': False,
            'DATA_TYPE': 0,
            'EXTRA': '-of GTiff',
            'INPUT': parameters['RawDEM'],
            'NODATA': None,
            'TARGET_CRS': None,
            'OUTPUT': os.path.join(tempPath, 'inputDEM.tif'),
            'SILENT': True
        }
        
        temp = processing.run('gdal:translate', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        # Identify working resolution and convert threshold in number of cells 
        dem = temp['OUTPUT']
        raster = gdal.Open(dem)
        gt =raster.GetGeoTransform()
        pixelSizeX = abs(gt[1])
        pixelSizeXm = abs(degrees2kilometers(gt[1])*1000)
        Xres = degrees2kilometers(pixelSizeX)
        n_cell = int(parameters['streamdefinitionkm2']/(Xres*Xres))
        # ----------------------------------------------------------------------
        
        # Hydrological derivatives definition ----------------------------------
        # Basic watershed analysis
        feedback.setProgressText(" --> Starting watershed analysis")
        feedback.setProgressText("Basic watershed analysis")
        alg_params = {
            '-4': False,
            '-a': True,
            '-b': False,
            '-m': False,
            '-s': True,
            'GRASS_RASTER_FORMAT_META': '',
            'GRASS_RASTER_FORMAT_OPT': '',
            'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
            'GRASS_REGION_PARAMETER': dem,
            'blocking': None,
            'convergence': 6,
            'depression': None,
            'disturbed_land': None,
            'elevation': parameters['CondDEM'],
            'flow': None,
            'max_slope_length': None,
            'memory': 300,
            'threshold': n_cell,
            'accumulation': os.path.join(tempPath, 'accumulation.tif'),
            'drainage': os.path.join(tempPath,'pnt.tif'),
            'stream': os.path.join(tempPath,'stream.tif')
        }
        outputs['Rwatershed'] = processing.run('grass7:r.watershed', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(1)
        if feedback.isCanceled():
            return {}
            
        # Make output mask
        if parameters['UseMask'] is True:
            alg_params = {
                'GRASS_RASTER_FORMAT_META': '',
                'GRASS_RASTER_FORMAT_OPT': '',
                'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
                'GRASS_REGION_PARAMETER': dem,
                'coordinates': parameters['OutletSection'],
                'input': os.path.join(tempPath,'pnt.tif'),
                'output': os.path.join(tempPath,'mask.tif')
            }
            mask = processing.run('grass7:r.water.outlet', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

            feedback.setCurrentStep(2)
            if feedback.isCanceled():
                return {}
        else:
            alg_params = {
                'GRASS_RASTER_FORMAT_META': '',
                'GRASS_RASTER_FORMAT_OPT': '',
                'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
                'GRASS_REGION_PARAMETER': dem,
                'a': parameters['CondDEM'],
                'b': parameters['CondDEM'],
                'c': None,
                'd': None,
                'e': None,
                'expression': 'A/B',
                'f': None,
                'output': os.path.join(tempPath,'mask.tif')
            }
            mask = processing.run('grass7:r.mapcalc.simple', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

            feedback.setCurrentStep(2)
            if feedback.isCanceled():
                return {}
        
        # Export DEM
        varName = 'dem'
        outputs[varName] = maskOutput(dem,mask['output'],gt[1],dem,tempPath)
        results[varName] = convertAAIGrid(varName, outputs[varName]['output'],2,'Float32', outPath, domain)
        replaceNoData(outPath,domain,varName,'Float32')
        
        feedback.setCurrentStep(3)
        if feedback.isCanceled():
            return {}
            
        # Export area
        varName = 'area'
        outputs[varName] = maskOutput(outputs['Rwatershed']['accumulation'],mask['output'],gt[1],parameters['CondDEM'],tempPath)
        results[varName] = convertAAIGrid(varName, outputs[varName]['output'],0,'Float32', outPath, domain)
        replaceNoData(outPath,domain,varName,'Float32')
        
        feedback.setCurrentStep(4)
        if feedback.isCanceled():
            return {}
        
        # Compute masked choice
        varName = 'choice'
        alg_params = {
            'GRASS_RASTER_FORMAT_META': '',
            'GRASS_RASTER_FORMAT_OPT': '',
            'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
            'GRASS_REGION_PARAMETER': parameters['CondDEM'],
            'a': outputs['Rwatershed']['accumulation'],
            'b': mask['output'],
            'c': None,
            'd': None,
            'e': None,
            'expression': 'if(B==1,if(A>' + str(n_cell) + ',1,0),null())',
            'f': None,
            'output': os.path.join(tempPath, 'choice.tif')
        }
        outputs[varName] = processing.run('grass7:r.mapcalc.simple', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results[varName] = convertAAIGrid(varName, outputs[varName]['output'],0,'Int16', outPath, domain)
        replaceNoData(outPath,domain,varName,'Int16')
        
        feedback.setCurrentStep(5)
        if feedback.isCanceled():
            return {}
            
        # Compute pointers
        varName = 'pnt'
        alg_params = {
                'GRASS_RASTER_FORMAT_META': '',
                'GRASS_RASTER_FORMAT_OPT': '',
                'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
                'GRASS_REGION_PARAMETER': parameters['CondDEM'],
                'input': outputs['Rwatershed']['drainage'],
                'rules': '',
                'txtrules': '1 -1=9\n2 -2=8\n3 -3=7\n4 -4=4\n5 -5=1\n6 -6=2\n7 -7=3\n8 -8=6\n*=-9999',
                'output': os.path.join(tempPath, 'pntHMC.tif')
        }
        pnt = processing.run('grass7:r.reclass', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        outputs[varName] = maskOutput(pnt['output'],mask['output'],gt[1],parameters['CondDEM'],tempPath)
        results[varName] = convertAAIGrid(varName, outputs[varName]['output'],0,'Int16', outPath, domain)
        replaceNoData(outPath,domain,varName,'Int16')
        
        feedback.setCurrentStep(6)
        if feedback.isCanceled():
            return {}
            
        # Export mask
        varName = 'mask'
        results[varName] = convertAAIGrid(varName, mask['output'],0,'Int16', outPath, domain)
        replaceNoData(outPath,domain,varName,'Int16')
        
        feedback.setCurrentStep(7)
        if feedback.isCanceled():
            return {}
            
        # ----------------------------------------------------------------------

        # Extract geographical features ----------------------------------------
        # Compute latitude
        varName = 'lat'
        alg_params = {
            'GRASS_RASTER_FORMAT_META': '',
            'GRASS_RASTER_FORMAT_OPT': '',
            'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
            'GRASS_REGION_PARAMETER': parameters['CondDEM'],
            'a': parameters['CondDEM'],
            'b': None,
            'c': None,
            'd': None,
            'e': None,
            'expression': 'y()',
            'f': None,
            'output': os.path.join(tempPath, 'lat.tif')
        }
        outputs[varName] = processing.run('grass7:r.mapcalc.simple', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        # Export latitude
        results[varName] = convertAAIGrid(varName, outputs[varName]['output'],4,'Float32', outPath, domain)
        
        feedback.setCurrentStep(8)
        if feedback.isCanceled():
            return {}
            
        # Compute longitude
        varName = 'lon'
        alg_params = {
            'GRASS_RASTER_FORMAT_META': '',
            'GRASS_RASTER_FORMAT_OPT': '',
            'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
            'GRASS_REGION_PARAMETER': parameters['CondDEM'],
            'a': parameters['CondDEM'],
            'b': None,
            'c': None,
            'd': None,
            'e': None,
            'expression': 'x()',
            'f': None,
            'output': os.path.join(tempPath, 'lon.tif')
        }
        outputs[varName] = processing.run('grass7:r.mapcalc.simple', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        # Export longitude
        results[varName] = convertAAIGrid(varName, outputs[varName]['output'],4,'Float32', outPath, domain)
        
        feedback.setCurrentStep(9)
        if feedback.isCanceled():
            return {}        
        
        # Compute Area Cell
        varName = 'areacell'
        alg_params = {
            'GRASS_RASTER_FORMAT_META': '',
            'GRASS_RASTER_FORMAT_OPT': '',
            'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
            'GRASS_REGION_PARAMETER': parameters['CondDEM'],
            'a': os.path.join(tempPath, 'lon.tif'),
            'b': None,
            'c': None,
            'd': None,
            'e': None,
            'expression': 'area()',
            'f': None,
            'output': os.path.join(tempPath, 'areacell.tif')
        }
        outputs[varName] = processing.run('grass7:r.mapcalc.simple', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(10)
        if feedback.isCanceled():
            return {}  
            
        # Export areacell
        results[varName] = convertAAIGrid(varName, outputs[varName]['output'],4,'Float32', outPath, domain)
        
        # ----------------------------------------------------------------------
        
        # Alpha and beta -------------------------------------------------------
        # Sink drainage route detection
        alg_params = {
            'ELEVATION': parameters['CondDEM'],
            'THRESHOLD       ': False,
            'THRSHEIGHT': 100,
            'SINKROUTE': os.path.join(tempPath, 'sink_route.sdat')
        }
        outputs['SinkDrainageRouteDetection'] = processing.run('saga:sinkdrainageroutedetection', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(11)
        if feedback.isCanceled():
            return {}

        # Sink removal
        alg_params = {
            'DEM': parameters['CondDEM'],
            'METHOD': 1,
            'SINKROUTE': outputs['SinkDrainageRouteDetection']['SINKROUTE'],
            'THRESHOLD        ': False,
            'THRSHEIGHT': 100,
            'DEM_PREPROC': os.path.join(tempPath, 'DEM_preproc.sdat')
        }
        outputs['SinkRemoval'] = processing.run('saga:sinkremoval', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(12)
        if feedback.isCanceled():
            return {}
            
        makeAlphaBeta(outPath, tempPath, domain)
        
        feedback.setCurrentStep(13)
        if feedback.isCanceled():
            return {}
        # ----------------------------------------------------------------------
        
        # Cleaning system ------------------------------------------------------        
        # Delete temporary folder
        shutil.rmtree(tempPath)
        [os.remove(os.path.join(outPath,f)) for f in os.listdir(outPath) if f.endswith(".xml")]
        [os.remove(os.path.join(outPath,f)) for f in os.listdir(outPath) if f.endswith(".prj")]
    
        #try:
        #    [os.remove(os.path.join(outPath,f)) for f in os.listdir(outPath) if f.endswith(".tfw")]
        #except:
        #    pass
        # ----------------------------------------------------------------------
        
        return results

    def name(self):
        return 'Make hydro derivatives'

    def displayName(self):
        return 'Make hydro derivatives'

    def group(self):
        return 'fp_geo_hmc'

    def groupId(self):
        return 'fp_geo_hmc'

    def createInstance(self):
        return FpMakeHydrologicalDerivatives()
