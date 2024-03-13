"""
Q-HyDE : FP - Hydrological Conditioning
__date__ = '20230220'
__version__ = '2.0.0'
__author__ =
        'Andrea Libertino (andrea.libertino@cimafoundation.org',
__library__ = 'fp_geo_hmc'

Version(s):
20230220 (2.0.0) --> Optimized for Continuum 3.1.6
20201023 (1.0.0) --> Beta release
"""
# ------------------------------------------------------------------------------

from qgis.core import QgsProcessing
from qgis.core import QgsProcessingAlgorithm
from qgis.core import QgsProcessingMultiStepFeedback
from qgis.core import QgsProcessingParameterRasterLayer
from qgis.core import QgsProcessingParameterVectorLayer
from qgis.core import QgsProcessingParameterNumber
from qgis.core import QgsProcessingParameterBoolean
from qgis.core import QgsProcessingParameterRasterDestination
from qgis.core import QgsProcessingParameterDefinition
from qgis.core import QgsProcessingParameterFile
from qgis.core import QgsProcessingParameterString
from qgis.core import QgsRasterLayer
from qgis.core import QgsProject
from osgeo import gdal
import processing
import tempfile
import os
import math
import shutil


class FpHydrologicalConditioning(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        
        # Input / output -------------------------------------------------------
        # Mandatory
        # Input DEM
        self.addParameter(QgsProcessingParameterRasterLayer('DEM', 'DEM', defaultValue=None))
        # Stream definition threshold
        self.addParameter(QgsProcessingParameterNumber('streamdefinitionkm2', 'Stream definition threshold (km2)', type=QgsProcessingParameterNumber.Double, defaultValue=0))
        # Output directory
        self.addParameter(QgsProcessingParameterFile('OutputDirectory', 'Output Directory', behavior=QgsProcessingParameterFile.Folder, fileFilter='All files (*.*)', defaultValue=None))
        # Domain name
        self.addParameter(QgsProcessingParameterString('domain', 'Domain Name', defaultValue=None))
        
        # Additional
        # Flag to resample the DEM to the output resolution
        param = QgsProcessingParameterBoolean('ResampleDEM', 'Resample DEM', optional=True, defaultValue=False)
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        #Output resolution
        param = QgsProcessingParameterNumber('Outputresolutionkm', 'Output resolution (km)', optional=True, type=QgsProcessingParameterNumber.Double, minValue=-1.79769e+308, maxValue=1.79769e+308, defaultValue=None)
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        # Flag to carve the DEM
        param = QgsProcessingParameterBoolean('CarveDEM', 'Carve DEM', optional=True, defaultValue=False)
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        # Specify carve depth
        param = QgsProcessingParameterNumber('Carvingdepthm', 'Carving depth (m)', optional=True, type=QgsProcessingParameterNumber.Double, defaultValue=None)
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        # Carving bluelines
        param = QgsProcessingParameterVectorLayer('Bluelines', 'Blue lines', optional=True, types=[QgsProcessing.TypeVectorLine], defaultValue=None)
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        # Flag to wall the DEM
        param = QgsProcessingParameterBoolean('WallDEM', 'Use walls', optional=True, defaultValue=False)
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        # Walls shapefile
        param = QgsProcessingParameterVectorLayer('Walls', 'Walls', optional=True, types=[QgsProcessing.TypeVectorLine], defaultValue=None)
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        # ----------------------------------------------------------------------

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        
        # Additional functions -------------------------------------------------
        # Harvesine formula for approximate degrees to kms and viceversa  
        def kilometers2degrees(kilometer, radius=6371):
            return kilometer / (2.0 * radius * math.pi / 360.0)
        
        def degrees2kilometers(gradi, radius=6371):
            return gradi * (2.0 * radius * math.pi / 360.0)
    
        # Create working paths and initialize variables      
        tempPath = tempfile.mkdtemp()
        outPath = parameters['OutputDirectory']
        domain = parameters['domain']
        feedback = QgsProcessingMultiStepFeedback(9, model_feedback) 
        results = {}
        outputs = {}
        # ----------------------------------------------------------------------

        
        # Resample DEM ---------------------------------------------------------
        if parameters['ResampleDEM'] is True:
        
            # Convert output resolution in degree
            ResOut = kilometers2degrees(parameters['Outputresolutionkm'])
            
            # Reproject 
            alg_params = {
                'DATA_TYPE': 0,
                'EXTRA': '',
                'INPUT': parameters['DEM'],
                'MULTITHREADING': False,
                'NODATA': -9999,
                'OPTIONS': '',
                'RESAMPLING': 2,
                'SOURCE_CRS': None,
                'TARGET_CRS': None,
                'TARGET_EXTENT': None,
                'TARGET_EXTENT_CRS': None,
                'TARGET_RESOLUTION': ResOut,
                'OUTPUT': os.path.join(tempPath,'DEMTemp.tif')
            }
            outputs['Reprojectfunc'] = processing.run('gdal:warpreproject', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
                        
            feedback.setCurrentStep(1)
            if feedback.isCanceled():
                return {}
        
        else:
            
            # Import DEM as it
            alg_params = {
            'GRASS_RASTER_FORMAT_META': '',
            'GRASS_RASTER_FORMAT_OPT': '',
            'GRASS_REGION_CELLSIZE_PARAMETER': 0,
            'GRASS_REGION_PARAMETER': '',
            'a': parameters['DEM'],
            'b': None,
            'c': None,
            'd': None,
            'e': None,
            'expression': 'A*1',
            'f': None,
            'output': os.path.join(tempPath,'DEMTemp.tif')
            }
            outputs['Reprojectfunc'] = processing.run('grass7:r.mapcalc.simple', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
            
            feedback.setCurrentStep(1)
            if feedback.isCanceled():
                return {}
        
        # Save copy of the dem pre-hydrological conditioning
        os.makedirs(outPath, exist_ok=True)
        shutil.copyfile(os.path.join(tempPath,'DEMTemp.tif'), os.path.join(outPath, domain + '_DEMraw.tif'))

        # Identify working resolution and convert threshold in number of cells 
        raster = gdal.Open(os.path.join(tempPath,'DEMTemp.tif'))
        gt =raster.GetGeoTransform()
        pixelSizeX = abs(gt[1])
        Xres = degrees2kilometers(pixelSizeX)
        n_cell = int(parameters['streamdefinitionkm2']/(Xres*Xres))
        # ----------------------------------------------------------------------
        
        # Carving streams ------------------------------------------------------
        if parameters['CarveDEM'] is True:
            
            # rasterize BlueLines
            alg_params = {
                'GRASS_MIN_AREA_PARAMETER': 0.0001,
                'GRASS_RASTER_FORMAT_META': '',
                'GRASS_RASTER_FORMAT_OPT': '',
                'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
                'GRASS_REGION_PARAMETER': os.path.join(tempPath,'DEMTemp.tif'),
                'GRASS_SNAP_TOLERANCE_PARAMETER': -1,
                'attribute_column': '',
                'input': parameters['Bluelines'],
                'label_column': '',
                'memory': 300,
                'rgb_column': '',
                'type': [1],
                'use': 2,
                'value': 1,
                'where': '',
                'output': os.path.join(tempPath,'BlueLinesTemp.tif')
            }
            outputs['Rasterizebluelinesfunc'] = processing.run('grass7:v.to.rast', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

            feedback.setCurrentStep(4)
            if feedback.isCanceled():
                return {}

            # carving
            alg_params = {
                'GRASS_RASTER_FORMAT_META': '',
                'GRASS_RASTER_FORMAT_OPT': '',
                'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
                'GRASS_REGION_PARAMETER': os.path.join(tempPath,'DEMTemp.tif'),
                'a': os.path.join(tempPath,'DEMTemp.tif'),
                'b': os.path.join(tempPath,'BlueLinesTemp.tif'),
                'c': None,
                'd': None,
                'e': None,
                'expression': 'if(isnull(B),A,A-B*' +str(parameters['Carvingdepthm']) +')',
                'f': None,
                'output': os.path.join(tempPath,'DEMTemp.tif')
            }
            outputs['carvingFunc'] = processing.run('grass7:r.mapcalc.simple', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
            
            feedback.setCurrentStep(5)
            if feedback.isCanceled():
                return {}
            
        else:
                        
            feedback.setCurrentStep(5)
            if feedback.isCanceled():
                return {}
        # ----------------------------------------------------------------------
        
        # Insert walls ---------------------------------------------------------

        if parameters['WallDEM'] is True:    

            # rasterize walls
            alg_params = {
                'GRASS_MIN_AREA_PARAMETER': 0.0001,
                'GRASS_RASTER_FORMAT_META': '',
                'GRASS_RASTER_FORMAT_OPT': '',
                'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
                'GRASS_REGION_PARAMETER': os.path.join(tempPath,'DEMTemp.tif'),
                'GRASS_SNAP_TOLERANCE_PARAMETER': -1,
                'attribute_column': '',
                'input': parameters['Walls'],
                'label_column': '',
                'memory': 300,
                'rgb_column': '',
                'type': [1],
                'use': 2,
                'value': 5000,
                'where': '',
                'output': os.path.join(tempPath,'WallsTemp.tif')
            }
            outputs['Rasterizewalls'] = processing.run('grass7:v.to.rast', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

            feedback.setCurrentStep(6)
            if feedback.isCanceled():
                return {}

            # insert walls
            alg_params = {
                'GRASS_RASTER_FORMAT_META': '',
                'GRASS_RASTER_FORMAT_OPT': '',
                'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
                'GRASS_REGION_PARAMETER': os.path.join(tempPath,'DEMTemp.tif'),
                'a': os.path.join(tempPath,'DEMTemp.tif'),
                'b': os.path.join(tempPath,'WallsTemp.tif'),
                'c': None,
                'd': None,
                'e': None,
                'expression': 'if(isnull(B),A,A+B)',
                'f': None,
                'output': os.path.join(tempPath,'DEMTemp.tif')
            }
            outputs['Wallingfunc'] = processing.run('grass7:r.mapcalc.simple', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
            feedback.setCurrentStep(7)
            if feedback.isCanceled():
                return {}
                
        else:
            
            feedback.setCurrentStep(7)
            if feedback.isCanceled():
                return {}
        # ----------------------------------------------------------------------
        
        # Finalize outputs -----------------------------------------------------
        # Convert and save DEM
        alg_params = {
            'COPY_SUBDATASETS': False,
            'DATA_TYPE': 0,
            'EXTRA': '-of GTiff',
            'INPUT': os.path.join(tempPath,'DEMTemp.tif'),
            'NODATA': None,
            'OPTIONS': '',
            'TARGET_CRS': None,
            'OUTPUT': os.path.join(outPath, domain + '_DEMconditioned.tif')
        }
        
        outputs['TranslateConvertFormat'] = processing.run('gdal:translate', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['Output'] = outputs['TranslateConvertFormat']['OUTPUT']
        
        feedback.setCurrentStep(8)
        if feedback.isCanceled():
            return {}
        
        # Extract streams
        alg_params = {
            '-4': False,
            '-a': True,
            '-b': False,
            '-m': False,
            '-s': True,
            'GRASS_RASTER_FORMAT_META': '',
            'GRASS_RASTER_FORMAT_OPT': '',
            'GRASS_REGION_CELLSIZE_PARAMETER': gt[1],
            'GRASS_REGION_PARAMETER': outputs['TranslateConvertFormat']['OUTPUT'],
            'blocking': None,
            'convergence': 6,
            'depression': None,
            'disturbed_land': None,
            'elevation': outputs['TranslateConvertFormat']['OUTPUT'],
            'flow': None,
            'max_slope_length': None,
            'memory': 300,
            'threshold': n_cell,
            'accumulation': os.path.join(tempPath, 'accumulation.tif'),
            'drainage': os.path.join(tempPath,'pnt.tif'),
            'stream': os.path.join(tempPath,'stream.tif')
        }
        outputs['Rwatershed'] = processing.run('grass7:r.watershed', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        feedback.setCurrentStep(9)
        if feedback.isCanceled():
            return {}
        
        # Convert and save streams
        alg_params = {
            'COPY_SUBDATASETS': False,
            'DATA_TYPE': 0,
            'EXTRA': '-of GTiff',
            'INPUT': os.path.join(tempPath,'stream.tif'),
            'NODATA': None,
            'OPTIONS': '',
            'TARGET_CRS': None,
            'OUTPUT': os.path.join(outPath, domain + '_streams.tif')
        }
        
        outputs['TranslateConvertFormat'] = processing.run('gdal:translate', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        results['OutputStream'] = {'stream':outputs['TranslateConvertFormat']['OUTPUT']}
        # ----------------------------------------------------------------------
        
        # Cleaning system -----------------------------------------------------        
        shutil.rmtree(tempPath)
        
        return results

    def name(self):
        return 'Land preprocessing'

    def displayName(self):
        return 'Land preprocessing'

    def group(self):
        return 'fp_geo_hmc'

    def groupId(self):
        return 'fp_geo_hmc'

    def createInstance(self):
        return FpHydrologicalConditioning()
