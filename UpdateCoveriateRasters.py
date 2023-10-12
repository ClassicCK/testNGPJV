#NGP Tool Script #1
#Covariate Raster Generator
#8/4/2022
        # Tool inputs:       
                    #        RAP HAGB Layer
                    #        RAP Veg Cover Layer
                    #        NASS Layer
                    #        Included Model Layers - should be a gdb that contains the following 1k rasters:
                                # HAGB Average
                                # Wetland Area
                                # Wetland Count
                                # Lat
                                # Lon
                                # VRM
                                # ESD High Med Low
                                # Road length
                                # UTM grid to standardize everything
                                
#Scratch Workspace / Environment Settings
import arcpy, os, math, time
from arcpy import env
from arcpy.sa import *

#p = arcpy.mp.ArcGISProject("CURRENT")

#scratch_workspace = str(p.defaultGeodatabase)

arcpy.env.overwriteOutput = True

    #Set Environment Coordinate System
Output_CS = 'PROJCS["NAD_1983_UTM_Zone_13N",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["Meter",1]]'
arcpy.env.outputCoordinateSystem = Output_CS

#Inputs: 
#        Model Layers (user defined, geodatabase containing the input layers):

ModelLayersGDB = arcpy.GetParameterAsText(0)

#ModelLayersGDB = r"C:\Users\eric.chabot\Documents\NGPJV_spatial_tool\Test_Data_acquisition_tool\ModelLayers.gdb"
# For testing

arcpy.management.Copy(ModelLayersGDB, arcpy.Describe(ModelLayersGDB).path+"\ModelLayers_updated", "Workspace", None)
scratch_workspace = str(arcpy.Describe(ModelLayersGDB).path + "\ModelLayers_updated.gdb")
study_area_z13 = scratch_workspace+"\Study_area_z13"
HAGBMean1k = ModelLayersGDB+"\HAGBMean1k"
WetlandArea = ModelLayersGDB+"\WetlandAreaM2_1k"
WetlandCount = ModelLayersGDB+"\WetlandCount1k"
Lat = ModelLayersGDB+"\Lat"
Lon = ModelLayersGDB+"\Lon"
VRM = ModelLayersGDB+"\VRM1k"
Roads2010 = ModelLayersGDB+"\Roads_2010_1k"
Roads2021 = ModelLayersGDB+"\Roads_2021_1k"
ESDPolys = ModelLayersGDB+"\ESD_NGP"
ESD_prod_class = ModelLayersGDB+"\ESD_Class"
UTMGrid = ModelLayersGDB+"\Grid_1k_NGP"
CropReclass = ModelLayersGDB+"\CDL_reclass"
DevCropH20 = ModelLayersGDB+"\DevCropH20"



#        Study area polygon (User defined)

#study_area = r"C:\Users\eric.chabot\Documents\NGPJV_spatial_tool\Test_Data_acquisition_tool\TestRun06172022\Test_Study_Area\StudyArea.shp"
# For testing

study_area = arcpy.GetParameterAsText(1)




#        RAP Veg Cover layer (from 'Acquire New Data' Script)

#VegCover = r"C:\Users\eric.chabot\Documents\NGPJV_spatial_tool\Test_Data_acquisition_tool\TestRun06172022\RAP_VegCover_2021.tif"
# For testing

VegCover = arcpy.GetParameterAsText(2)
    # This comes as a multi-band raster
    
    
    
    
#        RAP HAGB layer (from 'Acquire New Data' Script)

#HAGB = r"C:\Users\eric.chabot\Documents\NGPJV_spatial_tool\Test_Data_acquisition_tool\TestRun06172022\RAP_HAGB_2021.tif"
# For testing

HAGB = arcpy.GetParameterAsText(3)



#        NASS layers (from 'Acquire New Data' Script)

#NASS = r"C:\Users\eric.chabot\Documents\NGPJV_spatial_tool\Test_Data_acquisition_tool\TestRun06172022\NASS_2021.tif"
# For testing

NASS = arcpy.GetParameterAsText(4)
Year = str(VegCover)[-8:-4]


# Define some more paths for later use:
study_area_grid = scratch_workspace+"/study_area_grid"
Crop9k = scratch_workspace+"\Crop9k_1k"+Year
CropRaster = scratch_workspace+"\Cropland"
RapTemp = scratch_workspace+"\RapTemp"
CropMask = scratch_workspace+"\CropMask"
NASS_null = scratch_workspace+"\SetNull_NASS"
AGFCCover = scratch_workspace+"\AGFC"
LITRCover = scratch_workspace+"\LITR"
PGFCCover = scratch_workspace+"\PGFC"
SHRBCover = scratch_workspace+"\SHRB"
TREECover = scratch_workspace+"\TREE"
Rangeland = scratch_workspace+"\\Rangeland"
HAGB_masked = scratch_workspace+"\HAGB_masked"
SnapRaster = scratch_workspace+"/Snap_Raster_1k"
zonaltable = scratch_workspace+"/zonaltable"
MeanPGFC = scratch_workspace+"/PGFC_1k"+Year
MeanAGFC = scratch_workspace+"/AGFC_1k"+Year
MeanSHRB = scratch_workspace+"/SHRB_1k"+Year
MeanTREE = scratch_workspace+"/TREE_1k"+Year
MeanLITR = scratch_workspace+"/LITR_1k"+Year
MeanHAGB = scratch_workspace+"/HAGB_1k"+Year
Range_1k = scratch_workspace+"/Range_1k"+Year
dHAGB = scratch_workspace+"/dHAGB_1k"+Year
RangeMask = scratch_workspace+"\\RangeMask"
GrassMask = scratch_workspace+"\\GrassMask"
ShrubMask = scratch_workspace+"\\ShrubMask"
study_area_ESD_Class = scratch_workspace+"/study_area_ESD_Class"
HiProdESD1k = scratch_workspace+"/HiProdESD_1k"+Year
MidProdESD1k = scratch_workspace+"/MidProdESD_1k"+Year
LowProdESD1k = scratch_workspace+"/LowProdESD_1k"+Year
NewGDB_WetlandArea = scratch_workspace+"\WetlandAreaM2_1k"
NewGDB_WetlandCount = scratch_workspace+"\WetlandCount1k"
NewGDB_Lat = scratch_workspace+"\Lat"
NewGDB_Lon = scratch_workspace+"\Lon"
NewGDB_VRM = scratch_workspace+"\VRM1k"
NewGDB_Roads_2010 = scratch_workspace+"\Roads_2010_1k"
NewGDB_Roads_2021 = scratch_workspace+"\Roads_2021_1k"
    
# Mask out cropland
    # Create cropland raster (binary, 1/0)
outraster = arcpy.sa.ReclassByTable(NASS, CropReclass,"VALUE", "VALUE", "reclass", "NODATA", ); outraster.save(CropRaster)
    
    # Create Crop/Dev/H2O mask (binary, 1/0)
outraster = arcpy.ia.SetNull(NASS, NASS, "VALUE = 0"); outraster.save(NASS_null)
out_raster = arcpy.sa.ReclassByTable(NASS_null,\
                                     DevCropH20, "VALUE", "VALUE", "reclass", "NODATA"); \
                                     out_raster.save(RapTemp)
out_raster = arcpy.ia.IsNull(RapTemp); out_raster.save(CropMask)
    
    # Set RAP values to NoData in Croplands, Development, Water
raplayer = arcpy.management.MakeRasterLayer(VegCover, "layer", '', "", 1) # AGFC is band 1
out_raster = arcpy.ia.SetNull(CropMask, \
                              raplayer, \
                              "Value = 0"); \
                              out_raster.save(AGFCCover)
raplayer = arcpy.management.MakeRasterLayer(VegCover, "layer", '', "", 2) # LITR is band 2
outraster = arcpy.ia.Con(CropMask, \
                          raplayer, \
                          None, "Value = 1");\
                          outraster.save(LITRCover)
raplayer = arcpy.management.MakeRasterLayer(VegCover, "layer", '', "", 3)# PGFC is band 3
outraster = arcpy.ia.Con(CropMask, \
                          raplayer, \
                          None, "Value = 1");\
                          outraster.save(PGFCCover)
raplayer = arcpy.management.MakeRasterLayer(VegCover, "layer", '', "", 4)# SHRB is band 4
outraster = arcpy.ia.Con(CropMask, \
                          raplayer, \
                          None, "Value = 1");\
                          outraster.save(SHRBCover)
raplayer = arcpy.management.MakeRasterLayer(VegCover, "layer", '', "", 5)# TREE is band 5
outraster = arcpy.ia.Con(CropMask, \
                          raplayer, \
                          None, "Value = 1");\
                          outraster.save(TREECover)
outraster = arcpy.ia.Con(CropMask, \
                          HAGB, \
                          None, "Value = 1");\
                          outraster.save(HAGB_masked)

# Create 1 km grid for summarizing the Covariates
    # Export study area polygon so that the coordinate system is NAD83 UTM Z13N
arcpy.FeatureClassToFeatureClass_conversion(in_features=study_area, \
                                            out_path=scratch_workspace, \
                                            out_name="Study_area_z13", \
                                            where_clause="", \
                                            field_mapping="", config_keyword="") 
    
    # Extract the Polygon Fishnet for the study area
GridSelection = arcpy.management.SelectLayerByLocation(UTMGrid, \
                                       "HAVE_THEIR_CENTER_IN", \
                                       study_area_z13, None, "NEW_SELECTION", "NOT_INVERT")
arcpy.conversion.FeatureClassToFeatureClass(GridSelection, scratch_workspace, "study_area_grid", '', '', '')


# Cropland - 9k landscape summary
    #Zonal statistics - zones = each 1 km polygon; stats = sum on the cropland raster
arcpy.sa.ZonalStatisticsAsTable(study_area_grid, \
                                "OBJECTID", \
                                CropRaster, \
                                zonaltable, "DATA", "SUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
    #Multiply this by the pixel area from the in hectares
arcpy.management.CalculateField(zonaltable, "SUM_hec", "!SUM! * 0.09", "PYTHON3", '', "DOUBLE", "NO_ENFORCE_DOMAINS")
    #Join this back to the polygon fishnet
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "SUM_hec")
    #Convert this back to a raster
arcpy.conversion.PolygonToRaster(study_area_grid, \
                                     "SUM_hec", \
                                     RapTemp, \
                                     "CELL_CENTER", "NONE", \
                                     1000, \
                                     "BUILD")
    #Take the sum in the 3x3km landscape around each cell
crops = arcpy.sa.FocalStatistics(RapTemp, \
                                      "Rectangle 3 3 CELL", "SUM", "DATA", 90); \
                                      crops.save(Crop9k)

# 1km scale masked cover rasters
    #Zonal statistics - zones = each 1 km polygon; stats = mean on the PGFC raster
arcpy.sa.ZonalStatisticsAsTable(study_area_grid, \
                                "OBJECTID", \
                                PGFCCover, \
                                zonaltable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "MEAN")
    #Convert this back to a raster
arcpy.conversion.PolygonToRaster(study_area_grid, \
                                     "MEAN", \
                                     MeanPGFC, \
                                     "CELL_CENTER", "NONE", \
                                     1000, \
                                     "BUILD")
arcpy.management.DeleteField(study_area_grid, "MEAN", "DELETE_FIELDS")
    #Zonal statistics - zones = each 1 km polygon; stats = mean on the AGFC raster
arcpy.sa.ZonalStatisticsAsTable(study_area_grid, \
                                "OBJECTID", \
                                AGFCCover, \
                                zonaltable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "MEAN")
    #Convert this back to a raster
arcpy.conversion.PolygonToRaster(study_area_grid, \
                                     "MEAN", \
                                     MeanAGFC, \
                                     "CELL_CENTER", "NONE", \
                                     1000, \
                                     "BUILD")
arcpy.management.DeleteField(study_area_grid, "MEAN", "DELETE_FIELDS")
    #Zonal statistics - zones = each 1 km polygon; stats = mean on the SHRB raster
arcpy.sa.ZonalStatisticsAsTable(study_area_grid, \
                                "OBJECTID", \
                                SHRBCover, \
                                zonaltable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "MEAN")
    #Convert this back to a raster
arcpy.conversion.PolygonToRaster(study_area_grid, \
                                     "MEAN", \
                                     MeanSHRB, \
                                     "CELL_CENTER", "NONE", \
                                     1000, \
                                     "BUILD")
arcpy.management.DeleteField(study_area_grid, "MEAN", "DELETE_FIELDS")
    #Zonal statistics - zones = each 1 km polygon; stats = mean on the TREE raster
arcpy.sa.ZonalStatisticsAsTable(study_area_grid, \
                                "OBJECTID", \
                                TREECover, \
                                zonaltable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "MEAN")
    #Convert this back to a raster
arcpy.conversion.PolygonToRaster(study_area_grid, \
                                     "MEAN", \
                                     MeanTREE, \
                                     "CELL_CENTER", "NONE", \
                                     1000, \
                                     "BUILD")
arcpy.management.DeleteField(study_area_grid, "MEAN", "DELETE_FIELDS")
    #Zonal statistics - zones = each 1 km polygon; stats = mean on the LITR raster
arcpy.sa.ZonalStatisticsAsTable(study_area_grid, \
                                "OBJECTID", \
                                LITRCover, \
                                zonaltable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "MEAN")
    #Convert this back to a raster
arcpy.conversion.PolygonToRaster(study_area_grid, \
                                     "MEAN", \
                                     MeanLITR, \
                                     "CELL_CENTER", "NONE", \
                                     1000, \
                                     "BUILD")
arcpy.management.DeleteField(study_area_grid, "MEAN", "DELETE_FIELDS")
    
    #Zonal statistics - zones = each 1 km polygon; stats = mean on the HAGB raster
arcpy.sa.ZonalStatisticsAsTable(study_area_grid, \
                                "OBJECTID", \
                                HAGB_masked, \
                                zonaltable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "MEAN")
    #Convert this back to a raster
arcpy.conversion.PolygonToRaster(study_area_grid, \
                                     "MEAN", \
                                     MeanHAGB, \
                                     "CELL_CENTER", "NONE", \
                                     1000, \
                                     "BUILD")
arcpy.management.DeleteField(study_area_grid, "MEAN", "DELETE_FIELDS")
    # Create the dHAGB for the study area
out_raster = Raster(MeanHAGB) - Raster(HAGBMean1k);\
             out_raster.save(dHAGB)
    
    
    
    
# Create the ESD Class Rasters - Set it Null in Crops, Dev, Water, Etc
    # Base ESD class raster to extract, mask out crops
outraster = arcpy.ia.Con(CropMask, \
                          ESD_prod_class, \
                          None, "Value = 1");\
                          outraster.save(study_area_ESD_Class)
    # Mask out trees >4
out_raster = arcpy.ia.SetNull(TREECover,\
                              study_area_ESD_Class, "Value >= 5"); \
                              out_raster.save(RapTemp)
    # Mask out areas with PGFC < 30 & SHRB < 5
        # Create the mask
out_raster = arcpy.sa.Reclassify(PGFCCover, \
                                 "Value", \
                                 "0 30 0;31 101 1", "DATA");\
                                 out_raster.save(GrassMask)
out_raster = arcpy.sa.Reclassify(SHRBCover, \
                                 "Value", \
                                 "0 4 0;5 101 1", "DATA");\
                                 out_raster.save(ShrubMask)
out_raster = Raster(ShrubMask) + Raster(GrassMask);
out_raster.save(RangeMask)
        # Create a rangeland percentage raster
out_raster = arcpy.sa.SetNull(RangeMask, 1, "Value = 0"); \
             out_raster.save(Rangeland)
arcpy.ia.ZonalStatisticsAsTable(study_area_grid,\
                                "OBJECTID", Rangeland, \
                                zonaltable,\
                                "DATA", "SUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "AREA")
with arcpy.da.UpdateCursor(study_area_grid, ["AREA"]) as cursor:
    for row in cursor:
        if row[0] == None:
            row[0] = 0
            cursor.updateRow(row)
arcpy.management.CalculateField(study_area_grid, "Range_pct", "!AREA! / 10000", "PYTHON3", '', "DOUBLE", "NO_ENFORCE_DOMAINS")
arcpy.conversion.PolygonToRaster(study_area_grid, "Range_pct", Range_1k, "CELL_CENTER", "NONE", 1000, "BUILD")
        # Do the mask
out_raster = arcpy.ia.SetNull(RangeMask,\
                              RapTemp, "Value = 0"); \
                              out_raster.save(study_area_ESD_Class)
        
        # Split it into the three class rasters, and create the zonal stats raster at the 1k level
        
            # Low Productivity ESDs
            
out_raster = arcpy.sa.Reclassify(study_area_ESD_Class, \
                                 "Value", \
                                 "0 1 1", "NODATA");\
                                 out_raster.save(RapTemp)
arcpy.sa.ZonalStatisticsAsTable(study_area_grid, \
                                "OBJECTID", \
                                RapTemp, \
                                zonaltable, \
                                "DATA", "SUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "AREA")
arcpy.management.AlterField(study_area_grid, "AREA", "Area_Low_ESD_m2", "Area_Low_ESD_m2", "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")
            # Mid Productivity ESDs
            
out_raster = arcpy.sa.Reclassify(study_area_ESD_Class, \
                                 "Value", \
                                 "2 2 1", "NODATA");\
                                 out_raster.save(RapTemp)
arcpy.sa.ZonalStatisticsAsTable(study_area_grid, \
                                "OBJECTID", \
                                RapTemp, \
                                zonaltable, \
                                "DATA", "SUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "AREA")
arcpy.management.AlterField(study_area_grid, "AREA", "Area_Mid_ESD_m2", "Area_Mid_ESD_m2", "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")
            # High Productivity ESDs
            
out_raster = arcpy.sa.Reclassify(study_area_ESD_Class, \
                                 "Value", \
                                 "3 3 1", "NODATA");\
                                 out_raster.save(RapTemp)
arcpy.sa.ZonalStatisticsAsTable(study_area_grid, \
                                "OBJECTID", \
                                RapTemp, \
                                zonaltable, \
                                "DATA", "SUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
arcpy.management.JoinField(study_area_grid, "OBJECTID", zonaltable, "OBJECTID_1", "AREA")
arcpy.management.AlterField(study_area_grid, "AREA", "Area_High_ESD_m2", "Area_High_ESD_m2", "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")
    
            # Set the Null Values to Zero, values >100 to 100, convert to pct
with arcpy.da.UpdateCursor(study_area_grid, ["Area_High_ESD_m2", "Area_Mid_ESD_m2","Area_Low_ESD_m2"]) as curU:  
    for row in curU:  
        rowU = row  
        for field in range(3):  
            if rowU[field] == None:  
                rowU[field] = "0"  
            if rowU[field] != None:
                rowU[field] = int(rowU[field])/10000
            if int(rowU[field]) >= 1000000:
                rowU[field] = 100
        curU.updateRow(rowU)
del curU
            # Create the Rasters
arcpy.conversion.PolygonToRaster(study_area_grid, \
                                     "Area_High_ESD_m2", \
                                     HiProdESD1k, \
                                     "CELL_CENTER", "NONE", \
                                     1000, \
                                     "BUILD")
arcpy.conversion.PolygonToRaster(study_area_grid, \
                                     "Area_Mid_ESD_m2", \
                                     MidProdESD1k, \
                                     "CELL_CENTER", "NONE", \
                                     1000, \
                                     "BUILD")
arcpy.conversion.PolygonToRaster(study_area_grid, \
                                     "Area_Low_ESD_m2", \
                                     LowProdESD1k, \
                                     "CELL_CENTER", "NONE", \
                                     1000, \
                                     "BUILD")
    
# Clip all the static rasters
    
arcpy.management.Clip(Roads2010,\
                      "", NewGDB_Roads_2010,\
                      study_area_grid, "1.79e+308", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
arcpy.management.Clip(Roads2021,\
                      "", NewGDB_Roads_2021,\
                      study_area_grid, "1.79e+308", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
arcpy.management.Clip(Lat,\
                      "", NewGDB_Lat,\
                      study_area_grid, "1.79e+308", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
arcpy.management.Clip(Lon,\
                      "", NewGDB_Lon,\
                      study_area_grid, "1.79e+308", "ClippingGeometry", "NO_MAINTAIN_EXTENT")    
arcpy.management.Clip(VRM,\
                      "", NewGDB_VRM,\
                      study_area_grid, "1.79e+308", "ClippingGeometry", "NO_MAINTAIN_EXTENT")      
arcpy.management.Clip(WetlandArea,\
                      "", NewGDB_WetlandArea,\
                      study_area_grid, "1.79e+308", "ClippingGeometry", "NO_MAINTAIN_EXTENT")    
arcpy.management.Clip(WetlandCount,\
                      "", NewGDB_WetlandCount,\
                      study_area_grid, "1.79e+308", "ClippingGeometry", "NO_MAINTAIN_EXTENT")


# Delete the intermediate layers to de-clutter and save space
arcpy.management.Delete(scratch_workspace+"\AGFC;"\
                        +scratch_workspace+"\Cropland;"\
                        +scratch_workspace+"\CropMask;"\
                        +scratch_workspace+"\GrassMask;"\
                        +scratch_workspace+"\HAGB_1k;"\
                        +scratch_workspace+"\HAGB_masked;"\
                        +scratch_workspace+"\LITR_1k;"\
                        +scratch_workspace+"\LITR;"\
                        +scratch_workspace+"\PGFC;"\
                        +scratch_workspace+"\Rangeland;"\
                        +scratch_workspace+"\RangeMask;"\
                        +scratch_workspace+"\RapTemp;"\
                        +scratch_workspace+"\SetNull_NASS;"\
                        +scratch_workspace+"\SHRB;"\
                        +scratch_workspace+"\ShrubMask;"\
                        +scratch_workspace+"\TREE;"\
                        +scratch_workspace+"\zonaltable", '')
