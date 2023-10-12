### Modules, Environment Settings, Etc ###
import time
import arcpy
import numpy as np
import math
import pandas as pd
import os
import matplotlib.pyplot as plt
from arcpy import env
from arcpy.sa import *
from arcpy import mp
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("spatial")


#Parameters:

    # Species
Spp = arcpy.GetParameterAsText(0)
#Spp = 'BAIS' # For testing

    # ModelLayersGDB
ModelLayersGDB = arcpy.GetParameterAsText(1)
#ModelLayersGDB = r'C:\Users\eric.chabot\Documents\NGPJV_spatial_tool\Tool_Input_Rasters\ModelLayers.gdb' # For testing
arcpy.env.workspace = ModelLayersGDB

ResultsGDB = arcpy.GetParameterAsText(2)
if ResultsGDB == '':
    ResultsGDB_path = arcpy.Describe(ModelLayersGDB).path
    ResultsGDB_name = "Results_"+Spp+".gdb"
    arcpy.management.CreateFileGDB(ResultsGDB_path, ResultsGDB_name, "CURRENT")
    ResultsGDB = ResultsGDB_path +"\\" + ResultsGDB_name
    
    # Area of Interest
AOI_fc = arcpy.GetParameterAsText(3)
#AOI_fc = r'C:\Users\eric.chabot\Documents\NGPJV_spatial_tool\Tool_Input_Rasters\ModelLayers.gdb\NGP_studyarea_USA' # For testing
#AOI_fc = r'C:\Users\eric.chabot\Documents\NGPJV_spatial_tool\Test_Data_acquisition_tool\TestRun08042022\drive-download-20220804T194643Z-001\StudyArea.shp'

AOI = "StudyArea"
arcpy.management.MakeFeatureLayer(AOI_fc,\
                                      AOI, \
                                      '', None,\
                                      "")

############
    # Management levels to measure as an effect on future abundance
    
#AGFCpct = 0 # For testing
AGFCRas = arcpy.GetParameterAsText(4)
#PGFCpct = 0 # For testing
PGFCRas = arcpy.GetParameterAsText(5)
#SHRBpct = 200 # For testing
SHRBRas = arcpy.GetParameterAsText(6)
#TREEpct = 5000 # For testing
TREERas = arcpy.GetParameterAsText(7)

GZRas = arcpy.GetParameterAsText(8)
#EApct = 0 # For testing
EARas = arcpy.GetParameterAsText(9)

    # Number of years to predict beyond the latest data
#PredYrs = 5 # For testing
PredYrs = arcpy.GetParameterAsText(10)
PredictionYears = int(PredYrs)


    # Confidence Interval?
CIsize = arcpy.GetParameterAsText(11)
#CIsize = 95
per1 = int(CIsize)
per2 = 100 - per1


# Bring in the file containing all the parameter estimates:

#coef_dataframe = table_to_data_frame(ModelLayersGDB+'\\'+Spp+'_pars') # For testing
coef_dataframe = pd.read_csv(arcpy.Describe(ModelLayersGDB).path+"\\Pars_"+Spp+".csv")
#coef_dataframe = pd.read_csv(r'C:\Users\eric.chabot\Documents\NGPJV_spatial_tool\Create_Prediction_Map_tool\Test_prediction'+"\\Pars_"+Spp+".csv") # For testing


# Define the expit function used in prediction

def expit(x):
    return np.exp(x)/(1+np.exp(x))


# Create a dictionary of mean and SD values for each covariate, used to scale and center them

Scale_dataframe = pd.read_csv(arcpy.Describe(ModelLayersGDB).path+"\\Pars_scale.csv")
rightside = zip(Scale_dataframe.mn, Scale_dataframe.sd)
full_dict = zip(Scale_dataframe.Cov, list(rightside))
ZDict = dict(full_dict)

# Create a dictionary for each species' area correction factor
Spp_correction_dataframe = pd.read_csv(arcpy.Describe(ModelLayersGDB).path+"\\spp_area_correction.csv")
Spp_dict = dict(zip(Spp_correction_dataframe.Spp_code, list(Spp_correction_dataframe.correction_factor)))
Spp_Area_Correction = Spp_dict[Spp] # Set the correction factor


# Create a temporary model layers workspace where clipped rasters can be stored, and clip them to the study area
arcpy.management.CreateFileGDB(arcpy.Describe(ModelLayersGDB).path, "TempModelLayers", "CURRENT")
TempModelLayersGDB = arcpy.Describe(ModelLayersGDB).path+'\TempModelLayers.gdb'
with arcpy.EnvManager(scratchWorkspace=ModelLayersGDB,\
                      workspace=ModelLayersGDB):
    CovariateRasters = arcpy.ListRasters()

arcpy.AddMessage("Clipping Rasters to Study Area...")
ras_count = len(CovariateRasters)
arcpy.SetProgressor("step", "Clipping Rasters...",
                    0, ras_count, 1)
for ras in CovariateRasters:
    out_raster = os.path.join(TempModelLayersGDB, os.path.basename(ras))
    in_raster = os.path.join(ModelLayersGDB, os.path.basename(ras))
    arcpy.management.Clip(in_raster,\
                          "", \
                          out_raster,\
                          AOI,\
                          "1.79e+308", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    arcpy.SetProgressorPosition()
    
arcpy.ResetProgressor()
ModelLayersGDB = TempModelLayersGDB

# Calculate median rasters for projection beyond when data is available
out_raster = arcpy.sa.CellStatistics(\
    ModelLayersGDB + '\AGFC_1k2010;'+\
    ModelLayersGDB + '\AGFC_1k2011;'+\
    ModelLayersGDB + '\AGFC_1k2012;'+\
    ModelLayersGDB + '\AGFC_1k2013;'+\
    ModelLayersGDB + '\AGFC_1k2014;'+\
    ModelLayersGDB + '\AGFC_1k2015;'+\
    ModelLayersGDB + '\AGFC_1k2016;'+\
    ModelLayersGDB + '\AGFC_1k2017;'+\
    ModelLayersGDB + '\AGFC_1k2018;'+\
    ModelLayersGDB + '\AGFC_1k2019;'+\
    ModelLayersGDB + '\AGFC_1k2020;'+\
    ModelLayersGDB + '\AGFC_1k2021',\
    "MEDIAN", "DATA", "SINGLE_BAND", 90, "AUTO_DETECT"); \
    out_raster.save(ModelLayersGDB + '\AGFC_median')

out_raster = arcpy.sa.CellStatistics(\
    ModelLayersGDB + '\PGFC_1k2010;'+\
    ModelLayersGDB + '\PGFC_1k2011;'+\
    ModelLayersGDB + '\PGFC_1k2012;'+\
    ModelLayersGDB + '\PGFC_1k2013;'+\
    ModelLayersGDB + '\PGFC_1k2014;'+\
    ModelLayersGDB + '\PGFC_1k2015;'+\
    ModelLayersGDB + '\PGFC_1k2016;'+\
    ModelLayersGDB + '\PGFC_1k2017;'+\
    ModelLayersGDB + '\PGFC_1k2018;'+\
    ModelLayersGDB + '\PGFC_1k2019;'+\
    ModelLayersGDB + '\PGFC_1k2020;'+\
    ModelLayersGDB + '\PGFC_1k2021',\
    "MEDIAN", "DATA", "SINGLE_BAND", 90, "AUTO_DETECT"); \
    out_raster.save(ModelLayersGDB + '\PGFC_median')

out_raster = arcpy.sa.CellStatistics(\
    ModelLayersGDB + '\TREE_1k2010;'+\
    ModelLayersGDB + '\TREE_1k2011;'+\
    ModelLayersGDB + '\TREE_1k2012;'+\
    ModelLayersGDB + '\TREE_1k2013;'+\
    ModelLayersGDB + '\TREE_1k2014;'+\
    ModelLayersGDB + '\TREE_1k2015;'+\
    ModelLayersGDB + '\TREE_1k2016;'+\
    ModelLayersGDB + '\TREE_1k2017;'+\
    ModelLayersGDB + '\TREE_1k2018;'+\
    ModelLayersGDB + '\TREE_1k2019;'+\
    ModelLayersGDB + '\TREE_1k2020;'+\
    ModelLayersGDB + '\TREE_1k2021',\
    "MEDIAN", "DATA", "SINGLE_BAND", 90, "AUTO_DETECT"); \
    out_raster.save(ModelLayersGDB + '\TREE_median')


out_raster = arcpy.sa.CellStatistics(\
    ModelLayersGDB + '\SHRB_1k2010;'+\
    ModelLayersGDB + '\SHRB_1k2011;'+\
    ModelLayersGDB + '\SHRB_1k2012;'+\
    ModelLayersGDB + '\SHRB_1k2013;'+\
    ModelLayersGDB + '\SHRB_1k2014;'+\
    ModelLayersGDB + '\SHRB_1k2015;'+\
    ModelLayersGDB + '\SHRB_1k2016;'+\
    ModelLayersGDB + '\SHRB_1k2017;'+\
    ModelLayersGDB + '\SHRB_1k2018;'+\
    ModelLayersGDB + '\SHRB_1k2019;'+\
    ModelLayersGDB + '\SHRB_1k2020;'+\
    ModelLayersGDB + '\SHRB_1k2021',\
    "MEDIAN", "DATA", "SINGLE_BAND", 90, "AUTO_DETECT"); \
    out_raster.save(ModelLayersGDB + '\SHRB_median')

out_raster = arcpy.sa.CellStatistics(\
    ModelLayersGDB + '\dhagb_1k2010;'+\
    ModelLayersGDB + '\dhagb_1k2011;'+\
    ModelLayersGDB + '\dhagb_1k2012;'+\
    ModelLayersGDB + '\dhagb_1k2013;'+\
    ModelLayersGDB + '\dhagb_1k2014;'+\
    ModelLayersGDB + '\dhagb_1k2015;'+\
    ModelLayersGDB + '\dhagb_1k2016;'+\
    ModelLayersGDB + '\dhagb_1k2017;'+\
    ModelLayersGDB + '\dhagb_1k2018;'+\
    ModelLayersGDB + '\dhagb_1k2019;'+\
    ModelLayersGDB + '\dhagb_1k2020;'+\
    ModelLayersGDB + '\dhagb_1k2021',\
    "MEDIAN", "DATA", "SINGLE_BAND", 90, "AUTO_DETECT"); \
    out_raster.save(ModelLayersGDB + '\dhagb_median')

out_raster = arcpy.sa.CellStatistics(\
    ModelLayersGDB + '\HiProdESD_1k2010;'+\
    ModelLayersGDB + '\HiProdESD_1k2011;'+\
    ModelLayersGDB + '\HiProdESD_1k2012;'+\
    ModelLayersGDB + '\HiProdESD_1k2013;'+\
    ModelLayersGDB + '\HiProdESD_1k2014;'+\
    ModelLayersGDB + '\HiProdESD_1k2015;'+\
    ModelLayersGDB + '\HiProdESD_1k2016;'+\
    ModelLayersGDB + '\HiProdESD_1k2017;'+\
    ModelLayersGDB + '\HiProdESD_1k2018;'+\
    ModelLayersGDB + '\HiProdESD_1k2019;'+\
    ModelLayersGDB + '\HiProdESD_1k2020;'+\
    ModelLayersGDB + '\HiProdESD_1k2021',\
    "MEDIAN", "DATA", "SINGLE_BAND", 90, "AUTO_DETECT"); \
    out_raster.save(ModelLayersGDB + '\HiProdESD_median')

out_raster = arcpy.sa.CellStatistics(\
    ModelLayersGDB + '\MidProdESD_1k2010;'+\
    ModelLayersGDB + '\MidProdESD_1k2011;'+\
    ModelLayersGDB + '\MidProdESD_1k2012;'+\
    ModelLayersGDB + '\MidProdESD_1k2013;'+\
    ModelLayersGDB + '\MidProdESD_1k2014;'+\
    ModelLayersGDB + '\MidProdESD_1k2015;'+\
    ModelLayersGDB + '\MidProdESD_1k2016;'+\
    ModelLayersGDB + '\MidProdESD_1k2017;'+\
    ModelLayersGDB + '\MidProdESD_1k2018;'+\
    ModelLayersGDB + '\MidProdESD_1k2019;'+\
    ModelLayersGDB + '\MidProdESD_1k2020;'+\
    ModelLayersGDB + '\MidProdESD_1k2021',\
    "MEDIAN", "DATA", "SINGLE_BAND", 90, "AUTO_DETECT"); \
    out_raster.save(ModelLayersGDB + '\MidProdESD_median')





# Figure out when the latest year of data is, in case additional years of raster data are added using the update tool

templist = []
for i in CovariateRasters:
    templist.extend([i[-4:]])
templist2 = numpy.unique(templist).tolist()
result = [val for val in templist2 if val.isdigit()]
yearlist = [eval(i) for i in result]
yearlist.sort()


# Create Future-year rasters to predict management effects

    # Create a year list including future years set by the prediction time:
future_years = [max(yearlist)+1]
for _ in range(PredictionYears-1):
    future_years.append(max(future_years)+1)
all_years = yearlist + future_years
RefYear = int(max(yearlist))
    

arcpy.AddMessage("Generating Management Rasters...")
arcpy.SetProgressor("step", "Generating Management Rasters...",
                    0, 100, 12)
if TREERas != '':
    out_ras = Raster(TREERas)
    out_ras.save(ModelLayersGDB + '\\TREE_1k_mgmt')
    arcpy.SetProgressorPosition()
else:
    out_ras = ModelLayersGDB + '\TREE_median'
    out_ras.save(ModelLayersGDB + '\\TREE_1k_mgmt')
    arcpy.SetProgressorPosition()


if AGFCRas != '':
    out_ras = Raster(AGFCRas)
    out_ras.save(ModelLayersGDB + '\\AGFC_1k_mgmt')
    arcpy.SetProgressorPosition()
else:
    out_ras = ModelLayersGDB + '\AGFC_median'
    out_ras.save(ModelLayersGDB + '\\AGFC_1k_mgmt')
    arcpy.SetProgressorPosition()
    
if SHRBRas != '':
    out_ras = Raster(SHRBRas)
    out_ras.save(ModelLayersGDB + '\\SHRB_1k_mgmt')
    arcpy.SetProgressorPosition()
else:
    out_ras = ModelLayersGDB + '\SHRB_median'
    out_ras.save(ModelLayersGDB + '\\SHRB_1k_mgmt')
    arcpy.SetProgressorPosition()
    
if PGFCRas != '':
    out_ras = Raster(PGFCRas)
    out_ras.save(ModelLayersGDB + '\\PGFC_1k_mgmt')
    arcpy.SetProgressorPosition()
else:
    out_ras = ModelLayersGDB + '\PGFC_median'
    out_ras.save(ModelLayersGDB + '\\PGFC_1k_mgmt')
    arcpy.SetProgressorPosition()


if GZRas != '':
    out_ras = Raster(GZRas)
    out_ras.save(ModelLayersGDB + '\\Trt_GZ_mgmt')
    arcpy.SetProgressorPosition()
else:
    out_ras = ModelLayersGDB + '\Trt_GZ'
    out_ras.save(ModelLayersGDB + '\\Trt_GZ_mgmt')
    arcpy.SetProgressorPosition()

if EARas != '':
    out_ras = Raster(EARas)
    out_ras.save(ModelLayersGDB + '\\Trt_EA_mgmt')
    arcpy.SetProgressorPosition()
else:
    out_ras = ModelLayersGDB + '\Trt_EA'
    out_ras.save(ModelLayersGDB + '\\Trt_EA_mgmt')
    arcpy.SetProgressorPosition()


#             arcpy.SetProgressorPosition()

arcpy.ResetProgressor()

arcpy.AddMessage("Scaling and Centering Covariates...")
arcpy.SetProgressor("step", "Scaling and Centering Covariates...",
                    0, len(all_years), 1)
# y is year, scale and center covariate rasters and convert these to numpy arrays for faster processing
    # Static covariates:
Trt_EA = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\Trt_EA').astype('float32')-ZDict["Trt_EA"][0])/ZDict["Trt_EA"][1]

Trt_EA_mgmt = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\Trt_EA_mgmt')-ZDict["Trt_EA"][0])/ZDict["Trt_EA"][1]

Trt_GZ = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\Trt_GZ').astype('float32')-ZDict["Trt_GZ"][0])/ZDict["Trt_GZ"][1]

Trt_GZ_mgmt = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\Trt_GZ_mgmt')-ZDict["Trt_GZ"][0])/ZDict["Trt_GZ"][1]

Wetland_Count = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\WetlandCount1k').astype('float32')-ZDict["Wetland_Count"][0])/ZDict["Wetland_Count"][1]
Wetland_Area = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\WetlandAreaM2_1k').astype('float32')-ZDict["Wetland_Area"][0])/ZDict["Wetland_Area"][1]
VRM = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\VRM1k').astype('float32')-ZDict["VRM"][0])/ZDict["VRM"][1]
Road_2010 = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\Roads_2010_1k').astype('float32')-ZDict["Road_km"][0])/ZDict["Road_km"][1]
Road_2021 = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\Roads_2021_1k').astype('float32')-ZDict["Road_km"][0])/ZDict["Road_km"][1]
Lat = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\Lat').astype('float32')-ZDict["Lat"][0])/ZDict["Lat"][1]
Lat2 = Lat*Lat
Lon = (arcpy.RasterToNumPyArray(ModelLayersGDB + '\Lon').astype('float32')-ZDict["Lon"][0])/ZDict["Lon"][1]
Lon2 = Lon*Lon    
    # Covariates that change each year
        # Each year's covariate array will be stored in a list for that covariate
            #first initialize the empty lists:
Crop_pct = []
Range_pct = []
agfc = []
pgfc = []
shrb = []
tree = []
agfc_mgmt = []
pgfc_mgmt = []
shrb_mgmt = []
tree_mgmt = []
devprd = []
ESDmid = []
ESDhi =[]
pgfcXdevprd = []
pgfcXdevprd_mgmt = []
ESD_mid_prodXdevprd = []
ESD_high_prodXdevprd = []
            #Then fill em up
for y in yearlist:
    Crop_pct.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\Crop9k_1k'+str(y)).astype('float32'))-ZDict["Crop_pct"][0])/ZDict["Crop_pct"][1])
    Range_pct.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\Range_1k'+str(y)).astype('float32'))-ZDict["RANG_pct"][0])/ZDict["RANG_pct"][1])
    agfc.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\AGFC_1k'+str(y)).astype('float32'))-ZDict["agfc"][0])/ZDict["agfc"][1])
    pgfc.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\PGFC_1k'+str(y)).astype('float32'))-ZDict["pgfc"][0])/ZDict["pgfc"][1])
    shrb.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\SHRB_1k'+str(y)).astype('float32'))-ZDict["shrb"][0])/ZDict["shrb"][1])
    tree.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\TREE_1k'+str(y)).astype('float32'))-ZDict["tree"][0])/ZDict["tree"][1])
    agfc_mgmt.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\AGFC_1k'+str(y)).astype('float32'))-ZDict["agfc"][0])/ZDict["agfc"][1])
    pgfc_mgmt.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\PGFC_1k'+str(y)).astype('float32'))-ZDict["pgfc"][0])/ZDict["pgfc"][1])
    shrb_mgmt.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\SHRB_1k'+str(y)).astype('float32'))-ZDict["shrb"][0])/ZDict["shrb"][1])
    tree_mgmt.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\TREE_1k'+str(y)).astype('float32'))-ZDict["tree"][0])/ZDict["tree"][1])    
    devprd.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\dhagb_1k'+str(y)).astype('float32'))-ZDict["devprd"][0])/ZDict["devprd"][1])
    ESDmid.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\MidProdESD_1k'+str(y)).astype('float32'))-ZDict["ESD_mid_prod"][0])/ZDict["ESD_mid_prod"][1])
    ESDhi.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\HiProdESD_1k'+str(y)).astype('float32'))-ZDict["ESD_high_prod"][0])/ZDict["ESD_high_prod"][1])

    pgfcXdevprd.append(\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\dhagb_1k'+str(y)).astype('float32')-ZDict["devprd"][0])/ZDict["devprd"][1])*\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\PGFC_1k'+str(y)).astype('float32')-ZDict["pgfc"][0])/ZDict["pgfc"][1])\
        )
    
    pgfcXdevprd_mgmt.append(\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\dhagb_1k'+str(y)).astype('float32')-ZDict["devprd"][0])/ZDict["devprd"][1])*\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\PGFC_1k'+str(y)).astype('float32')-ZDict["pgfc"][0])/ZDict["pgfc"][1])\
        )


    ESD_mid_prodXdevprd.append(\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\MidProdESD_1k'+str(y)).astype('float32')-ZDict["ESD_mid_prod"][0])/ZDict["ESD_mid_prod"][1])*\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\dhagb_1k'+str(y)).astype('float32')-ZDict["devprd"][0])/ZDict["devprd"][1])\
        )
    
    ESD_high_prodXdevprd.append(\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\HiProdESD_1k'+str(y)).astype('float32')-ZDict["ESD_high_prod"][0])/ZDict["ESD_high_prod"][1])*\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\dhagb_1k'+str(y)).astype('float32')-ZDict["devprd"][0])/ZDict["devprd"][1])\
        )
    arcpy.SetProgressorPosition()
# Extend the lists using the management rasters, or the existing median rasters for those without estimated management effects
for y in future_years:
    
    # No estimated management effect:
    Crop_pct.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\Crop9k_1k'+str(max(yearlist))).astype('float32'))-ZDict["Crop_pct"][0])/ZDict["Crop_pct"][1])
    Range_pct.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\Range_1k'+str(max(yearlist))).astype('float32'))-ZDict["RANG_pct"][0])/ZDict["RANG_pct"][1])
    devprd.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\dhagb_median').astype('float32'))-ZDict["devprd"][0])/ZDict["devprd"][1])
    ESDmid.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\MidProdESD_median').astype('float32'))-ZDict["ESD_mid_prod"][0])/ZDict["ESD_mid_prod"][1])
    ESDhi.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\HiProdESD_median').astype('float32'))-ZDict["ESD_high_prod"][0])/ZDict["ESD_high_prod"][1])


    
    ESD_mid_prodXdevprd.append(\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\MidProdESD_median').astype('float32')-ZDict["ESD_mid_prod"][0])/ZDict["ESD_mid_prod"][1])*\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\dhagb_median').astype('float32')-ZDict["devprd"][0])/ZDict["devprd"][1])\
        )
    
    ESD_high_prodXdevprd.append(\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\HiProdESD_median').astype('float32')-ZDict["ESD_high_prod"][0])/ZDict["ESD_high_prod"][1])*\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\dhagb_median').astype('float32')-ZDict["devprd"][0])/ZDict["devprd"][1])\
        )
    
    # Estimated management effect on these covariates
    agfc_mgmt.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\AGFC_1k_mgmt').astype('float32'))-ZDict["agfc"][0])/ZDict["agfc"][1])
    pgfc_mgmt.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\PGFC_1k_mgmt').astype('float32'))-ZDict["pgfc"][0])/ZDict["pgfc"][1])
    shrb_mgmt.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\SHRB_1k_mgmt').astype('float32'))-ZDict["shrb"][0])/ZDict["shrb"][1])
    tree_mgmt.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\TREE_1k_mgmt').astype('float32'))-ZDict["tree"][0])/ZDict["tree"][1])
    pgfcXdevprd_mgmt.append(\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\dhagb_median').astype('float32')-ZDict["devprd"][0])/ZDict["devprd"][1])*\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\PGFC_1k_mgmt').astype('float32')-ZDict["pgfc"][0])/ZDict["pgfc"][1])\
        )
    
    # Baseline values for covariates with a management effect:
    agfc.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\AGFC_median').astype('float32'))-ZDict["agfc"][0])/ZDict["agfc"][1])
    pgfc.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\PGFC_median').astype('float32'))-ZDict["pgfc"][0])/ZDict["pgfc"][1])
    shrb.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\SHRB_median').astype('float32'))-ZDict["shrb"][0])/ZDict["shrb"][1])
    tree.append(((arcpy.RasterToNumPyArray(ModelLayersGDB + '\TREE_median').astype('float32'))-ZDict["tree"][0])/ZDict["tree"][1])
    pgfcXdevprd.append(\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\dhagb_median').astype('float32')-ZDict["devprd"][0])/ZDict["devprd"][1])*\
        ((arcpy.RasterToNumPyArray(ModelLayersGDB + '\PGFC_median').astype('float32')-ZDict["pgfc"][0])/ZDict["pgfc"][1])\
        )
    arcpy.SetProgressorPosition()
arcpy.ResetProgressor()
start = time.time()
    
mapping2010 = {"Trt_EA":Trt_EA,
           "Trt_GZ":Trt_GZ,
           "Lat":Lat,
           "Lon":Lon,
           "Lat2":Lat2,
           "Lon2":Lon2,
           "Wetland_Count":Wetland_Count,
           "Wetland_Area":Wetland_Area,
           "VRM":VRM,
           "Road_km":Road_2010,
           "agfc":agfc[0],
           "pgfc":pgfc[0],
           "shrb":shrb[0],
           "tree":tree[0],
           "ESD_mid_prod":ESDmid[0],
           "ESD_high_prod":ESDhi[0],
           "Crop_pct":Crop_pct[0],
           "RANG_pct":Range_pct[0],
           "devprd":devprd[0],
           "pgfcXdevprd":pgfcXdevprd[0],
           "ESD_mid_prodXdevprd":ESD_mid_prodXdevprd[0],
           "ESD_high_prodXdevprd":ESD_high_prodXdevprd[0]
           }

def create_mapping(y):
    global mapping
    mapping = {"Trt_EA":Trt_EA,
           "Trt_GZ":Trt_GZ,
           "Lat":Lat,
           "Lon":Lon,
           "Lat2":Lat2,
           "Lon2":Lon2,
           "Wetland_Count":Wetland_Count,
           "Wetland_Area":Wetland_Area,
           "VRM":VRM,
           "Road_km":Road_2021,
           "agfc":agfc[y-2010],
           "pgfc":pgfc[y-2010],
           "shrb":shrb[y-2010],
           "tree":tree[y-2010],
           "ESD_mid_prod":ESDmid[y-2010],
           "ESD_high_prod":ESDhi[y-2010],
           "Crop_pct":Crop_pct[y-2010],
           "RANG_pct":Range_pct[y-2010],
           "devprd":devprd[y-2010],
           "pgfcXdevprd":pgfcXdevprd[y-2010],
           "ESD_mid_prodXdevprd":ESD_mid_prodXdevprd[y-2010],
           "ESD_high_prodXdevprd":ESD_high_prodXdevprd[y-2010]
           }
    
def create_mgmt_mapping(y):
    global mapping
    mapping = {"Trt_EA":Trt_EA_mgmt,
           "Trt_GZ":Trt_GZ_mgmt,
           "Lat":Lat,
           "Lon":Lon,
           "Lat2":Lat2,
           "Lon2":Lon2,
           "Wetland_Count":Wetland_Count,
           "Wetland_Area":Wetland_Area,
           "VRM":VRM,
           "Road_km":Road_2021,
           "agfc":agfc_mgmt[y-2010],
           "pgfc":pgfc_mgmt[y-2010],
           "shrb":shrb_mgmt[y-2010],
           "tree":tree_mgmt[y-2010],
           "ESD_mid_prod":ESDmid[y-2010],
           "ESD_high_prod":ESDhi[y-2010],
           "Crop_pct":Crop_pct[y-2010],
           "RANG_pct":Range_pct[y-2010],
           "devprd":devprd[y-2010],
           "pgfcXdevprd":pgfcXdevprd_mgmt[y-2010],
           "ESD_mid_prodXdevprd":ESD_mid_prodXdevprd[y-2010],
           "ESD_high_prodXdevprd":ESD_high_prodXdevprd[y-2010]
           }
    

# Define the function to calculate annual abundance and trend

# x = model iteration (row), y = year    



def calculate_persistence(x):
    coefs = coef_dataframe.loc[x,coef_dataframe.columns.str.startswith('ETA')]
    Temp = [s.replace('ETA.',"") for s in coefs.keys().values.tolist()]
    cov_names = Temp[1:Temp.index('M')]
    values = list(coefs)
    b = [mapping[s] for s in cov_names]
    a = values[1:Temp.index('M')]
    return expit((values[0]+sum([x*y for x,y in zip(a,b)]))+ (values[Temp.index('M')]*prev_year_N))

def calculate_colonization(x):
    coefs = coef_dataframe.loc[x,coef_dataframe.columns.str.startswith('DELTA')]
    Temp = [s.replace('DELTA.',"") for s in coefs.keys().values.tolist()]
    cov_names = Temp[1:]
    values = list(coefs)
    b = [mapping[s] for s in cov_names]
    a = values[1:]
    return expit(values[0]+sum([x*y for x,y in zip(a,b)]))

def calculate_occ(x):
    P = calculate_persistence(x)
    C = calculate_colonization(x)
    occ = (P*prev_year_occ)+(C*(1-prev_year_occ))
    return occ

def calculate_growth(x):
    coefs = coef_dataframe.loc[x,coef_dataframe.columns.str.startswith('delta')]
    Temp = [s.replace('delta.',"") for s in coefs.keys().values.tolist()]
    cov_names = Temp[1:]
    values = list(coefs)
    b = [mapping[s] for s in cov_names]
    a = values[1:]
    return np.exp(values[0]+sum([x*y for x,y in zip(a,b)]))


def calculate_abundance(x):
    global prev_year_NC
    global prev_year_occ
    occupancy = calculate_occ(x)
    r = calculate_growth(x)
    NC = r*prev_year_NC
    N = NC*occupancy
    
    NArrayList.append(N)
    prev_year_NC = NC
    prev_year_occ = occupancy
    
    
# Create additional functions for persistence, colonization, growth, occ, abundance, but using the management covariates instead of the baseline ones


def calculate_mgmt_persistence(x):
    coefs = coef_dataframe.loc[x,coef_dataframe.columns.str.startswith('ETA')]
    Temp = [s.replace('ETA.',"") for s in coefs.keys().values.tolist()]
    cov_names = Temp[1:Temp.index('M')]
    values = list(coefs)
    b = [mapping[s] for s in cov_names]
    a = values[1:Temp.index('M')]
    return expit((values[0]+sum([x*y for x,y in zip(a,b)]))+ (values[Temp.index('M')]*prev_year_N_mgmt))

def calculate_mgmt_occ(x):
    P = calculate_mgmt_persistence(x)
    C = calculate_colonization(x)
    occ = (P*prev_year_occ_mgmt)+(C*(1-prev_year_occ_mgmt))
    return occ


def calculate_mgmt_abundance(x):
    global prev_year_NC_mgmt
    global prev_year_occ_mgmt
    occupancy = calculate_mgmt_occ(x)
    r = calculate_growth(x)
    NC = r*prev_year_NC_mgmt
    N = NC*occupancy
    
    NArrayList_mgmt.append(N)
    prev_year_NC_mgmt = NC
    prev_year_occ_mgmt = occupancy




arcpy.AddMessage("Calculating Abundance...")
# Loop across years within model iterations to generate abundance and trend estimates
Spp_Area_Correction = Spp_dict[Spp]
# For table output:
IterationAbundancelist = [] # List (to become dataframe) for final year total abundance in each post. sample, median will be final total abundance estimate
IterationAbundancelist_mgmt = [] # List (to become dataframe) for final year total abundance in each post. \
                                # sample under management scenario, median will be final total abundance under management

# For Raster output
RefYearNArrayList = [] # List to fill with arrays of the reference year abundance, median will be ref year abundance raster
FinalYearNArrayList= [] # List to fill with Arrays of final year abundance, median will be final year abundance raster
FinalYearNArrayList_mgmt= [] # List to fill with Arrays of final year abundance under management, median will be final year abundance raster



TrendList = [] # List to fill with Arrays of trend from each posterior sample, median will be trend raster
TrendList_mgmt = [] # List to fill with Arrays of Trend, median will be trend raster


row_count = len(coef_dataframe)
arcpy.SetProgressor("step", "Processing Posterior Samples...",
                    0, row_count, 1)
for row in range(0,len(coef_dataframe)): # For each posterior sample...
#for row in range(0,5):
    
    # Create blank lists for trend arrays for each year to live in, for this posterior sample
    AnnualTrendList = [] # Create blank lists for trend arrays for each year to live in, for this posterior sample 
    AnnualTrendList_mgmt = []
    
    Nlist = [] # Create blank list for total abundance for each year to live in, for this posterior sample
    Nlist_mgmt = []
    
    NArrayList = [] # Create blank lists for abundance arrays to live in, for this posterior sample
    NArrayList_mgmt = []
    
    # Define initial occupancy
    mapping = mapping2010

    coefs = coef_dataframe.loc[row,coef_dataframe.columns.str.startswith('alpha')]
    Temp = [s.replace('alpha.',"") for s in coefs.keys().values.tolist()]
    cov_names = Temp[1:]
    values = list(coefs)
    b = [mapping2010[s] for s in cov_names]
    a = values[1:]
    #Ψt=2010 = expit{alpha0 + sum(alpha.x × X2010)} #Initial occupancy for each grid cell.
    occ = expit(values[0]+sum([x*y for x,y in zip(a,b)]))
    
    
    
    
    coefs = coef_dataframe.loc[row,coef_dataframe.columns.str.startswith('beta')]
    Temp = [s.replace('beta.',"") for s in coefs.keys().values.tolist()]
    cov_names = Temp[1:]
    values = list(coefs)
    b = [mapping2010[s] for s in cov_names]
    a = values[1:]

    # Calculate initial abundance
    NC2010 = (np.exp(values[0] + sum([x*y for x,y in zip(a,b)]))/Spp_Area_Correction)
    N2010 = NC2010*(occ)
    Nlist.append(np.nansum(N2010))
    NArrayList.append(N2010)
    prev_year_NC = NC2010                                  # Define the initial abundance for use in the following year
    prev_year_N = NArrayList[0]
    prev_year_occ = occ
    

    # Let it rip for the rest of the study period
    for y in yearlist[1:]:
        create_mapping(y) # define which year's rasters to use
        calculate_abundance(row) # Calculate abundance
        
        prev_year_N = NArrayList[-1]
        
        Nlist.append(np.nansum(NArrayList[-1]))
        if y > RefYear:
            TrendArray = NArrayList[-1]/NArrayList[-2]
            AnnualTrendList.append(TrendArray)
    
    
    RefYearNArrayList.append(NArrayList[yearlist.index(RefYear)]) # Grab the Reference Year Abundance Array for this posterior sample
    
    # Fill up the  lists with the previous years' predictions, for future years to compare baseline vs. management 

    Nlist_mgmt = [x for x in Nlist]
    NArrayList_mgmt = [x for x in NArrayList]
    prev_year_NC_mgmt = prev_year_NC      # Define the prev abundance and occupancy for use in the following year for management scenario
    prev_year_N_mgmt = NArrayList[-1]
    prev_year_occ_mgmt = prev_year_occ
    
    #############################################################
    AnnualTrendList_mgmt = [x for x in AnnualTrendList]
##############################################################
    
    ####### Issue: trend for projection years ONLY!!
    
    ###### add some comments in here
    for y in future_years:
        create_mapping(y) # define which year's rasters to use
        calculate_abundance(row) # Calculate abundance
        
              # Define the prev abundance for use in the following year
        prev_year_N = NArrayList[-1]

        create_mgmt_mapping(y)
        calculate_mgmt_abundance(row)
              # Define the prev abundance for use in the following year
        prev_year_N_mgmt = NArrayList_mgmt[-1]
        
        Trend = NArrayList[-1]/NArrayList[-2] # divide the last year's trend array/raster by the previous year's abundance to calculate trend
        Trend_mgmt = NArrayList_mgmt[-1]/NArrayList_mgmt[-2]
        
        AnnualTrendList.append(Trend) # add that year's trend array/raster to the list for the posterior sample
        AnnualTrendList_mgmt.append(Trend_mgmt)
        
        
        Nlist.append(np.nansum(NArrayList[-1]))
        Nlist_mgmt.append(np.nansum(NArrayList_mgmt[-1]))
    
    
    SampleTrend = np.mean(AnnualTrendList, axis=0) # Collapse the annual trend arrays into one array for the posterior sample
    SampleTrend_mgmt = np.mean(AnnualTrendList_mgmt, axis=0)
    
    
    TrendList.append(SampleTrend) 
    TrendList_mgmt.append(SampleTrend_mgmt)
    
    
    FinalYearNArrayList.append(NArrayList[-1])
    FinalYearNArrayList_mgmt.append(NArrayList_mgmt[-1])
    
    IterationAbundancelist.append(Nlist)
    IterationAbundancelist_mgmt.append(Nlist_mgmt)
    arcpy.SetProgressorPosition()
    

TrendDifferenceList = [np.subtract(x,y) for x,y in zip(TrendList_mgmt,TrendList)]

arcpy.ResetProgressor()

#Create the abundance raster and CIs
mx = arcpy.Raster(ModelLayersGDB + '\Lat').extent.XMin
my = arcpy.Raster(ModelLayersGDB + '\Lat').extent.YMin

FinalYearAbundance = np.nanpercentile(FinalYearNArrayList,50,axis=0)
with arcpy.EnvManager(snapRaster=ModelLayersGDB + '\Lat', outputCoordinateSystem = arcpy.Raster(ModelLayersGDB + '\Lat').spatialReference):
    AbunRas = arcpy.NumPyArrayToRaster(FinalYearAbundance, arcpy.Point(mx, my), 1000,1000)
AbunRas.save(ResultsGDB+"\\Abundance_"+str(max(all_years))+"_"+Spp)

FinalYearAbundance = np.nanpercentile(FinalYearNArrayList_mgmt,50,axis=0)
with arcpy.EnvManager(snapRaster=ModelLayersGDB + '\Lat', outputCoordinateSystem = arcpy.Raster(ModelLayersGDB + '\Lat').spatialReference):
    AbunRas = arcpy.NumPyArrayToRaster(FinalYearAbundance, arcpy.Point(mx, my), 1000,1000)
AbunRas.save(ResultsGDB+"\\Abundance_"+str(max(all_years))+"_mgmt_"+Spp)

RefYearAbundance = np.nanpercentile(RefYearNArrayList,50,axis=0)
with arcpy.EnvManager(snapRaster=ModelLayersGDB + '\Lat', outputCoordinateSystem = arcpy.Raster(ModelLayersGDB + '\Lat').spatialReference):
    AbunRas = arcpy.NumPyArrayToRaster(RefYearAbundance, arcpy.Point(mx, my), 1000,1000)
AbunRas.save(ResultsGDB+"\\Abundance_"+str(RefYear)+"_"+Spp)


#Create the trend raster and CIs
    
MedianTrend = np.nanpercentile(TrendList,50,axis=0)
with arcpy.EnvManager(snapRaster=ModelLayersGDB + '\Lat', outputCoordinateSystem = arcpy.Raster(ModelLayersGDB + '\Lat').spatialReference):
    MedianTrendRas = arcpy.NumPyArrayToRaster(MedianTrend, arcpy.Point(mx, my), 1000,1000); MedianTrendRas.save(ResultsGDB+"\\Trend_"+Spp)


# Repeat for management scenario

MedianTrend_mgmt = np.nanpercentile(TrendList_mgmt,50,axis=0)
with arcpy.EnvManager(snapRaster=ModelLayersGDB + '\Lat', outputCoordinateSystem = arcpy.Raster(ModelLayersGDB + '\Lat').spatialReference):
    ManagementTrendRas = arcpy.NumPyArrayToRaster(MedianTrend_mgmt, arcpy.Point(mx, my), 1000,1000); ManagementTrendRas.save(ResultsGDB+"\\Trend_mgmt_"+Spp)

# And for the difference between management and no action

DiffTrend = np.nanpercentile(TrendDifferenceList,50,axis=0)
with arcpy.EnvManager(snapRaster=ModelLayersGDB + '\Lat', outputCoordinateSystem = arcpy.Raster(ModelLayersGDB + '\Lat').spatialReference):
    DiffTrendRas = arcpy.NumPyArrayToRaster(DiffTrend, arcpy.Point(mx, my), 1000,1000); DiffTrendRas.save(ResultsGDB+"\\DiffTrendMedian_"+Spp)

DiffTrendHiCI = np.nanpercentile(TrendDifferenceList,per1,axis=0)
with arcpy.EnvManager(snapRaster=ModelLayersGDB + '\Lat', outputCoordinateSystem = arcpy.Raster(ModelLayersGDB + '\Lat').spatialReference):
    DiffTrendRasHiCI = arcpy.NumPyArrayToRaster(DiffTrendHiCI, arcpy.Point(mx, my), 1000,1000)

DiffTrendLowCI = np.nanpercentile(TrendDifferenceList,per2,axis=0)
with arcpy.EnvManager(snapRaster=ModelLayersGDB + '\Lat', outputCoordinateSystem = arcpy.Raster(ModelLayersGDB + '\Lat').spatialReference):
    DiffTrendRasLowCI = arcpy.NumPyArrayToRaster(DiffTrendLowCI, arcpy.Point(mx, my), 1000,1000); 

StatStrongRaster = Con((DiffTrendRasLowCI>0),1,Con(DiffTrendRasHiCI<0,-1,0))
StatStrongRaster.save(ResultsGDB+"\\StatisticallyStrong_"+Spp)

# Create the csv of abundance / management abundance
arcpy.AddMessage("Generating Abundance csv...")
AbundanceDF = pd.DataFrame(IterationAbundancelist)
AbundanceDF_mgmt = pd.DataFrame(IterationAbundancelist_mgmt)
Median = AbundanceDF.quantile(q=0.5, axis=0)
Median_mgmt = AbundanceDF_mgmt.quantile(q=0.5, axis=0)
Median.name = Spp+" Predicted Abundance - No Action"
Median_mgmt.name = Spp+" Predicted Abundance - Management"
CI_05 = AbundanceDF.quantile(q=(per2/100), axis=0)
CI_05_mgmt = AbundanceDF_mgmt.quantile(q=(per2/100), axis=0)
CI_05.name = "CI_"+str(per2)
CI_05_mgmt.name = "CI_"+str(per2)
CI_95 = AbundanceDF.quantile(q=(per1/100), axis=0)
CI_95_mgmt = AbundanceDF_mgmt.quantile(q=(per1/100), axis=0)
CI_95.name = "CI_"+str(per1)
CI_95_mgmt.name = "CI_"+str(per1)
year = pd.Series(all_years)
year.name = "Year"
TotalAbundanceDF=pd.concat([year,CI_05,Median,CI_95,CI_05_mgmt,Median_mgmt,CI_95_mgmt],axis=1)
TotalAbundanceDF.to_csv(arcpy.Describe(ModelLayersGDB).path+'/'+Spp+'_AbundanceOverTime.csv')

end = time.time()
mins = (end - start) / 60
print(mins)



# Add the output layers to the active map and set their symbology
arcpy.addOutputsToMap = True



arcpy.SetParameter(12,ResultsGDB+"\\Abundance_"+str(RefYear)+"_"+Spp)

arcpy.SetParameter(13,ResultsGDB+"\\Abundance_"+str(max(all_years))+"_"+Spp)

arcpy.SetParameter(14,ResultsGDB+"\\Abundance_"+str(max(all_years))+"_mgmt_"+Spp)

arcpy.SetParameter(15,ResultsGDB+"\\Trend_"+Spp)

arcpy.SetParameter(16,ResultsGDB+"\\Trend_mgmt_"+Spp)

arcpy.SetParameter(17,ResultsGDB+"\\StatisticallyStrong_"+Spp)

outsymSS = arcpy.Describe(ModelLayersGDB).path+"\StatStrongTemplateMgmt.lyrx"
outsymTrend = arcpy.Describe(ModelLayersGDB).path+"\TrendTemplate.lyrx"
outsymAbun = arcpy.Describe(ModelLayersGDB).path+"\AbundanceTemplate.lyrx"

arcpy.SetParameterSymbology(12, outsymAbun)
arcpy.SetParameterSymbology(13, outsymAbun)
arcpy.SetParameterSymbology(14, outsymAbun)

arcpy.SetParameterSymbology(15, outsymTrend)
arcpy.SetParameterSymbology(16, outsymTrend)
arcpy.SetParameterSymbology(17, outsymSS)



