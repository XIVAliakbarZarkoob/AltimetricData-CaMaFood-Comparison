# Comparison with Altimetric Data
Contains MATLAB scripts for comparing altimetric river water level data with CaMa-Flood model outputs (sfcelv: Water Surface Elevation, and rivout: Discharge) for the Niger and Ganges-Brahmaputra river basins.

The comparison is performed in two ways:

  1. The first temporal components from altimetric data and CaMa-Flood outputs are extracted and compared.
  2. Time series from selected virtual stations are compared to the closest CaMa-Flood output.

The altimetric data is accessible with the link:
https://github.com/XIVAliakbarZarkoob/Altimetric-Data-Check/tree/main/Data

The CaMa-Flood output files required for running this script is stored of google drive accessible with the link:
https://drive.google.com/drive/folders/1gnALukgkgDRnXFPZQBUvw0sqKTq2AO-D?usp=sharing

## Latest MATLAB Scripts

V09: user can specify path to each data in the "Specify Paths" section

default paths are mentioned below where REGION is either Niger or Ganges-Brahmaputra:                                                                     
Dahiti Altimetry Data: './Dahiti/REGION Basin/'                                                                     
Hydroweb Altimetry Data: './Hydroweb/REGION Basin/'                                                                      
CLMS Altimetry data: './CLMS/REGION Basin/'                                                                      
WSE CaMa-Flood Outputs: './REGION_03min_W3RA_AAU/WSE/'                                                                      
Q CaMa-Flood Outputs: './REGION_03min_W3RA_AAU/Q/'

V10: opens a window for each data source to select desired data. a data source can be ignored by closing it's window without selecting any data.



# Runoff Regionalization
Contains MATLAB scripts for regionalizing global runoff data. The script trims global NetCDF runoff files on the Niger and Ganges-Brahmaputra basins and export the results as new NetCDF files.
