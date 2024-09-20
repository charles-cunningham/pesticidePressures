# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Read in spatial data
#
# Script Description:

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(terra)

### DATA MANAGEMENT ------------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Raw/"
# If working locally: "../Data/Raw/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Raw/"

# Manual dataset download folders

# Specify directory for pesticide data
pestDir <-  paste0(dataDir, "Pesticide_data/")
# Specify directory for fertiliser data
fertDir <-  paste0(dataDir, "Fertiliser_data/")
# Specify directory for surface flow data
flowDir <- paste0(dataDir, "Flow_data/")
# Specify directory for land cover data
landDir <-  paste0(dataDir, "Land_cover_data/")

# Automated dataset download folders

# Specify directory for catchment data
catchmentDir <- paste0(dataDir, "Catchment_data/")
# Specify directory for national boundaries data
countryDir <- paste0(dataDir, "Country_data/")
# Specify directory for farming archetype data
farmDir <-  paste0(dataDir, "Farm_archetype_data/")
# Specify directory for topographic data
topoDir <- paste0(dataDir, "Topographic_data/")

# Create directories if they do not exist
c(pestDir,
  fertDir,
  flowDir,
  landDir,
  catchmentDir,
  countryDir,
  farmDir,
  topoDir) %>%
  lapply(., function(x) {
    if (!file.exists(x)) {
      dir.create(x, recursive = TRUE)
    }
  })

### DOWNLOAD PESTICIDE AND FERTILISER DATA [MANUAL] ----------------------------

# Information on data here:l
# https://www.ceh.ac.uk/data/ukceh-land-cover-plus-fertilisers-and-pesticides

# This must be ordered through the EIDC data request:
# https://www.ceh.ac.uk/data-request-form

# Once downloaded:
# (i) Unzip both pesticide and fertiliser data
# (ii) Move to 'pestDir' and 'fertDir' respectively
# (iii) Each directory should be organised in such a way that the respective 
# data and documentation directories sit directly below 'pestDir' or 'fertDir'

### DOWNLOAD FLOW DATA [MANUAL] ------------------------------------------------

# This must be requested by raising a support ticket with Defra Data Services
# https://environment.data.gov.uk/explore/36e7f4d3-61b2-4e64-aaa2-2b85bceb61a9?download=true

# Dataset documentation available:
# https://www.data.gov.uk/dataset/c9dd994d-9649-4041-96d4-cdc0f1a53152/overland-flow-pathways

# Once downloaded:
# (i) .gpkg file is unzipped, then
# (ii) file is renamed 'Flow_data.gpkg', then 
# (iii) file is moved to flowDir [>file.copy("../Flow_data.gpkg", flowDir)]

### DOWNLOAD LAND COVER DATA [MANUAL] ------------------------------------------

# Information on data here:
# https://catalogue.ceh.ac.uk/documents/bb15e200-9349-403c-bda9-b430093807c7

# This must be ordered through the EIDC data request:
# https://order-eidc.ceh.ac.uk/resources/7RDVVS3B/order

# Once downloaded:
# (i) File is unzipped
# (ii) Contents are moved to LandDir
# (iii) Converted to flat structure (i.e. contents of all directories moved to 
# 'landDir' and empty directories removed)

### DOWNLOAD CATCHMENT DATA [AUTOMATED] ----------------------------------------

# Bulk downloads for entire dataset here:
# https://environment.data.gov.uk/catchment-planning

# Dataset documentation available here:
# https://environment.data.gov.uk/catchment-planning/help/usage#the-catchment-data-explorer

# Set download url
url <- "https://environment.data.gov.uk/catchment-planning/England/shapefile.zip"

# Download
download.file(url = url,
              destfile = paste0(catchmentDir, "Catchment_data.zip"),
              sep = "")

# Unzip
unzip(paste0(catchmentDir, "Catchment_data.zip"),
      exdir = catchmentDir)

# Delete .zip files as no longer needed
unlink(paste0(catchmentDir, "Catchment_data.zip"), recursive = TRUE)

### DOWNLOAD AND PROCESS ENGLAND BOUNDARY FOR ANALYSIS [AUTOMATED] -------------

# Information on data here:
# https://www.data.gov.uk/dataset/fac1d8f9-e8ab-47ac-b37a-5640a29634f9/countries-december-2023-boundaries-uk-bsc

# Set download url
url <- "https://open-geography-portalx-ons.hub.arcgis.com/api/download/v1/items/2a0af0a1ecfe473e98e13c7fb8457013/shapefile?layers=0"

# Download
download.file(url = url,
              destfile = paste0(countryDir, "England.zip"),
              sep = "")

# Unzip
unzip(paste0(countryDir, "England.zip"),
      exdir = paste0(countryDir, "EnglandUnzip"),
      overwrite = TRUE)

# Read in boundary
boundaryUK <- paste0(countryDir,
                     "EnglandUnzip/CTRY_DEC_2023_UK_BSC.shp") %>%
  vect

# Filter to England only
boundaryEngland <- boundaryUK[boundaryUK$CTRY23NM == "England",]

# Project to British National Grid
boundaryEngland <- project(boundaryEngland, "EPSG:27700")

# Save to England folder
saveRDS(boundaryEngland,
        file = paste0(countryDir, "England.Rds"))

# Remove files and directories no longer needed
unlink(paste0(countryDir, "England.zip"), recursive = TRUE)
unlink(paste0(countryDir, "EnglandUnzip"), recursive = TRUE)

### DOWNLOAD OTHER BOUNDARIES FOR PLOTTING [AUTOMATED] -------------------------
# N.B. Download all British isles for plotting (UK, Ireland, and Isle of Man)

# Create data frame with name of country, and GADM code
boundaries <- data.frame(country  = c( "UK",  "Ireland", "IsleOfMan"),
                         GADMcode = c( "GBR", "IRL",     "IMN"))

# Loop through each country
for (i in 1:NROW(boundaries)) {
  
  # Assign year and resolution for row i
  country <- boundaries[i,"country"]
  GADMcode <- boundaries[i,"GADMcode"]
  
  # Download country boundary
  boundary <- geodata::gadm(country = GADMcode, level = 0,
                            path = countryDir)
  
  # Project to British National Grid
  boundary <- project(boundary, "EPSG:27700")
  
  # Save country as vector
  saveRDS(boundary,
          file = paste0(countryDir, country, ".Rds"))
  
  # Unlink gadm file
  unlink(paste0(countryDir, "/gadm"),
         recursive = TRUE)
}

### DOWNLOAD FARMING ARCHETYPE DATA [AUTOMATED] --------------------------------
# Information on data here:
# https://catalogue.ceh.ac.uk/documents/3b44375a-cbe6-468c-9395-41471054d0f3

# Download tier 1 archetypes (GB)
download.file(url = "https://catalogue.ceh.ac.uk/datastore/eidchub/3b44375a-cbe6-468c-9395-41471054d0f3/tier1archetypesgb.tif",
              destfile = paste0(farmDir, "tier1.tif"),
              sep = "")

# Download tier 2 archetypes (GB)
download.file(url = "https://catalogue.ceh.ac.uk/datastore/eidchub/3b44375a-cbe6-468c-9395-41471054d0f3/tier2archetypesgb.tif",
              destfile = paste0(farmDir, "tier2.tif"),
              sep = "")

# Download tier 3 archetypes(England)
download.file(url = "https://catalogue.ceh.ac.uk/datastore/eidchub/3b44375a-cbe6-468c-9395-41471054d0f3/tier3archetypesengland.tif",
              destfile = paste0(farmDir, "tier3.tif"),
              sep = "")

### DOWNLOAD TOPOGRAPHIC DATA [AUTOMATED] --------------------------------------

# Bulk downloads for entire dataset here:
# https://osdatahub.os.uk/downloads/open/Terrain50

# Set download url
url <- "https://api.os.uk/downloads/v1/products/Terrain50/downloads?area=GB&format=ASCII+Grid+and+GML+%28Grid%29&redirect"

# Download
download.file(url = url,
              destfile = paste0(topoDir, "Topographic_data.zip"),
              sep = "")

# Unzip data directory and licence.txt to topoDir
unzip(paste0(topoDir, "Topographic_data.zip"),
      exdir = topoDir,
      overwrite = TRUE)

# Delete .zip file as no longer needed
unlink(paste0(topoDir, "Topographic_data.zip"), recursive = TRUE)

### PROCESS TOPOGRAPHIC DATA ---------------------------------------------------
# N.B. Warning: This step takes several hours to run

# List all .zip files
zipFiles <- list.files(paste0(topoDir, "data"),
                       full.names = TRUE,
                       pattern = "\\.zip$",
                       recursive = TRUE)

# Unzip to new folder
lapply(zipFiles, unzip, exdir = paste0(topoDir, "dataUnzipped"))

# List all .asc files in unzipped folder
allAscFiles <- list.files(paste0(topoDir, "dataUnzipped"),
                          full.names = TRUE,
                          pattern = "\\.asc$",
                          recursive = TRUE)

# Create a list of spatRasters from allAscFiles
topoGB_list <- lapply(allAscFiles, rast)

# Merge all spatRasters together
topoGB_R <- do.call("merge", topoGB_list)

# Save merged spatRaster
writeRaster(topoGB_R,
            filename = paste0(topoDir, "topoGB.asc"),
            overwrite = TRUE)

# Delete data and dataUnzipped folder, and remove redundant objects
# For unknown reason need to run unlink x3 (too many files or sub-directories?)
unlink(paste0(topoDir, "data"), recursive = TRUE)
unlink(paste0(topoDir, "data"), recursive = TRUE)
unlink(paste0(topoDir, "dataUnzipped"), recursive = TRUE)
unlink(paste0(topoDir, "dataUnzipped"), recursive = TRUE)
rm(topoGB_list, topoGB_R)
