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

# Specify directory for catchment data
catchmentDir <- paste0(dataDir, "Catchment_data/")

# Specify directory for surface flow data
flowDir <- paste0(dataDir, "Flow_data/")

# Specify directory for national boundaries data
countryDir <- paste0(dataDir, "Country_data/")

# Specify directory for topographic data
topoDir <- paste0(dataDir, "Topographic_data/")

# Create directories if they do not exist
c(catchmentDir, flowDir, countryDir, topoDir) %>%
  lapply(., function(x) {
    if (!file.exists(x)) {
      dir.create(x, recursive = TRUE)
    }
  })

### DOWNLOAD CATCHMENT DATA ----------------------------------------------------

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

### DOWNLOAD FLOW DATA ---------------------------------------------------------

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

### DOWNLOAD TOPOGRAPHIC DATA --------------------------------------------------

# Bulk downloads for entire dataset here:
# https://osdatahub.os.uk/downloads/open/Terrain50

# Set download url
url <- "https://api.os.uk/downloads/v1/products/Terrain50/downloads?area=GB&format=ASCII+Grid+and+GML+%28Grid%29&redirect"

# Download
download.file(url = url,
              destfile = paste0(topoDir, "Topographic_data.zip"),
              sep = "")

# Unzip
unzip(paste0(topoDir, "Topographic_data.zip"),
      exdir = topoDir,
      overwrite = TRUE)

# Delete .zip files as no longer needed
unlink(paste0(topoDir, "Topographic_data.zip"), recursive = TRUE)

### DOWNLOAD AND PROCESS ENGLAND BOUNDARY FOR ANALYSIS -------------------------

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
      exdir = paste0(countryDir, "England/EnglandUnzip"),
      overwrite = TRUE)

# Read in boundary
boundaryUK <- vect(paste0(countryDir, "England/EnglandUnzip/CTRY_DEC_2023_UK_BSC.shp"))

# Filter to England only
boundaryEngland <- boundaryUK[boundaryUK$CTRY23NM == "England",]

# Save to England folder
writeVector(boundaryEngland,
            filename = paste0(countryDir, "England/boundary.shp"))

# Remove files and directories no longer needed
unlink(paste0(countryDir, "England.zip"), recursive = TRUE)
unlink(paste0(countryDir, "England/EnglandUnzip"), recursive = TRUE)

list.files(paste0(countryDir, "England"))
### DOWNLOAD OTHER BOUNDARIES FOR PLOTTING -------------------------------------
# N.B. Download all British isles for plotting (UK, Ireland, and Isle of Man)

# NATIONAL BOUNDARIES -------------------------------
# (UK, Ireland, and Isle of Man)

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
                            path = paste0(countryDir,
                                          country ))
  
  # Reproject country to BNG
  boundary <- project(boundary, bng)
  
  # Assign boundary to corresponding country
  assign(country, boundary)
  
  # If file not already saved, then write to file
  if(!file.exists(paste0(countryDir,
                         country, "/", GADMcode, ".shp" ))) {
    
    # Save country as vector
    writeVector(boundary, 
                filename = paste0(countryDir,
                                  country ),
                layer = GADMcode,
                filetype = "ESRI Shapefile")
  }
}

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

# Delete dataUnzipped folder, and remove redundant objects
unlink(paste0(topoDir, "dataUnzipped"), recursive = TRUE)
rm(topoGB_list, topoGB_R)
