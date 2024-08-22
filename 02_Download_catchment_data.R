# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process EA catchment data
#
# Script Description:

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(terra)

### DATA MANAGEMENT ------------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Raw/SpatialData/CatchmentData/"
# If working locally: "../Data/Raw/Spatial_data/CatchmentData/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Raw/SpatialData/CatchmentData/"

# Create directory if it doesn't exist
if (!file.exists(dataDir)) {
  dir.create(dataDir, recursive = TRUE)
}

### DOWNLOAD CATCHMENT DATA -------------------------------------------------------

# Bulk downloads for entire dataset here:
# https://environment.data.gov.uk/catchment-planning

# Dataset documentation available here:
# https://environment.data.gov.uk/catchment-planning/help/usage#the-catchment-data-explorer

# Set download url
url <- "https://environment.data.gov.uk/catchment-planning/England/shapefile.zip"

# Set file name
file_name <- "CatchmentData"

# Download
download.file(url = url, destfile = paste0(dataDir, file_name, ".zip"), sep = "")

# Unzip
unzip(paste0(dataDir, "/", file_name, ".zip"), exdir = dataDir)

# Delete .zip files as no longer needed
unlink(paste0(dataDir,
              "CatchmentData.zip"))

test <- vect(paste0(dataDir, "WFD_River_Water_Bodies_Cycle_3.shp"))

plot(test)
