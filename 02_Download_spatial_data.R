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

# Specify directory for topgraphic data
topoDir <- paste0(dataDir, "Topographic_data/")

# Create directories if they do not exist
c(catchmentDir, topoDir) %>%
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

### DOWNLOAD TOPOGRAPHIC DATA --------------------------------------------------

# Bulk downloads for entire dataset here:
# https://osdatahub.os.uk/downloads/open/Terrain50

url <- "https://api.os.uk/downloads/v1/products/Terrain50/downloads?area=GB&format=ASCII+Grid+and+GML+%28Grid%29&redirect"

# Set file name
file_name <- "CatchmentData"

# Download
download.file(url = url, destfile = paste0(dataDir, file_name, ".zip"), sep = "")

# Unzip
unzip(paste0(dataDir, "/", file_name, ".zip"), exdir = dataDir)

# Delete .zip files as no longer needed
unlink(paste0(dataDir,
              "CatchmentData.zip"))


