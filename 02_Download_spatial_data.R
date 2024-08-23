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

# Set download url
url <- "https://api.os.uk/downloads/v1/products/Terrain50/downloads?area=GB&format=ASCII+Grid+and+GML+%28Grid%29&redirect"

# Download
download.file(url = url,
              destfile = paste0(topoDir, "Topographic_data.zip"),
              sep = "")

# Unzip
unzip(paste0(topoDir, "Topographic_data.zip"),
      exdir = topoDir)

# Delete .zip files as no longer needed
unlink(paste0(topoDir, "Topographic_data.zip"), recursive = TRUE)

### PROCESS TOPOGRAPHIC DATA ---------------------------------------------------

# List all .zip files
zipFiles <- list.files(paste0(topoDir, "data"),
                   full.names = TRUE,
                   pattern = "\\.zip$",
                   recursive = TRUE)

# Unzip to new folder
lapply(zipFiles, unzip, exdir = paste0(topoDir, "dataUnzipped"))

# List all .asc files in unzipped folder


# Read and merge these files


# Save under topoDir


# Unlist dataUnzipped folder

list.files(topoDir)
