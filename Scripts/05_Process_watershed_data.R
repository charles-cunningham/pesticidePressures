# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process watershed data
#
# Script Description:

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(terra)
library(sf)

# Set terra options to speed up
terraOptions(memfrac = 0.9)

### DATA MANAGEMENT ------------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

# Create processed watershed data folder
lapply(paste0(dataDir, "Processed/Watersheds"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

# Read England spatVector, and convert to sf object for later processing
england <- readRDS(paste0(dataDir,
                          "Raw/Country_data/England.Rds")) %>%
  st_as_sf

### SEPARATE, FILTER AND SAVE WATERSHED DATA -----------------------------------

# Read in watershed data from geopackage, with query that selects only 
# 'River' water catchment type
watershedData <- read_sf(dsn = paste0(dataDir, "Raw/Flow_data/Flow_data.gpkg"),
                         query = "
                         SELECT *
                         FROM ea_detailed_watersheds
                         WHERE watercat = 'River'
                         ")

# Drop columns not needed, and set standardised crs
watershedData <- watershedData %>%
  select (!c(country, watercat, shape_length, shape_area)) %>%
  st_set_crs("EPSG:27700")

# Filter watersheds to those which intersect with England country boundary
# N.B. Cannot use watershedData data column "country" as unreliable at Welsh border
watershedData <- watershedData %>%
  filter(st_intersects(., england, sparse = FALSE)[,1])

# Create unique ID column as does not exist
watershedData$ID <- 1:NROW(watershedData)

# Save
saveRDS(watershedData,
        file = paste0(dataDir,
                      "Processed/Watersheds/Watershed_data.Rds"))

# Remove objects and clear memory
rm(watershedData)
gc()
