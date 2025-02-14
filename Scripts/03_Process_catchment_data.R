# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process catchment data
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

# Create processed catchement data folder
lapply(paste0(dataDir, "Processed/Catchments"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

# Read England spatVector, and convert to sf object for later processing
england <- readRDS(paste0(dataDir,
                          "Raw/Country_data/England.Rds")) %>%
  st_as_sf

### SEPARATE, FILTER AND SAVE CATCHMENT DATA -----------------------------------

# Read in catchment data from geopackage, with query that selects:
# columns needed, only River water catchment type, and only segments with 1km
# flow accumulation or greater
catchmentData <- read_sf(dsn = paste0(dataDir, "Raw/Flow_data/Flow_data.gpkg"),
                         query = "
                         SELECT *
                         FROM ea_detailed_watersheds
                         WHERE watercat = 'River'
                         ")

# Drop columns not needed, and set standardised crs
catchmentData <- catchmentData %>%
  select (!c(country, watercat, shape_length, shape_area)) %>%
  st_set_crs("EPSG:27700")

# Filter catchments to those which intersect with England country boundary
# N.B. Cannot use original data column "country" as unreliable at Welsh border
catchmentData <- catchmentData %>%
  filter(st_intersects(., england, sparse = FALSE)[,1])

# Save
saveRDS(catchmentData,
        file = paste0(dataDir,
                      "Processed/Catchments/Catchment_data_only.Rds"))

# Remove objects and clear memory
rm(catchmentData)
gc()
