# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Join Biosys and flow data
#
# Script Description:

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(sf)

### DIRECTORY MANAGEMENT -------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

# Create processed Biosys data folder
lapply(paste0(dataDir, "Processed/Biosys"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

### LOAD DATA ------------------------------------------------------------------

# Read Biosys data, and convert to sf object
invData_sf <- readRDS(file = paste0(dataDir, "Raw/Biosys/invData.Rds")) %>%
  st_as_sf(.,
           coords = c("FULL_EASTING", "FULL_NORTHING"),
           crs = 27700)

# Read processed flow data
flowChemData <- readRDS(file = paste0(dataDir,
                                      "/Processed/Flow/Flow_aggregated_data.Rds"))

### JOIN PESTICIDE SUMMARY FLOW DATA TO BIOSYS DATA ----------------------------

# Subset for testing 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
testData <- invData_sf[invData_sf$CATCHMENT %>% grep("DEE",.,
                                                     ignore.case = TRUE),]
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Add columns to populate to Biosys data
testData$pesticideLoad <- testData$pesticideDiv <- NA

# Loop through unique geometries (individual sites)
for(i in unique(st_geometry(testData))) {
  
  # Check if any flow segments are within 100m
  if (lengths(st_is_within_distance(i, flowChemData, dist = 100)) > 0) {
    
    # Find nearest feature
    nearestSegment <- flowChemData[st_nearest_feature(i,
                                                      flowChemData,
                                                      check_crs = FALSE), ]
    
    # Transfer aggregated pesticide values from nearest segment to site
    testData[lengths(st_equals(testData, i)) == 1, "pesticideLoad"] <-
      nearestSegment$pesticideLoad
    testData[lengths(st_equals(testData, i)) == 1, "pesticideDiv"] <-
      nearestSegment$pesticideDiv
    
  } # Else, leave as NA
}

### SAVE -----------------------------------------------------------------------

# Save processed invData ready for modelling
saveRDS(testData,
        file = paste0(dataDir,
                      "/Processed/Biosys/invDataSpatial.Rds"))
