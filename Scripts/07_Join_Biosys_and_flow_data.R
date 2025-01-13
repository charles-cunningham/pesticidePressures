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

# Add columns to populate to Biosys data
invData_sf$pesticideLoad <- invData_sf$pesticideDiv <- NA

# Loop through unique geometries (individual sites)
for(i in unique(st_geometry(invData_sf))) {

  # Find waterbody ID of i; first identify rows that match geometry
  iWaterbody <- invData_sf[lengths(st_equals(invData_sf, i)) == 1,
                           "WFD_WATERBODY_ID"] %>%  # 
    # Then take unique waterbody ID (summarise many rows with same ID)
    unique %>% 
    .[[1]]

  # Filter flowChemData to waterbody i (both datasets have this information)
  waterbodyFlowData <- flowChemData[flowChemData$ea_wb_id == iWaterbody,]
  
  # Find any segments nearby (within 100m)
  nearbyGeometry <- st_is_within_distance(i, waterbodyFlowData, dist = 100) %>%
    .[[1]] %>%
    flowChemData[.,]

  # If any flow segments are within 100m...
  if (NROW(nearbyGeometry) > 0) {
 
    # Find nearest feature
    nearestSegment <- nearbyGeometry[st_nearest_feature(i,
                                                        nearbyGeometry,
                                                        check_crs = FALSE), ]
    
    # Transfer aggregated pesticide values from nearest segment to site
    invData_sf[lengths(st_equals(invData_sf, i)) == 1, "pesticideLoad"] <-
      nearestSegment$pesticideLoad
    invData_sf[lengths(st_equals(invData_sf, i)) == 1, "pesticideDiv"] <-
      nearestSegment$pesticideDiv
    
  } # Else, leave as NA
}

### SAVE -----------------------------------------------------------------------

# Save processed invData ready for modelling
saveRDS(invData_sf,
        file = paste0(dataDir,
                      "/Processed/Biosys/invDataSpatial.Rds"))
