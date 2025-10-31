# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Join Biosys and upstream and site data
#
# Script Description: Runs overnight

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
invData_sf <- readRDS(file = paste0(dataDir, "Processed/Biosys/invData.Rds")) %>%
  st_as_sf(.,
           coords = c("FULL_EASTING", "FULL_NORTHING"),
           crs = 27700)

# Read processed flow data
flowData <- readRDS(file = paste0(dataDir,
                                      "/Processed/Flow/Flow_aggregated_data.Rds"))

# List flowData columns to join to invData_sf
upstreamData <- as_tibble(flowData) %>%
  select(!c("ea_wb_id", "opcatch", "mancatch", "ea_area", "wbname",               
                 "rbd", "meanslope", "startz", "endz", "permid",          
                 "maxflowacc", "shape", "withinEngland")) %>%
  names

### JOIN PESTICIDE SUMMARY FLOW DATA TO BIOSYS DATA ----------------------------

# Add columns to populate to Biosys data
invData_sf[, upstreamData] <- NA

# Find unique Site IDs
sites <- unique(invData_sf$SITE_ID)

# Loop through individual sites
for(site in sites) {
  
  # Find site rows in invData_sf
  siteRows <- which(invData_sf$SITE_ID == site)
  
  # Find site geometries for all site rows, then take unique values
  siteGeometry <- st_geometry(invData_sf[siteRows,]) %>%
    unique()
  
  # Make sure only one geometry per site (only one unique geometry value)
  if (length(siteGeometry) > 1) {
    stop("Site has more than one geometry")   
  }

  # Find waterbody ID for site rows...
  iWaterbody <- invData_sf[siteRows, "WFD_WATERBODY_ID"] %>%
    # ...then take unique values
    unique %>% 
    .[[1]]
  
  # Filter flowChemData to waterbodies for i (both datasets have this)
  waterbodyFlowData <- flowData[flowData$ea_wb_id %in% iWaterbody,]
  
  # Find any waterbody segments nearby (within 50m)
  nearbySegments <- st_is_within_distance(siteGeometry[[1]], 
                                          waterbodyFlowData,
                                          dist = 50)[[1]] %>%
    waterbodyFlowData[. ,]

  # If any flow segments are within 50m...
  if (NROW(nearbySegments) > 0) {
 
    # Find nearest feature from nearby segments
    nearestSegment <- nearbySegments[st_nearest_feature(siteGeometry[[1]],
                                                        nearbySegments,
                                                        check_crs = FALSE), ] %>%
      as_tibble()
    
    # Transfer aggregated upstreamData values from nearest segment to site
    invData_sf[siteRows, upstreamData] <- nearestSegment[, upstreamData]
    
    } # Else, leave as NA
  
  # PRINT UPDATE

  # Every 100 sites...
  if (which(sites %in% site) %% 100 == 0) {
    
    # Print update
    paste(which(sites %in% site),
          "of",
          length(sites),
          "complete") %>%
      print()
  }
}

### JOIN SITE DATA -------------------------------------------------------------

# Collect relevent site variable file names
siteFiles <- list.files(paste0(dataDir, "/Raw/Site_data/data"), full.names = TRUE) %>%
  grep("sitevariables.csv", ., value = TRUE)

# Read and bind files
siteData <- lapply(siteFiles, read.csv) %>% 
  bind_rows()

# Join data
invData_sf <- left_join(invData_sf, siteData,
                        by = c("SITE_ID" = "BIO_SITE_ID"))

### SAVE -----------------------------------------------------------------------

# Save processed invData ready for modelling
saveRDS(invData_sf,
        file = paste0(dataDir,
                      "/Processed/Biosys/invDataSpatial.Rds"))
