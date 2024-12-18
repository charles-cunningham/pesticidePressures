# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Validate overland flow data
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

# Create processed flow data folder
lapply(paste0(dataDir, "Processed/Screen"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

### LOAD DATA ------------------------------------------------------------------

# Load flow data
flowChemData <- readRDS(paste0(dataDir,
                               "Processed/Flow/Flow_chem_data.Rds"))

# Load LC-MS (liquid chromatography-mass spectrometry) screen data
lcms <- paste0(dataDir,
               "Raw/Monitoring_screen_data/LCMS Target and Non-Targeted Screening.csv") %>%
  read.csv()

# Load GC-MS (gas chromatography-mass spectrometry) screen data

gcms <- paste0(dataDir,
               "Raw/Monitoring_screen_data/GCMS Target and Non-Targeted Screening.csv") %>%
  read.csv()

# Join together for processing (can still be separated by 'method' column)
# N.B. From README metadata - "when overlap occurs (ie. same compounds can be 
# detected by both methods) the LC-MS is able to detect using a much lower limit 
# of detection than GC-MS."
screenData <- rbind(lcms, gcms)

# Remove redundant objects
rm(lcms, gcms)

### PROCESS SCREEN DATA --------------------------------------------------------

# FILTER DATA

# Filter:
# - sampling material to rivers and running freshwater
# - sampling type to freshwater only
# - year to only 2010-2019 to be in similar time frame to flow data
screenData <- screenData %>%
  filter(SMC_DESC == "RIVER / RUNNING SURFACE WATER") %>%
  filter(grepl("FRESHWATER",.$SPT_DESC)) %>%
  filter(year >=2010 & year < 2020)

# REMOVE COLUMNS NOT NEEDED
screenData <- screenData %>%
  select(!c(SMPT_TYPE,
            SPT_DESC,
            SAMP_MATERIAL,
            SMC_DESC,
            SAMP_PURPOSE_CODE,
            ARE_CODE,
            Latitude,
            Longitude,
            COUNTRY))

# CREATE SPATIAL OBJECT

screenData <- st_as_sf(screenData,
                       coords = c("SMPT_EASTING", "SMPT_NORTHING"),
                       crs = 27700)

### APPEND FLOW DATA TO SCREEN DATA ------------------------------------------

# Extract pesticide names
pesticides <- names(flowChemData) %>%
  .[grepl("pesticide_", .)] %>%
  gsub("pesticide_", "", .)

# Add columns to populate in screen data
screenData$MODELLED_APPLICATION <- screenData$MODELLED_CONCENTRATION <- NA

# Loop through every sampling site (unique screen data geometry)
for(i in unique(st_geometry(screenData))) {
  
  #i <- unique(st_geometry(screenData))[[102]]
  
  # Check if any flow segments are within 100m (only use these sites)
  if (lengths(st_is_within_distance(i, flowChemData, dist = 100)) > 0) {
    
    # Find nearest feature, and extract pesticide values
    iNearestSegment <- flowChemData[st_nearest_feature(i,
                                                      flowChemData,
                                                      check_crs = FALSE), ]
    
    # Find flow accumulation for later application per area calculation
    iMaxflowacc <- iNearestSegment$maxflowacc
    
    # Select pesticide layers
    iNearestSegment <- iNearestSegment %>%
      select(.,contains("pesticide_")) %>%
      tibble
    
    # Check if there are any integer values in nearest segment, i.e. not a 
    # maximum value river or a river that has upstream segment outside England
    if (any(!is.na(iNearestSegment))) {
      
      # Loop through every sample at site (screen data rows with same geometry)
      for (j in which(lengths(st_equals(screenData, i)) == 1)) {

        jPesticide <- grep(screenData[j,]$Compound_Name,
                           pesticides)
        
        if (length(jPesticide) > 0) {
          
        # Assign application
        screenData[j,
                   "MODELLED_APPLICATION"] <- iNearestSegment[,jPesticide]
        
        # Assign per area application
        screenData[j, 
                   "MODELLED_CONCENTRATION"] <- screenData[j,
                                                           "MODELLED_APPLICATION"] /
          iMaxflowacc
        }
      }
    }
    
  } # Else, leave as NA
}


