# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process screen data
#
# Script Description: Process screen data for use in validating flow data in
# next script

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
                               "Processed/Flow/Flow_data_all.Rds")) %>%
  # Remove sites with upstream area greater than 500km^2 (500,000,000m^2) 
  # to avoid increasing uncertainty with largest areas
  filter((totalArea * 25*25) < 500000000)

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

# Filter data:
# - sampling material to rivers and running freshwater
# - sampling type to freshwater only
# - year to only 2010-2019 to be in similar time frame to flow data
screenData <- screenData %>%
  filter(SMC_DESC == "RIVER / RUNNING SURFACE WATER") %>%
  filter(grepl("FRESHWATER",.$SPT_DESC)) %>%
  filter(year >= 2010 & year < 2020)

# Remove columns not needed
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

# Create spatial object
screenData <- st_as_sf(screenData,
                       coords = c("SMPT_EASTING", "SMPT_NORTHING"),
                       crs = 27700)

### CHECK FLOW DATA AND SCREEN DATA NAME MATCHES AND FILTER --------------------

# Create simplified pesticide names from flowChemData (remove punctuation)
flowPestNames <- names(flowChemData) %>%
  .[grepl("pesticide_", .)] %>%
  gsub("pesticide_", "", .) %>%
  gsub("\\.", "", .) 

# Create simplified screenData comound name column (remove punctuation)
screenData$nameCheck <- screenData$Compound_Name %>%
  gsub("\\(", "", .) %>%
  gsub("\\)", "", .) %>%
  gsub("\\[", "", .) %>%
  gsub("\\]", "", .) %>%
  gsub("\\-", "", .) %>%
  gsub("\\,", "", .) %>%
  gsub("\\.", "", .) %>%
  gsub("\\/", "", .) %>%
  gsub("\\&", "", .) %>%
  gsub("\\'", "", .) %>%
  gsub(" ", "", .)

# Check exact matches (only use these in analysis)
matchNames <- flowPestNames[ flowPestNames %in% screenData$nameCheck ] %>%
  unique %>%
  sort

# Filter screen data to only rows containing pesticides found in flowData
screenData <- screenData %>%
  filter(nameCheck %in% matchNames)

# Check what proportion of pesticides we are validating
(length(matchNames) / length(flowPestNames) * 100) %>%
  round(., digits = 2) %>%
  paste(., "% of total Land Cover Plus pesticides included in validation") %>%
  print

### APPEND FLOW DATA TO SCREEN DATA --------------------------------------------

# Add columns to populate with modelled data in screenData
screenData$PESTICIDE_MOD <- 
  screenData$APPLICATION_MOD <- 
  screenData$APP_PER_AREA_MOD <- NA

# Loop through every sampling site (use unique screenData geometry)
for(i in unique(st_geometry(screenData))) {
  
  # Find nearest flowData segment to screenData sampling site
  iNearestSegment <- flowChemData[st_nearest_feature(i,
                                                     flowChemData,
                                                     check_crs = FALSE),]
  
  # IF nearest feature is within 50m (only use these sites)
  if (lengths(st_is_within_distance(i, iNearestSegment, dist = 50)) == 1) {
    
    # Find flow accumulation for later per area calculation
    iMaxflowacc <- iNearestSegment$maxflowacc
    
    # Convert iNearestSegment to tibble, then subset to pesticide columns only
    # (to align with flowPestNames)
    iNearestSegment <- iNearestSegment %>%
      as_tibble %>%
      select(.,contains("pesticide_"))
    
    # IF there are any integer values in iNearestSegment, i.e. not a 
    # maximum flowaccl value river or river with upstream segment outside England
    if (any(!is.na(iNearestSegment))) {
      
      # Loop through every sample j at site i 
      # (screenData rows with same geometry as i)
      for (j in which(lengths(st_equals(screenData, i)) == 1)) {
        
        # Find which of the flowPestNames match, then extract column number
        # N.B. Order of flowPestNames matches pesticide columns, so can use
        # order to find column
        jPesticide <- flowPestNames %in% screenData$nameCheck[j] %>%
          which()
        
        # Assign pesticide name for sense check
        screenData[j,
                   "PESTICIDE_MOD"] <- flowPestNames[jPesticide]
        
        # Assign total estimated application per year
        screenData[j,
                   "APPLICATION_MOD"] <-
          iNearestSegment[, jPesticide]
        
        # Assign total estimated application per year per m^2
        screenData[j,
                   "APP_PER_AREA_MOD"] <-
          screenData[j,]$APPLICATION_MOD / iMaxflowacc
      }
    }
  }
} # All other values are left as NA

# Remove screenData nameCheck column
screenData$nameCheck <- NULL

# Save processed screen data
saveRDS(screenData,
        file = paste0(dataDir,
                      "/Processed/Screen/Screen_data.Rds"))
