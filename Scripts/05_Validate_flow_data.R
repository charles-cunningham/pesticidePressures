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

# Filter data:
# - sampling material to rivers and running freshwater
# - sampling type to freshwater only
# - year to only 2010-2019 to be in similar time frame to flow data
screenData <- screenData %>%
  filter(SMC_DESC == "RIVER / RUNNING SURFACE WATER") %>%
  filter(grepl("FRESHWATER",.$SPT_DESC)) %>%
  filter(year >=2010 & year < 2020)

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

### APPEND FLOW DATA TO SCREEN DATA --------------------------------------------

# Extract simplified pesticide names from flowChemData (remove punctuation)
pesticideNames <- names(flowChemData) %>%
  .[grepl("pesticide_", .)] %>%
  gsub("pesticide_", "", .) %>%
  gsub("\\.", "", .) 

# Add columns to populate in screen data
screenData$PESTICIDE_MOD <- 
  screenData$APPLICATION_MOD <- 
  screenData$CONCENTRATION_MOD <- NA

# Loop through every sampling site (unique screenData geometry)
for(i in unique(st_geometry(screenData))) {

  # IF any flow segments are within 100m (only use these sites)
  if (lengths(st_is_within_distance(i, flowChemData, dist = 100)) > 0) {
    
    # Find nearest feature
    iNearestSegment <- flowChemData[st_nearest_feature(i,
                                                       flowChemData,
                                                       check_crs = FALSE),]
    
    # Find flow accumulation for later 'concentration' calculation
    iMaxflowacc <- iNearestSegment$maxflowacc
    
    # Convert iNearestSegment to tibble, then subset to pesticide columns only
    # (to align with pesticideNames)
    iNearestSegment <- iNearestSegment %>%
      as_tibble %>%
      select(.,contains("pesticide_"))
    
    # IF there are any integer values in iNearestSegment, i.e. not a 
    # maximum value river or a river that has upstream segment outside England
    if (any(!is.na(iNearestSegment))) {
      
      # Loop through every sample j at site i 
      # (screenData rows with same geometry as i)
      for (j in which(lengths(st_equals(screenData, i)) == 1)) {

        # Remove all punctuation from sample j chemical
        nameCheck <- screenData[j,]$Compound_Name %>%
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
        
        # Find whether any of the pesticideNames match, then extract column
        jPesticide <- lapply(pesticideNames, function(x) {
          grepl(x, nameCheck, ignore.case = TRUE)
          }) %>% 
          unlist %>% 
          which

        # If one of the pesticideNames matches exactly one chemical name...
        if (length(jPesticide) == 1) {
          
          # Assign pesticideName
          screenData[j,
                     "PESTICIDE_MOD"] <- pesticideNames[jPesticide]
          # Assign application
          screenData[j,
                     "APPLICATION_MOD"] <- iNearestSegment[, jPesticide]
          # Assign per area application
          screenData[j, 
                     "CONCENTRATION_MOD"] <- screenData[ j, ]$APPLICATION_MOD /
          iMaxflowacc
        }
      }
    }
  }
} # All other values are left as NA

# Save processed screen data
saveRDS(screenData,
        file = paste0(dataDir,
                      "/Processed/Screen/Screen_data.Rds"))

### VALIDATE FLOW DATA WITH SCREEN DATA ----------------------------------------

# Load SCREEN data
flowChemData <- readRDS(paste0(dataDir,
                               "/Processed/Screen/Screen_data.Rds"))

# Considerations
# Filter bad matches post hoc
# LCMS and GCMS values - carry out lm separately
# Remove NAs? - ?

testData <- screenData[!is.na(screenData$APPLICATION_MOD), ]
testData <- testData[testData$method == "LCMS",]

ggplot(data = testData,
       aes(x = log(Concentration), 
           y = log(CONCENTRATION_MOD))) +
  geom_point()
testMod <- lm( Concentration ~ APPLICATION_MOD, data = testData)
summary(testMod)

testModMixed = lmer(Concentration ~ log(CONCENTRATION_MOD) + 
                      (1 | Compound_Name) +
                      (1 | OPCAT_NAME) +
                      (1 | OPCAT_NAME:Sample_Site_ID ),
                    data = testData) 


summary(testModMixed)
hist(log(testData$MODELLED_APPLICATION))
hist(log(testData$MODELLED_CONCENTRATION))
