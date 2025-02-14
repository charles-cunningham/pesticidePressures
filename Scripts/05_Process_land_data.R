# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process land cover data
#
# Script Description:

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(terra)
library(sf)
library(gstat)

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

### READ IN CATCHMENT AND LAND COVER DATA --------------------------------------

### READ IN CATCHMENT DATA

# Read catchment .Rds
catchmentData <- readRDS(paste0(dataDir,
                                "Processed/Catchments/Catchment_chem_data.Rds"))

### READ IN LAND COVER DATA

# Read land cover spatRast
lcm2015 <- rast(paste0(dataDir, "Raw/Land_cover_data/lcm2015gb25m.tif"))[[1]]

# Rename land cover spatRast
names(LCM2021) <- "Identifier"

# Specify the LCM classes, add on "No Data"
classLCM <- c(
  "Deciduous woodland",
  "Coniferous woodland",
  "Arable",
  "Improved grassland",
  "Neutral grassland",
  "Calcareous grassland",
  "Acid grassland",
  "Fen",
  "Heather",
  "Heather grassland",
  "Bog",
  "Inland rock",
  "Saltwater",
  "Freshwater",
  "Supralittoral rock",
  "Supralittoral sediment",
  "Littoral rock",
  "Littoral sediment",
  "Saltmarsh",
  "Urban",
  "Suburban",
  "No data" # This is added (not in original LCM Classes)
)

# Create LCM identifier/class data frame
# N.B. "No Data" identifier is 'NA'
LCM_df <- data.frame("Identifier" = c(1:(length(classLCM) - 1),
                                      NA),
                     "Class" = classLCM)

# EXTRACT COVERAGE -----------------------------------------


# Create progress bar
progressBar = txtProgressBar(
  min = 0,
  max = NROW(SSSI_df),
  initial = 0,
  style = 3
)

# Start loop iterating through every SSSI
for (j in 1:NROW(SSSI_df)) {
  # Find all SSSI vect objects with same unique ID (part of SSSI 'j')
  SSSI_j <- subset(SSSI,
                   subset = values(SSSI)[, uniqueID] == SSSI_df[j, uniqueID])
  
  # Extract all 25m cells for each land cover class present for SSSI 'j'
  SSSIcells <- terra::extract(LCM2021, SSSI_j)
  
  # Count number of cells for each class
  # N.B. some classes may not be included as count is 0
  SSSIcount <- count(SSSIcells, Identifier)
  
  # Create percentage cover column from count
  SSSIcount$PercentCover <- SSSIcount$n / sum(SSSIcount$n) * 100
  
  # Add class names by joining coverage values to LCM data frame
  SSSIcount <- full_join(LCM_df, SSSIcount, by = "Identifier") %>%
    replace_na(list(PercentCover = 0)) # Convert PercentCover NAs to 0
  
  # Add proportion coverage for each class to SSSI data frame (row j)
  SSSI_df[j, colNumsLCM] <- SSSIcount$PercentCover
  
  # Iterate progress bar
  setTxtProgressBar(progressBar, j)
  
}

# Close progress bar
close(progressBar)


# Save
saveRDS(catchmentLandData,
        file = paste0(dataDir,
                      "Processed/Catchments/Catchment__data.Rds"))

