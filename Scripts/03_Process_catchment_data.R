# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process catchment scale data
#
# Script Description:

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(terra)
library(sf)

### DATA MANAGEMENT ------------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Raw/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

# Create processed data folder
lapply(paste0(dataDir, "/Processed/Catchments"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

### SEPARATE, FILTER AND SAVE CATCHMENT DATA -----------------------------------

# Read in catchment data from geopackage, with query that selects:
# columns needed, only River water catchment type, and only segments with 1km
# flow accumulation or greater
catchmentData <- read_sf(dsn = paste0(dataDir, "Raw/Flow_data/Flow_data.gpkg"),
                         query = "
                         SELECT *
                         FROM ea_detailed_watersheds
                         WHERE watercat = 'River'
                         AND country = 'England'
                         ")

# Drop columns not needed
catchmentData <- catchmentData %>%
  select (!c(country, watercat, shape_length, shape_area))

# Save
saveRDS(catchmentData,
        file = paste0(dataDir,
                      "/Processed/Catchments/Catchment_data_only.Rds"))

# Remove objects and clear memory
rm(catchmentData)
gc()

# READ IN CATCHMENT, FERTILISER AND PESTICIDE DATA -----------------------------  

### READ IN CATCHMENT DATA

# Read catchment .Rds
catchmentData <- readRDS(paste0(dataDir,
                                "/Processed/Catchments/Catchment_data_only.Rds"))

### READ IN FERTILISER DATA

# List files (each contains two layers - (1) data layer, (2) uncertainty layer)
fertFiles <- list.files(paste0(dataDir, "Raw/Fertiliser_data/data"),
                        full.names = TRUE)

# Extract the data layer from each file
fertData <- lapply(fertFiles, function(x) { 
  
  # Create spatRaster from first file
  i_R <- rast(x, lyrs = 1)
  
  # Change layer name to remove "uncertainty"
  names(i_R) <- gsub("_prediction_uncertainty_1", "", names(i_R) )
  
  # Return
  return(i_R)

}) %>% 
  
  # Combine together
  rast

### READ IN PESTICIDE DATA

# List files (each contains two layers - (1) data layer, (2) uncertainty layer)
pestFiles <- list.files(paste0(dataDir, "Raw/Pesticide_data/data"),
                        full.names = TRUE)

# Extract the data layer from each file
pestData <- lapply(pestFiles, function(x) { 
  
  # Create spatRaster from first file...
  i_R <- rast(x, lyrs = 1)
  
  # Crop, then extend, to fertiliser data extent (some have different extents)
  i_R <- i_R %>%
    crop(ext(fertData)) %>%
    extend(ext(fertData))
  
  # Change layer names to add pesticide, and remove "_1" from end
  names(i_R) <- paste0("pesticide_", names(i_R) ) %>%
    substring(., 1, nchar(.) - 2)
  
  # Return
  return(i_R)

}) %>% 
  
  # Combine together
  rast

# EXTRACT FERTILISER DATA TO CATCHMENTS ----------------------------------------
# N.B. Warning: this runs overnight

# Weighted sum (missing data is treated as 0)
# N.B. Since each 1x1km data square value is estimated amount applied per 1x1km,
# a weighted sum extract function over each catchment results in estimated 
# amount applied within that catchment
system.time(catchmentFert <- extract(fertData, catchmentData, exact = TRUE, fun = sum,
                         na.rm = TRUE, ID = FALSE, bind = TRUE))

# Save
saveRDS(catchmentFert,
        file = paste0(dataDir,
                      "/Processed/Catchments/Catchment_fertiliser.Rds"))

# EXTRACT PESTICIDE DATA TO CATCHMENTS ----------------------------------------
