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
terraOptions(memmax = 0.0001)     
terraOptions()

### DATA MANAGEMENT ------------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Raw/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

# Create processed catchement data folder
lapply(paste0(dataDir, "Processed/Catchments"), function(x) {
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

# Drop columns not needed, and set standardised crs
catchmentData <- catchmentData %>%
  select (!c(country, watercat, shape_length, shape_area)) %>%
  st_set_crs("EPSG:27700")

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

### READ IN FERTILISER DATA TO MEMORY
# As fertData cannot be loaded into memory (reason unknown), use workaround
# below to load the spatRast into memory which massively speeds up later
# extract() funtion

# Save to temporary .Rds file
saveRDS(fertData,
        file = paste0(dataDir, "temp.Rds"))

# Read back in from temporary .Rds file
fertData <- readRDS(paste0(dataDir, "temp.Rds"))

# Delete temporary .Rds file
unlink(paste0(dataDir, "temp.Rds"))

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

# JOIN FERTILISER AND PESTICIDE DATA TOGETHER ----------------------------------

# Join fertliser and pesticide spatRasters into single combined spatRaster
chemData <- c(fertData, pestData)

# Remove separate objects
rm(fertData, pestData)
gc()

# EXTRACT DATA TO CATCHMENTS ---------------------------------------------------
# N.B. Warning: this runs overnight

# Weighted sum (missing data is treated as 0)
# N.B. Since each 1x1km data square value is estimated amount applied per 1x1km,
# a weighted sum extract function for each catchment results in estimated 
# amount applied within that catchment
system.time(
  
  catchmentChemData <- extract(chemData, catchmentData, exact = TRUE,
                              fun = sum, na.rm = TRUE,
                              ID = FALSE, bind = TRUE)
  )

# Save
saveRDS(catchmentChemData,
        file = paste0(dataDir,
                      "/Processed/Catchments/Catchment_chem_data.Rds"))
