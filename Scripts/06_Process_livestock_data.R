# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process livestock data
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

# Create processed watershed data folder
lapply(paste0(dataDir, "Processed/Watersheds"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

### READ IN WATERSHED AND LAND COVER DATA --------------------------------------

# Read in watershed data

watershedData <- readRDS(paste0(dataDir,
                                "Processed/Watersheds/Watershed_data.Rds"))

# Read in livestock data
cattle <- rast(paste0(dataDir,
                      "Raw/Livestock_data/APHA_LDDG_Cattle_Pop_2015.tif"))
pigs <- rast(paste0(dataDir,
                    "Raw/Livestock_data/APHA_LDDG_Pig_Pop_2016_2017.tif"))
poultry <- rast(paste0(dataDir,
                       "Raw/Livestock_data/APHA_LDDG_Poultry_Pop_2016.tif"))
sheep <- rast(paste0(dataDir,
                     "Raw/Livestock_data/APHA_LDDG_Sheep_Pop_2015_2016.tif"))

plot(cattle)

# Read in 1km grid
grid1km <- vect(paste0(dataDir,
                        "Raw/Livestock_data/1km_grid.shp")) %>%
  as_tibble()

# Read in livestock data
live_df <- read.table(paste0(dataDir,
                               "Raw/Livestock_data/ADAS_1km_2014_Livestock.txt"),
                        header = TRUE, sep = ",")

# Join geometry to livestock data
live_df <- right_join(grid1km,
                      live_df,
                      by = c("UNIQUE" = "Unique"))

# Convert livestock data to raster (via spatialPoints)
live_R <- vect(live_df, geom=c("X", "Y"), crs = "EPSG:27700") %>%
  rast(., type = "xyz")

### PROCESS LIVESTOCK DATA -----------------------------------------------------

# Create extent object to crop the livestock spatRasts to
cropExtent <- vect(watershedData) %>%
  rasterize(., live_R, touches = TRUE) %>%
  ext()

# Crop and extend livestock spatRasts to the cropExtent extent
cattle <- crop(live_R, cropExtent) %>%
  extend(., cropExtent)

# Also manually set extent to remove small inconsistencies in extent
ext(live_R) <- ext(cropExtent)


select and rename

# Merge into single spatRast object
livestock <- c(cattle, pigs, poultry, sheep)

### EXTRACT DATA TO WATERSHEDS -------------------------------------------------
# N.B. Warning: this runs overnight

# Weighted sum (missing data is treated as 0) of load
# N.B. Since each 1x1km livestock_R value has units of number of livestock,
# a weighted sum extract function for each watershed results in
# estimated number within that watershed
watershedLiveData <- terra::extract(livestock, watershedData,
                                    exact = TRUE,
                                    fun = sum,
                                    na.rm = TRUE,
                                    ID = TRUE) # N.B. Same as watershedData ID

# Save
saveRDS(watershedLiveData,
        file = paste0(dataDir,
                      "Processed/Watersheds/Watershed_live_data.Rds"))
