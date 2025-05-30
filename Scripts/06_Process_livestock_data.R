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

### READ IN WATERSHED AND LIVESTOCK DATA --------------------------------------

# Read in watershed data
watershedData <- readRDS(paste0(dataDir,
                                "Processed/Watersheds/Watershed_data.Rds"))

# Read in 1km livestock grid
grid1km <- vect(paste0(dataDir,
                        "Raw/Livestock_data/1km_grid.shp")) %>%
  as_tibble()

# Read in livestock data frame
live_df <- read.table(paste0(dataDir,
                               "Raw/Livestock_data/ADAS_1km_2014_Livestock.txt"),
                        header = TRUE, sep = ",")

### PROCESS LIVESTOCK DATA -----------------------------------------------------

# Join grid coordinates to livestock data frame
live_df <- right_join(grid1km,
                      live_df,
                      by = c("UNIQUE" = "Unique"))

# Convert livestock data frame to spatRast (via spatialPoints)
live_R <- vect(live_df, geom=c("X", "Y"), crs = "EPSG:27700") %>%
  rast(., type = "xyz")

# Drop columns not needed
live_R <- subset(live_R, c("UNIQUE", "AREASQKM"), negate = TRUE)

### CROP LIVESTOCK DATA --------------------------------------------------------

# Rasterise all cells that touch watersheds which intersect with England
watershed_R <- vect(watershedData) %>%
  rasterize(., live_R, touches = TRUE)

# Create extent object to crop the livestock spatRasts to
cropExtent <- ext(watershed_R)

# Crop and extend livestock spatRast to cropExtent 
live_R <- crop(live_R, cropExtent) %>%
  extend(., cropExtent)

### SUBSET LIVESTOCK DATA ---------------------------------------------------

# Select total numbers of livestock
live_R <- subset(live_R, c("K299", "L98", "M98", "N98"))

# Rename to inutuitive names
names(live_R) <- c("cattle", "pigs", "sheep", "poultry")

# Ensure spatRast is in memory
toMemory(live_R)

### INTERPOLATE 1KM LIVESTOCK DATA  --------------------------------------------
# Use gstat inverse weighting approach in case any gaps

# For each layer name, i.e. each chemical
liveInterp <- lapply(names(live_R), function(i) {
  
  # Create 'x,y,z' data frame of chemical i (coordinates and values) 
  liveData_df <- as.points(live_R[[i]]) %>%
    data.frame(crds(.), . )
  
  # Rename livestock column to standardized name to avoid errors- "z"
  names(liveData_df) <- c("x", "y", "z")
  
  # Create a basic inverse distance weighted gstat interpolation formula
  iInterpModel <- gstat(formula = z ~ 1, # Only use location
                        locations = ~ x + y, # Location based on x, y
                        data = liveData_df, 
                        nmax = 5) # Only use nearest 5 points
  
  # Interpolate using gstat formula
  iInterp <- interpolate(watershed_R,
                         iInterpModel,
                         na.rm = TRUE)[["var1.pred"]]
  
  # Change name back
  names(iInterp) <- i
  
  # Return
  return(iInterp)
  
}) %>% rast() # Join all layers together into one spatRast

### EXTRACT DATA TO WATERSHEDS -------------------------------------------------
# N.B. Warning: this runs overnight

# Weighted sum (missing data is treated as 0) of load
# N.B. Since each 1x1km livestock_R value has units of number of livestock,
# a weighted sum extract function for each watershed results in
# estimated number within that watershed
watershedLiveData <- terra::extract(liveInterp, watershedData,
                                    exact = TRUE,
                                    fun = sum,
                                    na.rm = TRUE,
                                    ID = TRUE) # N.B. Same as watershedData ID

# Save
saveRDS(watershedLiveData,
        file = paste0(dataDir,
                      "Processed/Watersheds/Watershed_live_data.Rds"))
