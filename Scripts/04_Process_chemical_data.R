# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process chemical application data
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

# Create processed watershed data folder
lapply(paste0(dataDir, "Processed/Watersheds"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

### READ IN WATERSHED, FERTILISER AND PESTICIDE DATA ---------------------------  

### READ IN WATERSHED DATA

# Read watershed .Rds
watershedData <- readRDS(paste0(dataDir,
                                "Processed/Watersheds/Watershed_data.Rds"))

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

### JOIN FERTILISER AND PESTICIDE DATA TOGETHER --------------------------------

# Join fertliser and pesticide spatRasters into single combined spatRaster
chemData <- c(fertData, pestData)

# Save layer names
chemNames <- names(chemData)

# Optional: Read spatRast to memory (speeds up later extraction)
chemData <- toMemory(chemData)

# Remove separate objects
rm(fertData, pestData)
gc()

### INTERPOLATE RASTER DATA  ---------------------------------------------------

# Rasterise all cells that touch watersheds which intersect with England
# N.B. There are some gaps in the Land Cover Plus data to interpolate such as
# watersheds along England/Wales border
watershed_R <- vect(watershedData) %>%
  rasterize(., chemData, touches = TRUE)

# For each layer name, i.e. each chemical
chemDataInterp <- lapply(chemNames, function(i) {
  
  # Create 'x,y,z' data frame of chemical i (coordinates and values) 
  chemData_df <- as.points(chemData[[i]]) %>%
    data.frame(crds(.), . )
  
  # Rename chemical column to standardized name to avoid errors- "z"
  names(chemData_df) <- c("x", "y", "z")
  
  # Create a basic inverse distance weighted gstat interpolation formula
  iInterpModel <- gstat(formula = z ~ 1, # Only use location
                        locations = ~ x + y, # Location based on x, y
                        data = chemData_df, 
                        nmax = 5) # Only use nearest 5 points
  
  # Interpolate using gstat formula for chemical i
  iInterp <- interpolate(watershed_R,
                         iInterpModel,
                         na.rm = TRUE)[["var1.pred"]]
  
  # Change name back
  names(iInterp) <- i
  
  # Return
  return(iInterp)
  
}) %>% rast() # Join all layers together into one spatRast

# Remove objects no longer needed
rm(chemData, watershed_R)
gc()

### REMOVE EXTREMELY LOW ESTIMATE APPLICATIONS ---------------------------------
# Interpolation will overestimate the presence of extremely low applications
# which are more likely to be genuine absence, i.e. in National Parks

# Conservatively reassign values < 0.001 (less than 1 gram/km^2/year) to NA
chemDataInterp <- clamp(chemDataInterp,
                        lower = 0.001, # 1 gram/km^2/year
                        values  = FALSE) # Set NAs below value

### EXTRACT DATA TO WATERSHEDS -------------------------------------------------
# N.B. Warning: this runs overnight

# Weighted sum (missing data is treated as 0)
# N.B. Since each 1x1km data square value has units of kg/year,
# a weighted sum extract function for each watershed results in estimated 
# kg/year within that watershed
watershedChemData <- terra::extract(chemDataInterp, watershedData,
                                    exact = TRUE,
                                    fun = sum, na.rm = TRUE,
                                    ID = TRUE) # N.B. Same as watershedData ID

# Save
saveRDS(watershedChemData,
        file = paste0(dataDir,
                      "Processed/Watersheds/Watershed_chem_data.Rds"))
