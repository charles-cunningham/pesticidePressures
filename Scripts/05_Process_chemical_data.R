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

### READ IN WATERSHED, LOAD, FERTILISER AND PESTICIDE DATA ---------------------  

### READ IN WATERSHED AND LOAD DATA

# Read watershed .Rds
watershedData <- readRDS(paste0(dataDir,
                                "Processed/Watersheds/Watershed_data.Rds"))

# Read in load data
exportLoad <- paste0(dataDir, "Raw/Pesticide_data/chemical_export.tif") %>%
  rast()

### READ IN LAND COVER DATA

# Read land cover data as spatRast file (first layer is land cover class)
lcm2015 <- paste0(dataDir, "Raw/Land_cover_data/lcm2015gb25m.tif") %>%
  rast(., lyr = 1)

# Rename lcm2015 spatRast layer
names(lcm2015) <- "Identifier"

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

### INTERPOLATE 1KM APPLICATION DATA  ------------------------------------------
# Use gstat inverse weighting approach

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

### INTERPOLATE 100 M EXPORT LOAD ----------------------------------------------
# Use while loop with terra focal approach as gstat too slow due to resolution.
# This approach iteratively fills any data gaps by slowly expanding the focal 
# window

# Disaggregate watershed_R to 100 m resolution
watershed_R_100 <- disagg(watershed_R, fact = 10)

# Change exportload extent to match other datasets
exportLoad <- extend(exportLoad, watershed_R_100)
exportLoad <- crop(exportLoad, watershed_R_100)

# Set up while loop
w <- 1 # window size set to 1 initially (will expand)
filled <- exportLoad # Assign exportLoad to object to be iteratively filled
to_fill <- TRUE # Set up initial logical value to get while loop started

# Start while loop
while(to_fill) {
  
  # Expand window by 2 (has to be odd number)
  w <- w + 2
  
  # Carry out interpolation for all data gap cells within w of a non-NA value
  filled <- focal(filled, w = w, fun = mean, na.policy = "only", na.rm = T)
  
  # Check to see if any data gaps remain - find if any original NA values remain
  remainingGaps <- mask(watershed_R_100, filled, inverse = TRUE)
  to_fill <- any(!is.na(values(remainingGaps)))
}

# Mask to watershed_R_100
exportLoadInterp <- mask(filled, watershed_R_100)

# Remove objects no longer needed
rm(chemData, filled, remainingGaps)
gc()

# CALCULATE LOAD FOR EACH CHEMICAL ---------------------------------------------

# Change land cover map extent to match other datasets
lcm2015 <- extend(lcm2015, watershed_R_100)
lcm2015 <- crop(lcm2015, watershed_R_100)

# Create single value arable spatRast
arable <- classify(lcm2015,
                   cbind(3, 1),
                   others = NA)

# Aggregate arable spatRast to 1km
arable_1km <- aggregate(arable, fact = 40, 
                        sum, na.rm = TRUE)

# Get proportion of arable land in each 100m
arable_prop <- aggregate(arable, fact = 4, sum, na.rm = TRUE) / 
  disagg(arable_1km, fact = 10)

# Multiply proportion by chemDataInterp to downscale to 100m
chemDataInterp_100m <- arable_prop * disagg(chemDataInterp, fact = 10) 

# Multiply chemDataInterp_100m with exportLoad to get total estimated export 
# for each chemical
chemLoad <-  exportLoadInterp * chemDataInterp_100m

# Remove objects no longer needed
rm(chemDataInterp, chemDataInterp_100m,
   exportLoad, exportLoadInterp,
   lcm2015, arable, arable_1km, arable_prop,
   watershed_R, watershed_R_100)
gc()

### EXTRACT DATA TO WATERSHEDS -------------------------------------------------
# N.B. Warning: this runs overnight

# Weighted sum (missing data is treated as 0) of load
# N.B. Since each 1x1km chemDataInterp value has units of kg/year,
# a weighted sum extract function for each watershed  of chemLoad results in
# estimated kg/year applied within that watershed that ends up being exported 
watershedChemData <- terra::extract(chemLoad, watershedData,
                                    exact = TRUE,
                                    fun = sum,
                                    na.rm = TRUE,
                                    ID = TRUE) # N.B. Same as watershedData ID

# Save
saveRDS(watershedChemData,
        file = paste0(dataDir,
                      "Processed/Watersheds/Watershed_chem_data.Rds"))
