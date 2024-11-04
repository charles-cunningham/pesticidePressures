# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process catchment data
#
# Script Description:

### TODO
# - Add other variables to extract to catchment, i.e. land cover

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
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Raw/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

# Create processed catchement data folder
lapply(paste0(dataDir, "Processed/Catchments"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

# Read England spatVector, and convert to sf object for later processing
england <- readRDS(paste0(dataDir,
                          "Raw/Country_data/England.Rds")) %>%
  st_as_sf

### SEPARATE, FILTER AND SAVE CATCHMENT DATA -----------------------------------

# Read in catchment data from geopackage, with query that selects:
# columns needed, only River water catchment type, and only segments with 1km
# flow accumulation or greater
catchmentData <- read_sf(dsn = paste0(dataDir, "Raw/Flow_data/Flow_data.gpkg"),
                         query = "
                         SELECT *
                         FROM ea_detailed_watersheds
                         WHERE watercat = 'River'
                         ")

# Drop columns not needed, and set standardised crs
catchmentData <- catchmentData %>%
  select (!c(country, watercat, shape_length, shape_area)) %>%
  st_set_crs("EPSG:27700")

# Filter catchments to those which intersect with England country boundary
# N.B. Cannot use original data column "country" as unreliable at Welsh border
catchmentData <- catchmentData %>%
  filter(st_intersects(., england, sparse = FALSE)[,1])

# Save
saveRDS(catchmentData,
        file = paste0(dataDir,
                      "Processed/Catchments/Catchment_data_only.Rds"))

# Remove objects and clear memory
rm(catchmentData)
gc()

### READ IN CATCHMENT, FERTILISER AND PESTICIDE DATA ---------------------------  

### READ IN CATCHMENT DATA

# Read catchment .Rds
catchmentData <- readRDS(paste0(dataDir,
                                "Processed/Catchments/Catchment_data_only.Rds"))

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

### JOIN FERTILISER AND PESTICIDE DATA TOGETHER --------------------------------

# Join fertliser and pesticide spatRasters into single combined spatRaster
chemData <- c(fertData, pestData)

# Save layer names
chemNames <- names(chemData)

# Remove separate objects
rm(fertData, pestData)
gc()

### INTERPOLATE RASTER DATA  ---------------------------------------------------

# Rasterise all cells that touch catchments which intersect with England
# N.B. There are some gaps in the Land Cover Plus data to interpolate such as
# catchments along England/Wales border
catchment_R <- vect(catchmentData) %>%
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
  iInterp <- interpolate(catchment_R,
                         iInterpModel,
                         na.rm = TRUE)[["var1.pred"]]
  
  # Change name back
  names(iInterp) <- i
  
  # Return
  return(iInterp)

}) %>% rast() # Join all layers together into one spatRast

# Remove objects no longer needed
rm(chemData, catchment_R)
gc()

### REMOVE EXTREMELY LOW ESTIMATE APPLICATIONS ---------------------------------
# Interpolation will overestimate the presence of extremely low applications
# which are more likely to be genuine absence, i.e. in National Parks

# Conservatively reassign values < 0.001 (less than 1 gram/km^2/year) to NA
chemDataInterp <- clamp(chemDataInterp,
                        0.001, # 1 gram/km^2/year
                        values  = FALSE) # Set NAs below value

### EXTRACT DATA TO CATCHMENTS -------------------------------------------------
# N.B. Warning: this runs overnight

# Weighted sum (missing data is treated as 0)
# N.B. Since each 1x1km data square value has units of kg/year,
# a weighted sum extract function for each catchment results in estimated 
# kg/year within that catchment
catchmentChemData <- terra::extract(chemDataInterp, catchmentData,
                                    exact = TRUE,
                                    fun = sum, na.rm = TRUE,
                                    ID = FALSE, bind = TRUE)

# Save
saveRDS(catchmentChemData,
        file = paste0(dataDir,
                      "Processed/Catchments/Catchment_chem_data.Rds"))
