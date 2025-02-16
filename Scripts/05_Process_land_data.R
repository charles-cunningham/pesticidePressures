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

### READ IN LAND COVER DATA

# Read land cover data as spatRast file (first layer is land cover class)
lcm2015 <- paste0(dataDir, "Raw/Land_cover_data/lcm2015gb25m.tif") %>%
  rast() %>%
  .[[1]]

# Rename lcm2015 spatRast layer
names(lcm2015) <- "Identifier"

# Specify the LCM classes, add on "No Data" for NA values
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

### READ IN WATERSHED DATA

# Read watershed .Rds
watershedData <- readRDS(paste0(dataDir,
                                "Processed/Watersheds/Watershed_data.Rds"))

# Create data frame from watershedData
# N.B Assigning values later using this tibble speeds up significantly
watershed_tibble <-  matrix(ncol = length(classLCM),
                            nrow = NROW(watershedData)) %>%
  data.frame() %>%
  setNames(., classLCM)


# Find columns that match to LCM classes
colNumsLCM <- names(watershed_tibble) %in% classLCM %>%
  which(.)

# Add total area column (in 25x25m cells)
watershed_tibble[, "totalArea"] <- NA

# Add ID column
watershed_tibble[, "ID"] <- 1:NROW(watershed_tibble)

# EXTRACT COVERAGE (IN 25x25M CELLS) -------------------------------------------

# Create progress bar
progressBar = txtProgressBar(
  min = 0,
  max = NROW(watershedData),
  initial = 0,
  style = 3
)

# Start loop iterating through every watershed_tibble row
for (i in 1:NROW(watershed_tibble)) { # (Same row numbers as watershedData)

  # Extract all 25x25m cells for each land cover class present for watershed i
  # N.B. This is how rows are connected
  watershedCells <- terra::extract(lcm2015, watershedData[i,])
  
  # Count number of cells for each class
  # N.B. some classes may not be included as count is 0
  watershedCount <- count(watershedCells, Identifier, name = "Cover")

  # Add class names by joining coverage values to LCM data frame
  watershedCount <- full_join(LCM_df, watershedCount, by = "Identifier") %>%
    replace_na(list(Cover = 0)) # Convert 'Cover' NAs to 0
  
  # Add coverage for each class to watershed_tibble (row i)
  watershed_tibble[i, colNumsLCM] <- watershedCount$Cover
  
  # Add total cover
  watershed_tibble[i, "totalArea"] <- sum(watershedCount$Cover)

  # Iterate progress bar
  setTxtProgressBar(progressBar, i)
  
}

# Close progress bar
close(progressBar)

# Save
saveRDS(watershed_tibble,
        file = paste0(dataDir,
                      "Processed/Watersheds/Watershed_land_data.Rds"))
