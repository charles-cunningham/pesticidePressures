# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process overland flow data
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

# Create processed flow data folder
lapply(paste0(dataDir, "Processed/Flow"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

# Read England spatVector
england <- readRDS(paste0(dataDir,
                          "Raw/Country_data/England.Rds"))

# Convert England to sf object for later processing
england <- england %>%
  st_as_sf %>%
  st_transform(., "EPSG:27700")

# Read catchment fertiliser data
catchmentFert <- readRDS(paste0(dataDir,
                                "/Processed/Catchments/Catchment_fertiliser.Rds"))

# Read catchment pesticide data
# catchmentPest <- readRDS(paste0(dataDir,
#                                 "/Processed/Catchments/Catchment_pesticides.Rds"))

### FILTER AND SEPARATE FLOW DATA ----------------------------------------------

# FLOW DATA

# Read in flow data from geopackage, with query that selects:
# only River water catchment type, 
# and only segments with 1km flow accumulation or greater
flowData <- read_sf(dsn = paste0(dataDir, "Raw/Flow_data/Flow_data.gpkg"),
                    query = "
                SELECT *
                FROM ea_probable_overland_flow_pathways
                WHERE water_cat >= 'River'
                AND ( flowacccl = '1Km' OR flowacccl = '10Km' )
                ")

# Drop columns not needed
flowData <- flowData %>%
  select (!c(water_cat, ponding, zdiff, flowacccl, shape_length))

# Save
saveRDS(flowData,
        file = paste0(dataDir, "/Processed/Flow/Flow_data_only.Rds"))

# Remove objects and clear memory
rm(flowData)
gc()

# SET FLOW AND BASIN DATA ------------------------------------------------------

# Read in flow data
flowData <- readRDS(file = paste0(dataDir, "/Processed/Flow/Flow_data_only.Rds"))

# Set river basin districts
basins <- flowData$rbd %>% unique

# Set basins entirely within England
allEnglishBasins <- c("Anglian",
                      "Humber",
                      "North West",
                      "South East",
                      "South West",
                      "Thames")

# Add column indicating whether upsteam network is entirely within England 
flowData$withinEngland <- if_else(flowData$rbd %in% allEnglishBasins,
                                  "Yes",
                                  NA)

# EXTRACT CATCHMENT DATA TO CATCHMENTS -----------------------------------------

###!!!###
#testing#
###!!!###
testData <- flowData[flowData$opcatch == "Derwent",]
###!!!###
#testing#
###!!!###

### SET UP BASIN LOOP

# Loop through each basic separately to speed up
for(basin in "North West") { # for(i in basins) {
  
  # Subset to river basin
  basinFlow <- testData[testData$rbd == basin,] #~! change
  # basinFlow <- flowData[flowData$rbd == basin,]
  
  # Find all intersections between river segments within basin
  basinNetwork <- st_touches(basinFlow)
  
### SET UP SEGMENT LOOP AND IDENTIFY SMALLER STREAMS
  
  # For each segment in the river basin
  for(i in 1:NROW(basinFlow)) {
  
  # Find all segments with lower or equal to cumulative flow than segment i
  # (including i)
  higherSegments <- which(basinFlow$maxflowacc <= basinFlow$maxflowacc[i])
  
  # Subset the network to segment i and higherSegments, 
  higherNetwork <- basinNetwork[higherSegments]
  
  # Rename segments to correct number (number is lost on subsetting)
  names(higherNetwork) <- higherSegments
  
### ESTABLISH CONNECTED UPSTREAM NETWORK ---------------------------------------
  
  # Set starting current segment as i; set starting network as i
  iSegments <- iNetwork <- i

  # While there are still new segments to check, i.e. iSegments is not NULL ...
  while (is.numeric(iSegments)) {

    # Reset new segments to NULL, this will be filled in later 
    allNewSegments <- NULL
    
    # For every segment j (for first iteration this is 1 but can be >1)
    for (j in iSegments) {

      # Which higher segments is segment j connected to (can be >1)?
      jSegmentsTouches <- higherNetwork[as.character(j)][[1]]
      
      # Which new segments will be added to the network?
      jNewSegments <- setdiff(jSegmentsTouches, iNetwork)
      
      # Running tally of new added segments
      allNewSegments <- c(allNewSegments, jNewSegments)
      
      }
      
    # Add all new segments to the network
    iNetwork <- c(iNetwork, allNewSegments)
    
    # Assign new segments as iSegments for next iteration
    iSegments <- setdiff(allNewSegments, iSegments)
    
  }
  
### CHECK IF IN ENGLAND --------------------------------------------------------
# N.B. The st_interests funtion takes a while for a large stream network,
# hence prior if statements to reduce number of times this is needed
  
  # If entire river basin is not in England (NA), we need to check upstream...
  if (is.na(basinFlow$withinEngland[i])) {
    
    # Check if any upstream segments have already been classified as non-English
    if(any(basinFlow[inetwork, "withinEngland"] == "No")) {
      
      # If so, assign "No"
      basinFlow$withinEngland[i] <- "No"
      
      # Otherwise need to check manually
      } else {
        
        # Do all upstream segments intersect with England polygon?
        basinFlow$withinEngland[i] <- if_else(
          st_intersects(basinFlow[inetwork,],
                        england) %>%
            as.logical %>%
            all,
          "Yes",
          "No")
    }
    
### EXTRACT FLOW VARIABLES FROM CATCHMENTS -------------------------------------
  
  if (basinFlow$withinEngland[i] == TRUE) {
    
    
    extract columns = values
  
    } else {
    
    extract columns = NA
  }
  
  
}

  
  # Run function on entire network above segment
  testData[i, "test"] <- any(testData$startz[iNetwork] > 400 )


  
  

ggplot(st_as_sf(testData)) +
  geom_sf(aes(colour = test)) +
  theme_void()


saveRDS(testData, "test.Rds")

