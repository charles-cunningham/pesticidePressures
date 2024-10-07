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

### DIRECTORY MANAGEMENT -------------------------------------------------------

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

### LOAD DATA ------------------------------------------------------------------

# Read England spatVector
england <- readRDS(paste0(dataDir,
                          "Raw/Country_data/England.Rds"))

# Convert England to sf object for later processing
england <- england %>%
  st_as_sf 

# Read catchment fertiliser data
catchmentFert <- readRDS(paste0(dataDir,
                                "/Processed/Catchments/Catchment_fertiliser.Rds")) %>%
  st_as_sf

# List fertiliser layer names
fertLayers <- grep("fertiliser", names(catchmentFert), value = TRUE)

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

# Drop columns not needed, and set standardised crs
flowData <- flowData %>%
  select (!c(water_cat, ponding, zdiff, flowacccl, shape_length)) %>%
  st_set_crs("EPSG:27700")

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

# Important - sort data from smallest to largest flow to speed up later steps
flowData <- arrange(flowData, maxflowacc)

# EXTRACT CATCHMENT DATA TO RIVER SEGMENTS -------------------------------------

###!!!###
#testing#
###!!!###
#testData <- flowData[flowData$rbd == "Solway Tweed",]
#testData <- flowData[flowData$opcatch == "Esk and Irthing",]  
testData <- flowData[flowData$ea_wb_id == "GB102077074190",]
###!!!###
#testing#
###!!!###

# start system time
system.time(

### SET UP BASIN LOOP

# Loop through each basin separately to speed up
for(basin in "Solway Tweed") { # for(i in basins) {
  
  # Subset to river basin
  basinFlow <- testData[testData$rbd == basin,] #~! change
  # basinFlow <- flowData[flowData$rbd == basin,]
  
  # Find all intersections between river segments within basin
  segmentNetwork <- st_touches(basinFlow)
  
  # Set up list of upstream networks for each segment
  upstreamNetwork <- vector(mode = "list",
                            length = NROW(basinFlow))
  
  # Filter catchment data to basin only
  basinFert <- catchmentFert %>%
    filter(rbd == basin)
  
  # Add fertiliser layers
  basinFlow[, fertLayers] <- NA
  
### SET UP SEGMENT LOOP WITH IMMEDIATE NETWORK CHECK
  
  # For each segment i in the river basin
  for(i in 1:NROW(basinFlow)) {
    
    # Set starting network as i
    iNetwork <- i
 
    # Which segments is segment i connected to (can be >1)?
    checkSegments <- segmentNetwork[i][[1]]
    
    # Filter segments that are upstream of i (important initial check)
    checkSegments <- checkSegments[basinFlow$maxflowacc[checkSegments] <= 
                                           basinFlow$maxflowacc[i]]

### ESTABLISH ENTIRE CONNECTED UPSTREAM NETWORK

    # While there are still new segments to check
    while (!is.null(checkSegments)) {
      
      # Add segments to check to network
      iNetwork <- c(iNetwork, checkSegments)

      # Reset allNewSegments (new segments that are added)
      allNewSegments <- c()
      
      # For every j segment to check (one value for first iteration but can be >1)
      for (j in checkSegments) {
        
        # Which segments is segment j connected to (can be >1)?
        jSegmentsTouches <- segmentNetwork[j][[1]]
        
        # Make sure new segments are not already in network (so unidirectional)
        newSegments <- setdiff(jSegmentsTouches, iNetwork)
        
        # Running tally of all new added segments
        allNewSegments <- c(allNewSegments, newSegments)
        
      }
      
      # Make sure all new segments are unique
      allNewSegments <- unique(allNewSegments)
 
      # Reset checkSegments (segments to check on next iteration)
      checkSegments <- c()
      
      # For every new segment...
      for (k in allNewSegments) {
        
        # If upstream network for k has not been filled in yet...
        if (is.null(upstreamNetwork[[k]])) {

          # Add k to list of segments to check for next iteration
          checkSegments <- c(checkSegments, k)
          
          } else {
            
            # Add that upstream network to the existing network
            iNetwork <- c(iNetwork, upstreamNetwork[[k]]) %>%
              unique
          }
      }
    }

    # Save iNetwork in upstreamNetwork list
    upstreamNetwork[[i]] <- iNetwork

### CHECK IF WITHIN ENGLAND
# N.B. The st_interests function takes a while for a large stream network,
# hence prior if statements to reduce number of times this is needed

  # If entire river basin is not in England (NA), we need to check upstream...
  if (is.na(basinFlow$withinEngland[i])) {

    # Check if any upstream segments have already been classified as non-English
    if(any(basinFlow[iNetwork,]$withinEngland == "No", na.rm = TRUE)) {

      # If so, assign "No"
      basinFlow$withinEngland[i] <- "No"

      # Otherwise need to check manually
      } else {

        # Count number of iNetwork segments contained within England
        englishSegmentsN <- st_contains(england,
                                        basinFlow[iNetwork,],)[[1]] %>%
          length

        # Does number in englishSegmentsN match iNetwork (i.e. all English)
        basinFlow$withinEngland[i] <- if_else(
          englishSegmentsN == length(iNetwork),
          "Yes",
          "No")
      }
  }
    
### EXTRACT FLOW VARIABLES FROM CATCHMENTS

    # # If upstream of segment is entirely within England
    # if (basinFlow$withinEngland[i] == "Yes") {
    #   
    #   # Extract all catchments that intersect with upstream network
    #   iCatchment <- basinFert %>%
    #     filter(lengths(st_intersects(.,
    #                                  basinFlow[iNetwork,])) > 0) %>%
    #     as.tibble # Convert to tibble
    #   
    #   # For every fertiliser layer
    #   for (x in fertLayers) {
    # 
    #     # Sum all catchment values
    #     basinFlow[i, x] <- iCatchment[, x] %>%
    #       sum(., na.rm = TRUE)
    # 
    #   }
    # }
    
  }
}
)

ggplot(st_as_sf(basinFlow)) +
  geom_sf(aes(colour = withinEngland)) + 
  theme_void()

saveRDS(basinFlow, "test.Rds")
