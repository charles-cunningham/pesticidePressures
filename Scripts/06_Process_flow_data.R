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
library(parallel)

### CHECK CORES ----------------------------------------------------------------

# Check cores
nCores <- detectCores()

# Assign cores to most efficient number (speed vs memory)
options("mc.cores" = 4)

### DIRECTORY MANAGEMENT -------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

# Create processed flow data folder
lapply(paste0(dataDir, "Processed/Flow/Basins"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

### LOAD DATA ------------------------------------------------------------------

# Read England spatVector, and convert to sf object for later processing
england <- readRDS(paste0(dataDir,
                          "Raw/Country_data/England.Rds")) %>%
  st_as_sf

# Read watershed spatial data
watersheds <- readRDS(paste0(dataDir,
                                "Processed/Watersheds/Watershed_data.Rds")) %>%
  st_as_sf

# Read watershed chemical application data
watershedChemData <- readRDS(paste0(dataDir,
                                    "Processed/Watersheds/Watershed_chem_data.Rds"))

# Read watershed land cover data
watershedLandData <- readRDS(paste0(dataDir,
                                    "Processed/Watersheds/Watershed_land_data.Rds"))

# List chemical application and land cover layer names
chemLayers <- select(watershedChemData, -ID) %>%
  names()
landLayers <- select(watershedLandData, -ID) %>%
  names()

# Join watershed data together (drop ID column from watershedLandData as duplicate)
watershedData <- cbind(watershedChemData,
                       select(watershedLandData, -ID))

### FILTER AND SEPARATE FLOW DATA ----------------------------------------------

# Read in flow data from geopackage, with query that selects:
# - only River water catchment type,
# - only segments with 1km flow accumulation or greater,
# - only segments below the maximum value of flow accumulation (as this value is
# capped it causes values in the upstream search algorithm)
flowData <- read_sf(dsn = paste0(dataDir, "Raw/Flow_data/Flow_data.gpkg"),
                    query = "
                SELECT *
                FROM ea_probable_overland_flow_pathways
                WHERE water_cat = 'River'
                AND ( flowacccl = '1Km' OR flowacccl = '10Km' )
                AND maxflowacc < 2147483646
                ")

# Drop columns not needed, and set standardised crs
flowData <- flowData %>%
  select (!c(water_cat, ponding, zdiff, flowacccl, shape_length)) %>%
  st_set_crs("EPSG:27700")

# Save
saveRDS(flowData,
        file = paste0(dataDir, "Processed/Flow/Flow_data_only.Rds"))

# Remove objects and clear memory
rm(flowData)
gc()

### SET FLOW AND BASIN DATA ----------------------------------------------------

# Read in flow data
flowData <- readRDS(file = paste0(dataDir, "Processed/Flow/Flow_data_only.Rds"))

# Set river basin districts
# N.B. If using 4 cores, quicker to move largest basin (Anglian) to position 3
basins <- flowData$rbd %>% 
  unique %>%
  fct_relevel("Anglian", after = 2) %>%
  levels

# Set basins entirely within England (saves time later)
allEnglishBasins <- c("Anglian",
                      "Humber",
                      "North West",
                      "South East",
                      "South West",
                      "Thames")

# Add column indicating whether upstream network is entirely within England 
flowData$withinEngland <- if_else(flowData$rbd %in% allEnglishBasins,
                                  "Yes",
                                  NA)

# Important - sort data from smallest to largest flow to speed up later steps
flowData <- arrange(flowData, maxflowacc)

### EXTRACT WATERSHED DATA TO RIVER SEGMENTS (BASIN LOOP) ----------------------

### SET UP BASIN LOOP

# Apply function to each basin separately on separate cores to speed up
mclapply(basins, function(basin) {
  
  # Subset to river basin
  basinFlow <- flowData[flowData$rbd == basin,]
  
  # Add chemical application and land cover layers
  basinFlow[, c(chemLayers, landLayers)] <- NA
  
  # Filter watershed polygons to basin only
  basinWatersheds <- watersheds %>%
    filter(rbd == basin)
  
  # Set up list of upstream networks for each segment
  upstreamNetwork <- vector(mode = "list",
                            length = NROW(basinFlow))
  
  ### FIND ALL SEGMENT INTERSECTIONS WITHIN BASIN
  # N.B. Use st_is_within_distance rather than st_touches for any misalignment
  # errors in original flow data. This takes much longer and so a grid 
  # iteration workflow is needed to speed up
  
  # Create 20x20 hexagonal grid, and (2 x 20m search distance) buffer
  basinGrid <- st_make_grid(basinWatersheds, n=c(20,20), square = F) %>%
    st_buffer(., 40)
  
  # Intersect buffered grid with watersheds
  basinWatershedGrid <- basinWatersheds %>%
    st_join(basinGrid %>%
              st_as_sf() %>%
              tibble::rowid_to_column('Cell'),
            join=st_intersects, left = FALSE)
  
  # Intersect buffered grid with rivers
  basinFlowGrid <- basinFlow %>%
    st_join(basinGrid %>%
              st_as_sf() %>%
              tibble::rowid_to_column('Cell'),
            join=st_intersects, left = FALSE)
  
  # Use data.table to iterate st_is_within_distance over each grid cell
  segmentNetwork <- basinFlowGrid %>%
    group_nest(Cell, keep=FALSE) %>%
    mutate(segmentNeighbours = map(
      data,
      ~st_is_within_distance(.x,
                             remove_self = TRUE,
                             dist = 20)) # Search distance of 20m
    )
  
  ### SET UP SEGMENT LOOP
  
  # For each segment i in the river basin
  for(i in 1:NROW(basinFlow)) {
    
    # Set empty starting network
    iNetwork <- c()
    
    # Reset segments to check as i
    checkSegments <- i
    
    ### ESTABLISH ENTIRE CONNECTED UPSTREAM NETWORK
    
    # While there are still new segments to check
    while (!is.null(checkSegments)) {
      
      # Add segments to check to network
      iNetwork <- c(iNetwork, checkSegments)
      
      # Reset allNewSegments (all new segments that are added)
      allNewSegments <- c()
      
      # For every j segment to check (1 value for 1st iteration but can be >1)
      for (j in checkSegments) {
        
        # Reset jTouches (new segments that are added for each j)
        jTouches <- c()
        
        # Extract grid row number from original row number
        # N.B. these are not the same as some duplicates created by buffer
        jRows <- which(basinFlowGrid$permid == basinFlow$permid[j])
        
        # Loop through every matching row (in case >1 in buffer areas)
        for(jRow in jRows) {
          
          # Find cell number
          jCellNumber <- basinFlowGrid$Cell[jRow]
          
          # Find all row numbers (these are rows that st_is_within uses)
          jCellRows <- which(basinFlowGrid$Cell == jCellNumber)
          
          # Find relative position of focal row in jCellRows
          jDistanceNumber <- which(jCellRows == jRow)
          
          # Extract connected segments
          jConnectedNumbers <- segmentNetwork[segmentNetwork$Cell == jCellNumber,
                                              "segmentNeighbours"][[1]][[1]][[jDistanceNumber]]
          
          # Relate back to basinFlow
          basinFlowRows <- which(basinFlow$permid %in%
                                   basinFlowGrid$permid[jCellRows[jConnectedNumbers]])
          
          # Add to all segments to check
          jTouches <- c(jTouches, basinFlowRows)
          
        }
        
        # Make sure new segments are upstream of segment j -
        #  max flow accumulation is "<=" than segment j
        newSegments <-
          jTouches[basinFlow$maxflowacc[jTouches] <= basinFlow$maxflowacc[j] ]
        
        # Running tally of all new added segments
        allNewSegments <- c(allNewSegments, newSegments)
        
      }
      
      # Make sure new segments are unique, and not already in network
      allNewSegments <- unique(allNewSegments) %>%
        setdiff(., iNetwork)
      
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
    
    # If entire river basin is not in England (NA), we need to check upstream...
    if (is.na(basinFlow$withinEngland[i])) {
      
      # Check if any upstream segments have already been classified as non-English
      if(any(basinFlow[iNetwork,]$withinEngland == "No", na.rm = TRUE)) {
        
        # If so, assign "No"
        basinFlow$withinEngland[i] <- "No"
        
        # Otherwise, need to check manually
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
    
    ### EXTRACT FLOW VARIABLES FROM WATERSHEDS

    # If:
    # - upstream of segment is entirely within England,
    # - segment is above the highest astronomical tide in England (this is 8.4m
    #   above Ordnance Datum Newlyn [sea level] from correspondence with UK 
    #   Hydrographic Office) to remove high uncertainty near sea level segments
    if (basinFlow$withinEngland[i] == "Yes" &
        basinFlow$endz[i] > 8.4) {
      
      # Identify basinFlowGrid rows within iNetwork so that we can iterate cells
      basinGridNetwork <- basinFlowGrid %>%
        filter(permid %in% basinFlow[iNetwork,]$permid)
      
      # Iterate through all cells that iNetwork occupies
      cellWatersheds <- lapply(unique(basinGridNetwork$Cell), function(x) {
        
        # Filter basinWatershedGrid to cell x 
        basinWatershedGridCell <- basinWatershedGrid[basinWatershedGrid$Cell == x,]
        
        # Intersect all watersheds within cell x and iNetwork within cell x
        overlaps <- st_intersects(basinWatershedGridCell,
                                  basinGridNetwork[basinGridNetwork$Cell == x,]) %>%
          lengths # Convert to number of intersections per watershed 
        
        # Return basinChemGridCell rows with any overlaps
        return(basinWatershedGridCell[overlaps > 0,])
        
      }) %>% bind_rows # Join together (duplicates from overlapping cells included)
      
      # Use cellWatersheds row IDs to identify unique basinWatersheds
      iWatershedData <- watershedData[watershedData$ID %in% 
                                            cellWatersheds$ID,]

      # For every chemical application  and land cover layer 
      for (layer in c(chemLayers, landLayers)) {
        # Assign upstream application per area for segment i
        basinFlow[i, layer] <-
          iWatershedData[, layer] %>% # Collate upstream watershed layer values
          sum(., na.rm = TRUE) # Sum
      }
    }

    ### PRINT UPDATE
    
    # Every 1000 segments
    if (i %% 1000 == 0) {
      
      # Calculate percentage complete
      percentComplete <- ( i / NROW(basinFlow) ) %>%
        `*` (100) %>%
        round(., digits = 1) # Round to 1 decimal place
      
      # Print
      system(sprintf('echo "\n%s\n"', 
                     paste0(basin, " basin ", percentComplete, "% complete (",
                            i, " of ",  NROW(basinFlow), ")")))
    }
    
    ### END SEGMENT LOOP, SAVE, AND END BASIN LOOP
    
    # End segment loop
  }
  
  # Save basin
  saveRDS(basinFlow,
          file = paste0(dataDir,
                        "Processed/Flow/Basins/",
                        gsub(" ", "_", basin),
                        "_chem_data.Rds"))

  # Remove objects not needed
  rm(basinFlow, segmentNetwork, upstreamNetwork, basinWatersheds)
  gc()
  
})

### COMBINE BASINS INTO SINGLE FILE---------------------------------------------

# Create empty list
basinFlowData <- vector(mode = "list",
                        length = NROW(basins))

# Loop through basins
for (i in 1:length(basins)) {
  
  # Assign each processed basin object to single list 
  basinFlowData[[i]] <- readRDS(file = paste0(dataDir,
                                              "Processed/Flow/Basins/",
                                              gsub(" ", "_", basins[i]),
                                              "_chem_data.Rds"))
}

# Bind list of basin objects together
flowDataAll <- do.call(what = rbind,
                        args = basinFlowData)

# Save combined object
saveRDS(flowDataAll, 
        file = paste0(dataDir,
                      "Processed/Flow/Flow_data_all.Rds"))
