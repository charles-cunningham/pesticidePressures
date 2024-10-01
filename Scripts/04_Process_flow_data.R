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

# Create processed data folder
lapply(paste0(dataDir, "/Processed/Flow"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

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

# FILTER WATERCOURSES ----------------------------------------------------------

# Read in flow data
flowData <- readRDS(file = paste0(dataDir, "/Processed/Catchments/Flow_data_only.Rds"))

# Filter river basic districts entirely in England
# (no need to do complex river network filtering on these to speed up)
flowData$rbdAllEngland <- flowData$rbd %in% c("Anglian",
                                              "Humber",
                                              "North West",
                                              "South East",
                                              "South West",
                                              "Thames")


flowData[flowData$rbdAllEngland == FALSE,]


ggplot(st_as_sf(test)) +
  geom_sf(aes(fill = test$fertiliser_k)) +
  theme_void()



# extract fertiliser data onto catchments as a test
# crop to England larger catchments only, i.e. no part of a river should come from wales or scotland
# for each river segment, extract catchment data




testData <- flowData[flowData$ea_wb_id == "GB105033043250",]

ggplot(st_as_sf(testData)) +
  geom_sf(aes(colour = endz)) +
  theme_void()


# Find network
riverNetwork <- st_touches(testData)

for(i in 1:NROW(testData)) {
  
  iSegments <- iNetwork <- i
  
  higherSegments <- which(testData$maxflowacc <= testData$maxflowacc[i] |
                            testData$startz <= testData$endz[i])
  
  higherNetwork <- riverNetwork[higherSegments]
  names(higherNetwork) <- higherSegments
  
  while (is.numeric(iSegments)) {

    allNewSegments <- NULL
    
    for (j in iSegments) {

      jSegmentsTouches <- higherNetwork[as.character(j)][[1]]
      
      jNewSegments <- setdiff(jSegmentsTouches, iNetwork)
      
      allNewSegments <- c(allNewSegments, jNewSegments)
      
      }
      
    iNetwork <- c(iNetwork, allNewSegments)
    
    ## assign iSegments here
    iSegments <- setdiff(allNewSegments, iSegments)
    
  }
  
  testData[i, "test"] <- any(testData$startz[iNetwork] > 70 )
  
}



plot(testData[, "test"])


ggplot(st_as_sf(testData[iNetwork,])) +
  geom_sf(aes(colour = test)) +
  theme_void()

saveRDS(testData, "test.Rds")

