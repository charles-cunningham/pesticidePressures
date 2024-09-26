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
lapply(paste0(dataDir, "/Processed/Catchments"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

### FILTER AND SEPARATE FLOW AND CATCHMENT DATA --------------------------------

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
        file = paste0(dataDir, "/Processed/Catchments/Flow_data_only.Rds"))

# Remove objects and clear memory
rm(flowData)
gc()

# CATCHMENT DATA

# Read in catchment data from geopackage, with query that selects:
# columns needed, only River water catchment type, and only segments with 1km
# flow accumulation or greater
catchmentData <- read_sf(dsn = paste0(dataDir, "Raw/Flow_data/Flow_data.gpkg"),
                         query = "
                         SELECT *
                         FROM ea_detailed_watersheds
                         WHERE watercat = 'River'
                         AND country = 'England'
                         ")

# Drop columns not needed
catchmentData <- catchmentData %>%
  select (!c(country, watercat, shape_length, shape_area))

# Save
saveRDS(catchmentData,
        file = paste0(dataDir,
                      "/Processed/Catchments/Catchment_data_only.Rds"))

# Remove objects and clear memory
rm(catchmentData)
gc()

# EXTRACT FERTILISER AND PESTICIDE DATA TO CATCHMENTS --------------------------  

### READ IN CATCHMENT DATA

# Read catchment .Rds
catchmentData <- readRDS(paste0(dataDir,
                                "/Processed/Catchments/Catchment_data_only.Rds"))

### READ IN FERTILSIER DATA

# List files (each contains two layers - (1) data layer, (2) uncertainty layer)
fertFiles <- list.files(paste0(dataDir, "Raw/Fertiliser_data/data"),
                        full.names = TRUE)

# Extract the data layer from each file
fertData <- lapply(fertFiles, function(x) { 
  
  # Create spatRaster from first file...
  rast(x, lyrs = 1)
  
  # ... and combine together
}) %>% rast

# Change layer names to remove "uncertainty"
names(fertData) <- gsub("_prediction_uncertainty_1", "", names(fertData) )

### CALCULATE PER CATCHMENT FERTILISER APPLICATION

#### filter for testing!!!!!!!!!!!!!!!!!
catchmentData <- catchmentData %>%
  filter(wbname == "Heacham River" )
###!!!!!!

# Weighted sum (missing data is treated as 0)
# N.B. Since each 1x1km data square value is estimated amount applied per 1x1km,
# a weighted sum extract function over each catchment results in estimated 
# amount applied within that catchment
catchmentFert <- extract(fertData, catchmentData, exact = TRUE, fun = sum,
                         na.rm = TRUE, ID = FALSE, bind = TRUE)

# Save
saveRDS(catchmentFert,
        file = paste0(dataDir,
                      "/Processed/Catchments/Catchment_fertiliser.Rds"))



# FILTER WATERCOURSES

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

#########

# Find network (inlcude 10m buffer to not miss any spatial errors)
riverNetwork <- st_touches(testData, remove_self = TRUE, dist = 10)

for(i in 1:NROW(testData)) {

iSegments <- iNetwork <- i

isHigher = TRUE

while (is.numeric(iSegments)) {
  
  newSegments <- NULL

  for (j in iSegments) {
 
    iSegmentsTouches <- riverNetwork[j][[1]]
    
    isHigher <- (testData$startz[iSegmentsTouches] >= testData$startz[j]-10) |
      (testData$endz[iSegmentsTouches] >= testData$endz[j]-10)
                   
    if(any(isHigher)) {

      higherSegments <- iSegmentsTouches[isHigher]
      
      newSegments <- setdiff(higherSegments, iNetwork)
      
      iNetwork <- c(iNetwork, higherSegments) %>% unique
      
    }
    
  }
  ## assign iSegments here
  iSegments <- newSegments
  
}

print(iNetwork)
#test
testData[i, "test"] <- any(testData$startz[iNetwork] > 80 )

}


########
for(i in 1:NROW(testData)) {
  
  iSegments <- iNetwork <- i
  
  higherSegments <- which(testData$startz >= testData$startz[i])
  
  higherNetwork <- riverNetwork[higherSegments]
  names(higherNetwork) <- higherSegments
  
  while (is.numeric(iSegments)) {
    
    countSegments <- NULL
    
    for (j in iSegments) {
 
      iSegmentsTouches <- higherNetwork[as.character(j)][[1]]
      
      newSegments <- setdiff(iSegmentsTouches, iNetwork)
      
      iNetwork <- c(iNetwork, newSegments)
      
      countSegments <- c(countSegments, newSegments)
      
      }
      
    ## assign iSegments here
    iSegments <- setdiff(countSegments, iSegments)
    
  }
  
  print(iNetwork)
  #test
  testData[i, "test"] <- any(testData$startz[iNetwork] > 80 )
  
}





plot(testData[, "test"])


ggplot(st_as_sf(testData[iNetwork,])) +
  geom_sf(aes(colour = test)) +
  theme_void()

saveRDS(testData, "test.Rds")
