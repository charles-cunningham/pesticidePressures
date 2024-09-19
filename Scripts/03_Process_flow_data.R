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
lapply(paste0(dataDir, "/Processed"), function(x) {
  if (!file.exists(x)) {
    dir.create(x)
  }
})




test <- read_sf(dsn = paste0(dataDir, "Raw/Flow_data/Flow_data.gpkg"),
                query = 
                "SELECT
                *
                FROM
                Overland_flow_pathways — ea_probable_overland_flow_pathways
                WHERE
                (
                (flowacccl = '1Km') OR (flowacccl = '10Km')
                )
                AND
                (
                (endz >= 1000)
                )")
        
        
  
# read in river map 
# read in catchment map [ watercat = River AND country = England]
# extract fertiliser data onto catchments as a test
# for each river segment, extract catchment data
# crop to England larger catchments only,
# i.e. no part of a river should come from wales or scotland

