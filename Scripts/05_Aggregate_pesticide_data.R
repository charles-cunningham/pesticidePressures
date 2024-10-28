# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Aggregate pesticide data
#
# Script Description:

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(sf)

### DIRECTORY MANAGEMENT -------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

### LOAD DATA ------------------------------------------------------------------

# Read processed flow data
# flowChemData <- readRDS(paste0(dataDir,
#                                "/Processed/Flow/Flow_chem_data.Rds"))
flowChemData <- readRDS(paste0(dataDir,
                               "/Processed/Flow/basin_Dee_chem_data.Rds"))

# Find pesticide layers
pestLayers <- grep("pesticide", names(flowChemData), value = TRUE)

### CREATE TOTAL PESTICIDE LOAD METRIC -----------------------------------------
# N.B. This is a placeholder for now. Currently just sum, but need to introduce 
# weighted sum using toxicity

# Create pesticide load column as sum of all pesticide applications
flowChemData <- flowChemData %>% 
  mutate(pesticideLoad = rowSums(across(starts_with('pesticide_')),
         na.rm = TRUE))

### CREATE PESTICIDE DIVERSITY METRIC ------------------------------------------
# N.B. This is a placeholder for now. Currently just total number, but need 
# to introduce diversity metric describing evenness of pesticide applications

# Create pesticide diversity column as count of all pesticide applications
flowChemData <- flowChemData %>% 
  mutate(pesticideDiv = rowSums(across(starts_with('pesticide_'),
                                       ~ .x > 0),
                                na.rm = TRUE))

### SAVE -----------------------------------------------------------------------

#...


