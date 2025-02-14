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
flowChemData <- readRDS(paste0(dataDir,
                               "Processed/Flow/Flow_chem_data.Rds"))

# Read pesticide data
pestData <- read.csv(paste0(dataDir,
                            "Raw/Pesticide_data/Pesticide_properties.csv"))

### CREATE TOXICITY METRICS -----------------------------------------

# SUM
# Total upstream application of all pesticides in kg

# Create pesticideLoad column as sum of all pesticide application
flowChemData <- flowChemData %>%
  mutate(pesticideLoad =
           # If any pesticide layers are not NA ...
           # (needed otherwise all NA rows become 0, even with na.rm = TRUE)
           case_when(if_any(starts_with('pesticide_'),
                            complete.cases)
                     # Sum all pesticide values
                     ~ rowSums(across(starts_with('pesticide_')))))

# CHECK ORDER FOR WEIGHTED SUM
# (flowData columns have to be in same order as pestData rows to match up)

# pestData; Create simple pesticide name for matching
pestData$SimpleName <- pestData$Pesticide %>%
  gsub("\\(", "", .) %>%
  gsub("\\)", "", .) %>%
  gsub("\\[", "", .) %>%
  gsub("\\]", "", .) %>%
  gsub("\\-", "", .) %>%
  gsub("\\,", "", .) %>%
  gsub("\\.", "", .) %>%
  gsub("\\/", "", .) %>%
  gsub("\\&", "", .) %>%
  gsub("\\'", "", .) %>%
  gsub(" ", "", .)

# flowData: Create simple pesticide name for matching
pestData$FlowName <- names(flowChemData) %>%
  .[grepl("pesticide_", .)] %>%
  gsub("pesticide_", "", .) %>%
  gsub("\\.", "", .)

# Check that both match
if (identical(pestData$SimpleName, pestData$FlowName)) {
  print("pestData rows and flowData columns match")
} else {
  stop("pestData rows and flowData columns do not match")
}

# Remove columns as no longer needed
pestData$SimpleName <- pestData$FlowName <- NULL

# CREATE TOXICITY WEIGHTS

# Create toxicity column to be used for weighted sum  (inverse of NOEC)
toxicity <- 1/
  pestData$Toxicity_.Temperate_Freshwater_Aquatic_invertebrates_Chronic_21_day_NOEC_.mgl.1..

# Replace NAs with min toxicity value (max NOEC)
toxicity[is.na(toxicity)] <- min(toxicity,
                                 na.rm = TRUE)

# Set names of toxicity to flowChem pesticide column names
# N.B. previously checked order is correct
names(toxicity) <- grep("pesticide_",
                         names(flowChemData),
                         value = TRUE)

# WEIGHTED SUM

# Create pesticideToxicLoad column as sum of all pesticide application
flowChemData <- flowChemData %>%
  mutate(pesticideToxicLoad =
           # If any pesticide layers are not NA ...
           # (needed otherwise all NA rows become 0, even with na.rm = TRUE)
           case_when(if_any(starts_with('pesticide_'),
                            complete.cases)
                     # Sum all pesticide values
                     ~ rowSums(across(starts_with('pesticide_'),
                                      # Multiply by associated toxicity weights
                                      ~ .x * toxicity[cur_column()]),
                               na.rm = TRUE)))

### CREATE PESTICIDE DIVERSITY METRICS ------------------------------------------

# Create diversity indices
flowChemData <- flowChemData %>%
  
  # Create pesticideDiv column as Shannon diversity index
  mutate(pesticideShannon =
           vegan::diversity(across(starts_with('pesticide_')),
                            index = "shannon")) %>%
  
  # Create pesticideDiv column as Simpson diversity index
  mutate(pesticideSimpson =
           vegan::diversity(across(starts_with('pesticide_')),
                            index = "simpson"))

### SAVE -----------------------------------------------------------------------

# Save processed flow data
saveRDS(flowChemData,
        file = paste0(dataDir,
                      "/Processed/Flow/Flow_aggregated_data.Rds"))
