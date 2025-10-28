# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Aggregate upstream data
#
# Script Description: 

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(sf)
library(vegan)

### DIRECTORY MANAGEMENT -------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

### LOAD DATA ------------------------------------------------------------------

# Read processed flow data
flowDataAll <- readRDS(paste0(dataDir,
                               "Processed/Flow/Flow_data_all.Rds"))

# Read pesticide data
pestData <- read.csv(paste0(dataDir,
                            "Raw/Pesticide_data/Pesticide_properties.csv"))

### CREATE PESTICIDE METRICS ---------------------------------------------------

### SUM
# Total upstream application of all pesticides in kg

# Create pesticideLoad column as sum of all pesticide application
flowDataAll <- flowDataAll %>%
  mutate(pesticideLoad =
           # If any pesticide layers are not NA ...
           # (needed otherwise all NA rows become 0, even with na.rm = TRUE)
           case_when(if_any(starts_with('pesticide_'),
                            complete.cases)
                     # Sum all pesticide values
                     ~ rowSums(across(starts_with('pesticide_')))))

### WEIGHTED SUMS

# CHECK ORDER FOR WEIGHTS
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
pestData$FlowName <- names(flowDataAll) %>%
  .[grepl("pesticide_", .)] %>%
  gsub("pesticide_", "", .) %>%
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

# Check that both match
if (identical(pestData$SimpleName, pestData$FlowName)) {
  print("pestData rows and flowData columns match")
} else {
  stop("pestData rows and flowData columns do not match")
}

# Remove columns as no longer needed
pestData$SimpleName <- pestData$FlowName <- NULL

# CREATE TOXICITY WEIGHTS

# Create toxicity vector to be used for weighted sum  (inverse of NOEC)
toxicity <- 1/
  pestData$Toxicity_.Temperate_Freshwater_Aquatic_invertebrates_Chronic_21_day_NOEC_.mgl.1..

# Replace NAs with min toxicity value (max NOEC)
toxicity[is.na(toxicity)] <- min(toxicity,
                                 na.rm = TRUE)

# CREATE LEACHING WEIGHTS

# Create leachability vector
leachability <- pestData$Leachability..GUS_leaching_potential_index.

# Remove values under 0 (extremely low leaching potential) and NA values
leachability[leachability < 0 | is.na(leachability)] <- 0

# CREATE -CIDE WEIGHTS

insecticide <- pestData$Insecticide
herbicide <- pestData$Herbicide
fungicide <- pestData$Fungicide

# SET NAMES

# Set names of toxicity to flowChem pesticide column names
# N.B. previously checked order is correct
names(toxicity) <- names(leachability) <-
  names(insecticide) <- names(herbicide) <- names(fungicide) <-
  grep("pesticide_",
       names(flowDataAll),
       value = TRUE)

# CALCULATE WEIGHTED SUMS

# Create pesticideToxicLoad column as sum of all pesticide application
flowDataAll <- flowDataAll %>%
  mutate(pesticideToxicLoad =
           # If any pesticide layers are not NA ...
           # (needed otherwise all NA rows become 0, even with na.rm = TRUE)
           case_when(if_any(starts_with('pesticide_'),
                            complete.cases)
                     # Sum all pesticide values
                     ~ rowSums(across(starts_with('pesticide_'),
                                      # Multiply by associated toxicity weights
                                      ~ .x * 
                                        toxicity[cur_column()] *
                                        leachability[cur_column()] ),
                               na.rm = TRUE)))

# Create insecticideToxicLoad column as sum of insecticide application
flowDataAll <- flowDataAll %>%
  mutate(insecticideToxicLoad =
           # If any pesticide layers are not NA ...
           # (needed otherwise all NA rows become 0, even with na.rm = TRUE)
           case_when(if_any(starts_with('pesticide_'),
                            complete.cases)
                     # Sum all pesticide values
                     ~ rowSums(across(starts_with('pesticide_'),
                                      # Multiply by associated toxicity weights
                                      ~ .x * 
                                        toxicity[cur_column()] *
                                        leachability[cur_column()] *
                                        insecticide[cur_column()] ),
                               na.rm = TRUE)))

# Create insecticideToxicLoad column as sum of herbicide application
flowDataAll <- flowDataAll %>%
  mutate(herbicideToxicLoad =
           # If any pesticide layers are not NA ...
           # (needed otherwise all NA rows become 0, even with na.rm = TRUE)
           case_when(if_any(starts_with('pesticide_'),
                            complete.cases)
                     # Sum all pesticide values
                     ~ rowSums(across(starts_with('pesticide_'),
                                      # Multiply by associated toxicity weights
                                      ~ .x * 
                                        toxicity[cur_column()] *
                                        leachability[cur_column()] *
                                        herbicide[cur_column()] ),
                               na.rm = TRUE)))

# Create insecticideToxicLoad column as sum of all fungicide application
flowDataAll <- flowDataAll %>%
  mutate(fungicideToxicLoad =
           # If any pesticide layers are not NA ...
           # (needed otherwise all NA rows become 0, even with na.rm = TRUE)
           case_when(if_any(starts_with('pesticide_'),
                            complete.cases)
                     # Sum all pesticide values
                     ~ rowSums(across(starts_with('pesticide_'),
                                      # Multiply by associated toxicity weights
                                      ~ .x * 
                                        toxicity[cur_column()] *
                                        leachability[cur_column()] *
                                        fungicide[cur_column()] ),
                               na.rm = TRUE)))

### CREATE PESTICIDE DIVERSITY METRICS -----------------------------------------

# Create diversity indices
flowDataAll <- flowDataAll %>%
  
  # Create pesticideDiv column as Shannon diversity index
  mutate(pesticideShannon =
           diversity(across(starts_with('pesticide_')),
                            index = "shannon")) %>%
  
  # Create pesticideDiv column as Simpson diversity index
  mutate(pesticideSimpson =
           diversity(across(starts_with('pesticide_')),
                            index = "simpson"))

# REMOVE REDUNDANT PESTICIDE COLUMNS -------------------------------------------

# Select all invividual pesticide columns and remove as not needed (aggregated)
flowDataAll <- flowDataAll %>%
  select (-c(starts_with('pesticide_')))

### SAVE -----------------------------------------------------------------------

# Save processed flow data
saveRDS(flowDataAll,
        file = paste0(dataDir,
                      "/Processed/Flow/Flow_aggregated_data.Rds"))
