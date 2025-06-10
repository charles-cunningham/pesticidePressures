# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Join Biosys and site data
#
# Script Description: 

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)

### DIRECTORY MANAGEMENT -------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

### LOAD DATA ------------------------------------------------------------------

# BIOSYS

# Read Biosys data
invData_sf <- readRDS(file = paste0(dataDir,
                                    "/Processed/Biosys/invDataSpatial.Rds"))

# SITE DATA

# Collect relevent site variable file names
siteFiles <- list.files(paste0(dataDir, "/Raw/Site_data/data"), full.names = TRUE) %>%
  grep("sitevariables.csv", ., value = TRUE)

# Read and bind files
siteData <- lapply(siteFiles, read.csv) %>% 
  bind_rows()

### JOIN DATA ------------------------------------------------------------------

invData_sf <- left_join(invData_sf, siteData,
                  by = c("SITE_ID" = "BIO_SITE_ID"))

### SAVE -----------------------------------------------------------------------

# Save processed invData ready for modelling
saveRDS(invData_sf,
        file = paste0(dataDir,
                      "/Processed/Biosys/invDataSpatial.Rds"))
