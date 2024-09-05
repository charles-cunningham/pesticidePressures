# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Read in Biosys data from online repository
#
# Script Description: Using the Biosys website, this script reads in the separate
# files that comprise the Biosys dataset. It unzips them, and combines them into
# a single file, removing duplicated columns.

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)

### DATA MANAGEMENT ------------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Raw/Biosys/"
# If working locally: "../Data/Raw/Biosys/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Raw/Biosys/"

# Create directory if it doesn't exist
if (!file.exists(dataDir)) {
  dir.create(dataDir, recursive = TRUE)
}

### DOWNLOAD BIOSYS DATA -------------------------------------------------------

# Bulk downloads for entire dataset here:
# https://environment.data.gov.uk/ecology/explorer/downloads/

# Dataset documentation available here:
# https://www.data.gov.uk/dataset/3faf10d7-04bc-49e0-8377-61f75186d21d/freshwater-river-macroinvertebrate-surveys-biosys

# SURVEY DATA

# Set download url
url <- "https://environment.data.gov.uk/ecology/explorer/downloads/INV_OPEN_DATA.zip"

# Set file name
file_name <- "SurveyData"

# Download
download.file(url = url, destfile = paste0(dataDir, file_name, ".zip"), sep = "")

# Unzip (file names = "INV_OPEN_DATA_METRICS.csv",
# "INV_OPEN_DATA_SITE.csv", "INV_OPEN_DATA_TAXA.csv")
unzip(paste0(dataDir, "/", file_name, ".zip"), exdir = dataDir)

# TAXON INFO

# Set download url
url <- "https://environment.data.gov.uk/ecology/explorer/downloads/OPEN_DATA_TAXON_INFO.zip"

# Set file name
file_name <- "TaxonInfo"

# Download
download.file(url = url, destfile = paste0(dataDir, file_name, ".zip"), sep = "")

# Unzip (file name = "OPEN_DATA_TAXON_INFO.csv")
unzip(paste0(dataDir, "/", file_name, ".zip"), exdir = dataDir)

# Delete .zip files as no longer needed
unlink(paste0(dataDir,
              c("TaxonInfo.zip",
                "SurveyData.zip")))

### PROCESS BIOSYS DATA --------------------------------------------------------

# Read in separate files
invMetrics <- read.csv(paste0(dataDir, "INV_OPEN_DATA_METRICS.csv"))
invSite <- read.csv(paste0(dataDir, "INV_OPEN_DATA_SITE.csv"))
invTaxa <- read.csv(paste0(dataDir, "INV_OPEN_DATA_TAXA.csv"))
taxonInfo <- read.csv(paste0(dataDir, "OPEN_DATA_TAXON_INFO.csv"))

# !!! Add in this step while error in INV_OPEN_DATA_SITE.csv is fixed
# (duplicated SITE_ID)
invSite <- invSite %>%
  group_by(SITE_ID) %>%
  slice_head(n = 1) %>%
  ungroup

# Join files together
# N.B. Based on dataset documentation, these relationships can be used to join:
# INV_OPEN_DATA_SITE.SITE_ID = INV_OPEN_DATA_METRICS.SITE_ID
# INV_OPEN_DATA_METRICS.ANALYSIS_ID = INV_OPEN_DATA_TAXA.ANALYSIS_ID
# OPEN_DATA_TAXON_INFO.TAXON_LIST_ITEM_KEY = INV_OPEN_DATA_TAXA.TAXON_LIST_ITEM_KEY

invData <- 
  inner_join(invMetrics,
             invSite,
             # Use relationships above
             by = "SITE_ID", 
             # If duplicate columns, add ".y" to second duplicate
             suffix=c("",".y")) %>% 
  inner_join(invTaxa,
             .,
             by = "ANALYSIS_ID",
             suffix=c("",".y")) %>%
  inner_join(.,
             taxonInfo,
             by = "TAXON_LIST_ITEM_KEY",
             suffix=c("",".y")) %>%
  # Remove all colums with ".y" at end, i.e. duplicate columns
  select(-ends_with(".y")) 

# Remove large files no longer needed
rm(invMetrics, invSite, invTaxa, taxonInfo)
gc()

### SAVE BIOSYS DATA --------------------------------------------------------

# Save file
saveRDS(invData, file = paste0(dataDir, "invData.Rds"))
