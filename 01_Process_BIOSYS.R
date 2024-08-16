# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process BIOSYS data
#
# Script Description:
#

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)

### DATA MANAGEMENT ------------------------------------------------------------

# Set data directory
# If working on DASH: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Raw/Species"
# If working locally: "../Data/Raw/Species"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/"

# Create directory if it doesn't exist
if (!file.exists(dataDir)) {
  dir.create(dataDir, recursive = TRUE)
}

### DOWNLOAD BIOSYS DATA -------------------------------------------------------

# Dataset documentation available here:
# https://www.data.gov.uk/dataset/3faf10d7-04bc-49e0-8377-61f75186d21d/freshwater-river-macroinvertebrate-surveys-biosys

# Bulk downloads for entire dataset here:
# https://environment.data.gov.uk/ecology/explorer/downloads/

### DOWNLOAD TAXON INFO DATA FROM DATA.GOV.UK SITE

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

# Delete .zip files as no longer needed
unlink(paste0(dataDir,
              c("TaxonInfo.zip",
                "SurveyData.zip")))

# TAXON INFO

# Set download url
url <- "https://environment.data.gov.uk/ecology/explorer/downloads/OPEN_DATA_TAXON_INFO.zip"

# Set file name
file_name <- "TaxonInfo"

# Download
download.file(url = url, destfile = paste0(dataDir, file_name, ".zip"), sep = "")

# Unzip (file name = "OPEN_DATA_TAXON_INFO.csv")
unzip(paste0(dataDir, "/", file_name, ".zip"), exdir = dataDir)

### READ IN BIOSYS DATA --------------------------------------------------------

# Read in separate files
invMetrics <- read.csv(paste0(dataDir, "INV_OPEN_DATA_METRICS.csv"))
invSite <- read.csv(paste0(dataDir, "INV_OPEN_DATA_SITE.csv"))
invTaxa <- read.csv(paste0(dataDir, "INV_OPEN_DATA_TAXA.csv"))
taxonInfo <- read.csv(paste0(dataDir, "OPEN_DATA_TAXON_INFO.csv"))

# Join files together
# N.B. Based on dataset documentation, the following can be used to join:
# INV_OPEN_DATA_SITE.SITE_ID = INV_OPEN_DATA_METRICS.SITE_ID
# INV_OPEN_DATA_METRICS.ANALYSIS_ID = INV_OPEN_DATA_TAXA.ANALYSIS_ID
# OPEN_DATA_TAXON_INFO.TAXON_LIST_ITEM_KEY = INV_OPEN_DATA_TAXA.TAXON_LIST_ITEM_KEY

#test <- left_join(invMetrics, invSite, by = "SITE_ID")

