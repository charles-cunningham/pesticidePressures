# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Validate overland flow data
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

# Create processed flow data folder
lapply(paste0(dataDir, "Processed/Screen"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

### LOAD DATA ------------------------------------------------------------------

# Load flow data
flowChemData <- readRDS(paste0(dataDir,
                               "Processed/Flow/Flow_chem_data.Rds"))

# Load LC-MS (liquid chromatography-mass spectrometry) screen data
lcms <- paste0(dataDir,
               "Raw/Monitoring_screen_data/LCMS Target and Non-Targeted Screening.csv") %>%
  read.csv()

# Load GC-MS (gas chromatography-mass spectrometry) screen data

gcms <- paste0(dataDir,
               "Raw/Monitoring_screen_data/GCMS Target and Non-Targeted Screening.csv") %>%
  read.csv()

# Join together for processing (can still be separated by 'method' column)
# N.B. From README metadata - "when overlap occurs (ie. same compounds can be 
# detected by both methods) the LC-MS is able to detect using a much lower limit 
# of detection than GC-MS."
screenData <- rbind(lcms, gcms)

# Remove redundant objects
rm(lcms, gcms)

### PROCESS SCREEN DATA --------------------------------------------------------

# FILTER DATA

# Filter:
# - sampling material to rivers and running freshwater
# - sampling type to freshwater only
# - year to only 2010-2019 to be in similar time frame to flow data
screenData <- screenData %>%
  filter(SMC_DESC == "RIVER / RUNNING SURFACE WATER") %>%
  filter(grepl("FRESHWATER",.$SPT_DESC)) %>%
  filter(year >=2010 & year < 2020)

# REMOVE COLUMNS NOT NEEDED

Sample_Site_ID keep
SMPT_TYPE remove
SPT_DESC filter then remove
SAMP_ID keep
SAMP_MATERIAL remove
SMC_DESC filter and remove
SAMP_PURPOSE_CODE remove
PURP_DESC keep
MEAS_DETERMINAND_CODE
ARE_CODE remove 
ARE_DESC
Sample_datetime
year
Screening_Method_Details
CAS_Number
unit
Concentration
Spectral_Fit
Compound_Name
LOD
method
less_than
SMPT_LONG_NAME
SMPT_EASTING keep
SMPT_NORTHING keep 
Latitude remove
Longitude remove
OPCATNAME keep
COUNTRY remove (all england)
