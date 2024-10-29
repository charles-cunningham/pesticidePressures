# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Model Biosys data using inlabru
#
# Script Description: 

# INSTALL NON-CRAN PACKAGES -----------------------------------

#Run this code once

# # Install INLA
# install.packages("INLA",
#                  repos=c(getOption("repos"),INLA="http://inla.r-inla-download.org/R/stable"),
#                  dep=TRUE)
# 
# # Install inlabru
# install.packages("inlabru")
# 
# # Install BRCmap
# devtools::install_github("colinharrower/BRCmap")
# 
# # Needed to run on HPC
# inla.binary.install()
# 
# # Other packages as required from CRAN, i.e install.packages()

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(INLA)
library(terra)

### DATA MANAGEMENT ------------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Biosys/"
# If working locally: "../Data/Processed/Biosys/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Biosys/"

# Load Biosys data
invData <- readRDS(paste0(dataDir, "invDataSpatial.Rds"))

###        ---------------------------------------------------------------------


test <- invData %>%
  mutate(SAMPLE_DATE = as.POSIXct(SAMPLE_DATE, format = "%d/%m/%Y"))

test$year <- test$SAMPLE_DATE %>%
  format(., format="%Y") %>% 
  as.numeric(.)


