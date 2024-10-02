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

### DATA MANAGEMENT ------------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Raw/Biosys/"
# If working locally: "../Data/Raw/Biosys/"
dataDir <- "../Data/Raw/Biosys/"

# Load Biosys data
invData <- readRDS(paste0(dataDir, "invData.Rds"))

###        ---------------------------------------------------------------------
