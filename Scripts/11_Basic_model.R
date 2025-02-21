# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Model Biosys data using inlabru
#
# Script Description: 

# INSTALL INLA AND INLABRU PACKAGES --------------------------------------------

# # Update matrix package first if using R version 4.2.2
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz",
#                  repos=NULL, type="source")
# 
# # Install INLA
# install.packages("INLA",
#                  repos=c(getOption("repos"),INLA="http://inla.r-inla-download.org/R/stable"),
#                  dep=TRUE)
# 
# # Needed to run on LINUX machine
# library(INLA); inla.binary.install()
# 
# # Install inlabru
# install.packages("inlabru")
# 
# # Other packages as required from CRAN, i.e install.packages()

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(INLA)
library(inlabru)
library(sf)

### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Biosys/"
# If working locally: "../Data/Processed/Biosys/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Biosys/"

### LOAD DATA ------------------------------------------------------------------

# Load Biosys data
invData <- readRDS(paste0(dataDir, "invDataSpatial.Rds"))

#####
hist(invData$pesticideLoad)
hist(scale(invData$pesticideLoad))

### (TEMPORALLY) FILTER DATA ---------------------------------------------------

# PROCESS DATES

# Covert date to POSIX
invData <- invData %>%
  mutate(SAMPLE_DATE = as.POSIXct(SAMPLE_DATE, format = "%d/%m/%Y"))

# Create year column
invData$YEAR <- invData$SAMPLE_DATE %>%
  format(., format="%Y") %>% 
  as.numeric(.)

# Create month name column
invData$MONTH_NAME <- invData$SAMPLE_DATE %>%
  format(., format="%B")

# Create month integer column
invData$MONTH_NUM <- invData$SAMPLE_DATE %>%
  format(., format="%m") %>%
  as.numeric()

# Create week column
invData$WEEK <- invData$SAMPLE_DATE %>%
  format(., format="%V") %>%
  as.numeric() 
  
# FILTER DATES

# Filter years and season between 2015 and 2020
invData <- invData %>%
  # Filter years between 2010 and 2019 inclusive
  filter(YEAR >= 2010 & YEAR < 2020) %>%
  # Filter to spring and summer (March, April, May, June, July, August)
  filter(MONTH_NAME %in% c("March", "April", "May", "June", "July", "August"))

# Rescale months
invData$MONTH_NUM <- invData$MONTH_NUM - (min(invData$MONTH_NUM) - 1)

# FILTER ON OTHER VARIABLES ----------------------------------------------------

# Remove rows with no upstream data
invData <- invData %>%
  filter(!(is.na(pesticideLoad)))

# Remove rows with no abundance data
invData <- invData %>%
  filter(!(is.na(TOTAL_ABUNDANCE)))

### PROCESS TAXONOMY -----------------------------------------------------------

### todo: TO BE CREATED USING WILKES METHODOLOGY!!!!!!!!!!


### SEPARATE INTO SEPERATE SPECIES ---------------------------------------------

# Remove spaces from layer names



# Start loop here
for( i in unique(invData$TAXON_GROUP_NAME)[1]) {

  taxaData <- invData %>%
    filter(TAXON_GROUP_NAME == i)

  # Remove spaces from names
  names(taxaData) <- gsub(" ", "_", names(taxaData) )
  
### ZIP MODEL 

  # Priors for random effects
  seasonHyper <- list(theta = list(prior="pc.prec",
                                   param=c(0.5, 0.01)))
  
  comps <-  TOTAL_ABUNDANCE~ 
    pesticideDiv(pesticideShannon, model = "linear") +
    pesticideLoad(log(pesticideToxicLoad), model = "linear") +
    N(log(fertiliser_n), model = "linear") +
    P(log(fertiliser_p), model = "linear") +
    K(log(fertiliser_k), model = "linear") +
    upstream(log(totalArea), model = "linear") +
    arable(Arable/1000000, model = "linear") +
    grass(Improved_grassland/1000000, model = "linear") +
    urban(Urban/1000000, model = "linear") +
    month(main = MONTH_NUM,
         model = "seasonal",
         season.length = 6,
         hyper = seasonHyper) +
    wb(WATER_BODY, model = "iid") +
    Intercept(1)
  
  model <- bru( components = comps,
      family = "poisson", 
      data = taxaData,
      options=list(control.compute = list(waic = TRUE,
                                          dic = TRUE,
                                          cpo = TRUE),
                   verbose = TRUE))

}
taxaData$Arable %>% hist()

