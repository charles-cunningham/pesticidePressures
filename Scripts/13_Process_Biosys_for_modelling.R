# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process Biosys data for inlabru models
#
# Script Description: 

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(sf)
library(corrplot)

### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
plotDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Plots/"

# Create processed Biosys data folder
lapply(paste0(dataDir, "Processed/Species"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

### LOAD DATA ------------------------------------------------------------------

# Load Biosys data
invData <- readRDS(paste0(dataDir, "Processed/Biosys/invDataSpatial.Rds"))

# PROCESS DATA STRUCTURE -------------------------------------------------------

# Remove spaces from names for inlabru
names(invData) <- gsub(" ", "_", names(invData) )

### PROCESS DATES --------------------------------------------------------------

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

# FILTER VARIABLES -------------------------------------------------------------

# Filter temporally ( filter years after 2010)
invData <- invData %>%
  filter(YEAR > 2010)

# Rescale months and years
invData$MONTH_NUM <- invData$MONTH_NUM - (min(invData$MONTH_NUM) - 1)
invData$YEAR <- invData$YEAR - (min(invData$YEAR) - 1)

# Filter available  data
invData <- invData %>%
  # Remove rows with no upstream data
  filter(!(is.na(pesticideLoad))) %>%
  # Remove rows with no site data
  filter(!(is.na(HS_HMS_RSB_SubScore))) %>%
  # Remove_rows with missing BIOSYS site data
  filter(!(is.na(DEPTH)))%>%
  filter(!(is.na(SAND))) %>%
  filter(!(is.na(ALKALINITY))) %>%
  filter(!(is.na(SILT_CLAY))) 

### PROCESS TAXONOMY -----------------------------------------------------------

# Change species names to be file friendly
invData$TAXON_GROUP_NAME <- invData$TAXON_GROUP_NAME %>%
  # Remove "insect - " prefix
  gsub("insect - ", "", .) %>%
  # Remove spaces
  gsub(" ", "_", .)

# Change taxanomic group names to be file friendly
invData$TAXON <- invData$TAXON %>%
  # Remove spaces
  gsub(" ", "_", .) %>%
  # Remove slashes
  gsub("/", "-", .)

### AGGREGATE VARIABLES --------------------------------------------------------

# Woodland
invData$woodland <- invData$Deciduous_woodland + invData$Coniferous_woodland

# Residential
invData$residential <- invData$Urban + invData$Suburban

# Eutrophication risk
invData$eutroph <- mean(invData$fertiliser_n + invData$fertiliser_p)

# MODIFY UPSTREAM VARIABLES TO PER AREA VALUES----------------------------------

# Divide upstream variables by area (excluding diversity)
for(variable in c(
  "eutroph",
  "residential",
  "pesticideToxicLoad",
  "cattle",
  "pigs",
  "sheep",
  "poultry",
  "woodland")) {
  
  # Create new scaled column name
  colName <- paste0(variable, "_PerArea")
  
  # Assign scaled variable to new column
  invData[[colName]] <- invData[[ variable]] / invData[[ "totalArea"]]
  
}

### CONVERT SITE VARIABLES TO PCA ----------------------------------------------

sitePCA <- invData %>%
  as_tibble(.) %>%
  select(ALTITUDE,
         SLOPE,
         DIST_FROM_SOURCE,
         DISCHARGE,
         WIDTH,
         DEPTH,
         BOULDERS_COBBLES,
         PEBBLES_GRAVEL,
         SAND,
         SILT_CLAY,
         ALKALINITY) %>%
  prcomp(.)

summary(sitePCA)

invData <- cbind(invData, sitePCA$x)

### CORRELATION PLOTS ----------------------------------------------------------

# Create correlation data frame
corr_df <- invData %>%
  as_tibble(.) %>%
  select(pesticideShannon,
         pesticideToxicLoad_PerArea,
         eutroph_PerArea,
         residential_PerArea,
         woodland_PerArea,
         cattle_PerArea,
         pigs_PerArea,
         sheep_PerArea,
         poultry_PerArea,
         EDF_MEAN,
         HS_HMS_RSB_SubScore,
         HS_HQA,
         PC1,
         PC2,
         PC3,
         PC4,
         YEAR,
         MONTH_NUM)

# WITH WASTEWATER

# Filter out missing wastewater values and create correlation object
corPredictors <- filter(corr_df, !(is.na(EDF_MEAN))) %>%
  cor(.) 

png(filename = paste0(plotDir, 'corr_wastewater.png'),
    width = 40, height = 30, units = "cm", res = 600)
corrplot::corrplot(corPredictors,
                   type = "upper", order = "original", diag = FALSE,
                   method = "number", addCoef.col="white", tl.col = "black",
                   tl.srt = 45, tl.cex = 1)
dev.off()

# WITHOUT WASTEWATER

# Unselect wasterwater column and create correlation object
corPredictors <- select(corr_df, !(EDF_MEAN)) %>%
  cor(.) 

png(filename = paste0(plotDir, 'corr_noWastewater.png'),
    width = 40, height = 30, units = "cm", res = 600)
corrplot::corrplot(corPredictors,
                   type = "upper", order = "original", diag = FALSE,
                   method = "number", addCoef.col="white", tl.col = "black",
                   tl.srt = 45, tl.cex = 1)
dev.off()

# SCALE VARIABLES --------------------------------------------------------------

# List variables to be scaled
modelVariables <- c(
  # Upstream variables
  "eutroph_PerArea",
  "residential_PerArea",
  "woodland_PerArea",
  "pesticideShannon",
  "pesticideToxicLoad_PerArea",
  "cattle_PerArea",
  "pigs_PerArea",
  "sheep_PerArea",
  "poultry_PerArea",
  # Site variables
  "EDF_MEAN",
  "HS_HMS_RSB_SubScore",
  "HS_HQA",
  "PC1",
  "PC2",
  "PC3",
  "PC4")

# Create additional scaled column for each modelVariables
for(variable in modelVariables) {
  
  # Create new scaled column name
  colName <- paste0(variable, "_scaled")
  
  # Assign scaled variable to new column
  invData[[colName]] <- scale(invData[[variable]])[,1]
}

# Convert categorical variables for nested random effects to factors
# (Convert catchment and water body to numeric to save memory)
invData$BASIN_F <- as.factor(invData$REPORTING_AREA)
invData$CATCHMENT_F <- paste( invData$REPORTING_AREA,
                              invData$CATCHMENT) %>%
  as.factor() %>%
  as.numeric() %>%
  as.factor()
invData$WATER_BODY_F <-paste( invData$REPORTING_AREA,
                              invData$CATCHMENT,
                              invData$WATER_BODY) %>%
  as.factor() %>%
  as.numeric() %>%
  as.factor()

### SAVE DATASET ---------------------------------------------------------------

# Save file
saveRDS(invData, file = paste0(dataDir, "Processed/Biosys/invData_forModel.Rds"))
