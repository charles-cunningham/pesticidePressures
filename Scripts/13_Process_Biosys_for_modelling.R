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
library(terra)

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

# England boundary
england <- readRDS(paste0(dataDir, "Raw/Country_data/England.Rds"))

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
  # Remove sites with upstream area greater than 200km^2 (200,000,000m^2) 
  # to avoid increasing uncertainty with larger areas
  # area is in 25x25m cells
  filter((totalArea * 25*25) < 200000000) %>%
  # Filter rows to spring
  filter(MONTH_NUM >=3 & MONTH_NUM <= 5) %>%
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

# PROCESS FACTOR VARIABLES -----------------------------------------------------

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
invData$SITE_ID <- as.factor(invData$SITE_ID)

### CONVERT SAMPLING SITE ABIOTIC VARIABLES TO PCA -----------------------------

abioticVariables <- c(
  "ALTITUDE",
  "SLOPE",
  "DIST_FROM_SOURCE",
  "DISCHARGE",
  "WIDTH",
  "DEPTH",
  "BOULDERS_COBBLES",
  "PEBBLES_GRAVEL",
  "SAND",
  "SILT_CLAY",
  "ALKALINITY"
)

# Create additional scaled column for each modelVariables
for(variable in abioticVariables) {
  
  # Create new scaled column name
  colName <- paste0(variable, "_scaled")
  
  # Assign scaled variable to new column
  invData[[colName]] <- scale(invData[[variable]])[,1]
}

# Carry out PCA
sitePCA <- invData %>%
  as_tibble(.) %>%
  select(ALTITUDE_scaled,
         SLOPE_scaled,
         DIST_FROM_SOURCE_scaled,
         DISCHARGE_scaled,
         WIDTH_scaled,
         DEPTH_scaled,
         BOULDERS_COBBLES_scaled,
         PEBBLES_GRAVEL_scaled,
         SAND_scaled,
         SILT_CLAY_scaled,
         ALKALINITY_scaled
         ) %>%
  prcomp(.)

# Print summary
sitePCA; summary(sitePCA)

# Join to invData
invData <- cbind(invData, sitePCA$x)

###  PCA evaluation
# Code from https://en.proft.me/2016/11/15/principal-component-analysis-pca-r/)

evplot = function(ev) {
  # Broken stick model (MacArthur 1957)
  n = length(ev)
  bsm = data.frame(j=seq(1:n), p=0)
  bsm$p[1] = 1/n
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p = 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}
ev = sitePCA$sdev^2
evplot(ev)

### CORRELATION PLOTS WITH ALL VARIABLES ----------------------------------------------------------

# Create correlation data frame
corr_df <- invData %>%
  as_tibble(.) %>%
  select(pesticideShannon,
         pesticideToxicLoad,
         insecticideToxicLoad,
         herbicideToxicLoad,
         fungicideToxicLoad,
         fertiliser_k,
         fertiliser_n,
         fertiliser_p,
         Urban,
         Suburban,
         Deciduous_woodland,
         Coniferous_woodland,
         cattle,
         pigs,
         sheep,
         poultry,
         totalArea,
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

### AGGREGATE VARIABLES --------------------------------------------------------

invData <- invData %>%
  rowwise() %>%
  # Woodland
  mutate(woodland = Deciduous_woodland + Coniferous_woodland) %>%
  # Residential
  mutate(residential = Urban + Suburban) %>%
  # Eutrophication risk
  mutate(eutroph = mean(fertiliser_n, fertiliser_p))

# PCA CHEMICAL INPUT ----------------------------------------------------------

# # PCA eutrophication chemicals and pesticide load as very correlated
# chemPCA <- invData %>%
#   as_tibble(.) %>%
#   select(
#     pesticideToxicLoad_scaled,
#     eutroph
#   ) %>%
#   prcomp(.)
# 
# # Change names to aid interpretation
# dimnames(chemPCA$x)[[2]] <- c("chemicalApp", "lessPesticide")
# 
# # Print summary
# chemPCA; summary(chemPCA)
# 
# # Join to invData
# invData <- cbind(invData, chemPCA$x)

# SCALE VARIABLES --------------------------------------------------------------

# List variables to be scaled
modelVariables <- c(
  # Upstream variables
  "pesticideShannon",
  "pesticideToxicLoad",
  "eutroph",
  "woodland",
  "residential",
  "cattle",
  "pigs",
  "sheep",
  "poultry",
  "totalArea",
  # Site variables
  "EDF_MEAN",
  "HS_HMS_RSB_SubScore",
  "HS_HQA"
)

# Create additional scaled column for each modelVariables
for(variable in modelVariables) {
  
  # Create new scaled column name
  colName <- paste0(variable, "_scaled")
  
  # Assign scaled variable to new column
  invData[[colName]] <- scale(invData[[variable]])[,1]
}

### CORRELATION PLOTS WITH AGGREAGTED VARIABLES -------------------------------------------------------

# Create correlation data frame
corr_df <- invData %>%
  as_tibble(.) %>%
  select(pesticideShannon_scaled,
         pesticideToxicLoad_scaled,
         eutroph_scaled,
         residential_scaled,
         woodland_scaled,
         cattle_scaled,
         pigs_scaled,
         sheep_scaled,
         poultry_scaled,
         EDF_MEAN_scaled,
         HS_HMS_RSB_SubScore_scaled,
         HS_HQA_scaled,
         totalArea_scaled,
         PC1,
         PC2,
         PC3,
         PC4,
         YEAR,
         MONTH_NUM)

# WITH WASTEWATER

# Filter out missing wastewater values and create correlation object
corPredictors <- filter(corr_df, !(is.na(EDF_MEAN_scaled))) %>%
  cor(.) 

png(filename = paste0(plotDir, 'corr_wastewater_agg.png'),
    width = 40, height = 30, units = "cm", res = 600)
corrplot::corrplot(corPredictors,
                   type = "upper", order = "original", diag = FALSE,
                   method = "number", addCoef.col="white", tl.col = "black",
                   tl.srt = 45, tl.cex = 1)
dev.off()

# WITHOUT WASTEWATER

# Unselect wasterwater column and create correlation object
corPredictors <- select(corr_df, !(EDF_MEAN_scaled)) %>%
  cor(.) 

png(filename = paste0(plotDir, 'corr_noWastewater_agg.png'),
    width = 40, height = 30, units = "cm", res = 600)
corrplot::corrplot(corPredictors,
                   type = "upper", order = "original", diag = FALSE,
                   method = "number", addCoef.col="white", tl.col = "black",
                   tl.srt = 45, tl.cex = 1)
dev.off()

### PROCESS ENGLAND ------------------------------------------------------------

### Remove small islands

# Disaggregate
england <- disagg(england)

# Calculate area
england$area_sqkm <- expanse(england, unit = "km")

# Remove polygons with < 50km ^2 area
england <- england[england$area_sqkm > 50]

# Aggregate back
england <- aggregate(england)

### Smooth

# N.B. This is used to created mesh and functions as modelling boundary (domain)
englandSmooth <- st_as_sf(england) %>%
  smoothr::smooth(., method = "chaikin") %>%
  smoothr::fill_holes(., threshold = Inf)

### SAVE DATASETS ---------------------------------------------------------------

# Save file
saveRDS(invData, file = paste0(dataDir, "Processed/Biosys/invData_forModel.Rds"))
saveRDS(englandSmooth, file = paste0(dataDir, "Raw/Country_data/EnglandSmooth.Rds"))
