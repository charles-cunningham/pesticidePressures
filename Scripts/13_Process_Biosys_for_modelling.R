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

# Read land cover data as 1km spatRast file for raster template
lcm2015 <- paste0(dataDir, "Raw/Land_cover_data/lcm2015gb25m.tif") %>%
  rast() %>%
  .[[1]] %>%
  terra::aggregate(., fact = 200)

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

# PROCESS VARIABLES --------------------------------------------------------------

### SCALE VARIABLES

# List variables to be scaled
modelVariables <- c(
  # Upstream variables
  "pesticideShannon",
  "pesticideToxicLoad",
  "insecticideToxicLoad",
  "herbicideToxicLoad",
  "fungicideToxicLoad",
  "fertiliser_k",
  "fertiliser_n",
  "fertiliser_p",
  "Urban",
  "Suburban",
  "Deciduous_woodland",
  "Coniferous_woodland",
  "cattle",
  "pigs",
  "sheep",
  "poultry",
  "totalArea",
  # Site variables
  "EDF_MEAN",
  "HS_HMS_RSB_SubScore",
  "HS_HQA",
  #"ALTITUDE",
  #"SLOPE",
  #"DIST_FROM_SOURCE",
  #"DISCHARGE",
  #"WIDTH",
  #"DEPTH",
  "BOULDERS_COBBLES",
  "PEBBLES_GRAVEL",
  "SAND",
  "SILT_CLAY" #,
 # "ALKALINITY"
 )

# Create additional scaled column for each modelVariables
for(variable in modelVariables) {
  
  # Create new scaled column name
  colName <- paste0(variable, "_scaled")
  
  # Assign scaled variable to new column
  invData[[colName]] <- scale(invData[[variable]])[,1]
}

### CONVERT IID EFFECTS TO FACTORS

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

### CONVERT RW EFFECTS TO DISCRETE DATA

# Loop through site variables used for random effects
for (i in c("ALTITUDE",
            "SLOPE",
            "DIST_FROM_SOURCE",
            "DISCHARGE",
            "ALKALINITY")) {
# Group to 10 evenly sized bins and assign to new column
  invData[, paste0(i, "_GRP")] <- INLA::inla.group(invData[[i]],
                                                  n = 10, 
                                                  method = "cut")
  
}

### CONVERT SAMPLING SITE ABIOTIC VARIABLES TO PCA ----------------------------------------------

# Carry out PCA
sitePCA <- invData %>%
  as_tibble(.) %>%
  select(#ALTITUDE_scaled,
         #SLOPE_scaled,
         #DIST_FROM_SOURCE_scaled,
         #DISCHARGE_scaled,
         #WIDTH_scaled,
         #DEPTH_scaled,
         BOULDERS_COBBLES_scaled,
         PEBBLES_GRAVEL_scaled,
         SAND_scaled,
         SILT_CLAY_scaled,
         #ALKALINITY_scaled
         ) %>%
  prcomp(.)

# Print summary
sitePCA; summary(sitePCA)

# Join to invData
invData <- cbind(invData, sitePCA$x)

### CORRELATION PLOTS WITH ALL VARIABLES ----------------------------------------------------------

# Create correlation data frame
corr_df <- invData %>%
  as_tibble(.) %>%
  select(pesticideShannon_scaled,
         pesticideToxicLoad_scaled,
         insecticideToxicLoad_scaled,
         herbicideToxicLoad_scaled,
         fungicideToxicLoad_scaled,
         fertiliser_k_scaled,
         fertiliser_n_scaled,
         fertiliser_p_scaled,
         Urban_scaled,
         Suburban_scaled,
         Deciduous_woodland_scaled,
         Coniferous_woodland_scaled,
         cattle_scaled,
         pigs_scaled,
         sheep_scaled,
         poultry_scaled,
         totalArea_scaled,
         EDF_MEAN_scaled,
         HS_HMS_RSB_SubScore_scaled,
         HS_HQA_scaled,
         PC1,
         PC2,
         # PC3,
         # PC4,
         # PC5,
         # PC6,
         # PC8,
         YEAR,
         MONTH_NUM)

# WITH WASTEWATER

# Filter out missing wastewater values and create correlation object
corPredictors <- filter(corr_df, !(is.na(EDF_MEAN_scaled))) %>%
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
corPredictors <- select(corr_df, !(EDF_MEAN_scaled)) %>%
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
  mutate(woodland = Deciduous_woodland_scaled + Coniferous_woodland_scaled) %>%
  # Residential
  mutate(residential = Urban_scaled + Suburban_scaled) %>%
  # Eutrophication risk
  mutate(eutroph = mean(fertiliser_n_scaled, fertiliser_p_scaled))

# PCA CHEMICAL INPUT ----------------------------------------------------------

# PCA eutrophication chemicals and pesticide load as very correlated
chemPCA <- invData %>%
  as_tibble(.) %>%
  select(
    pesticideToxicLoad_scaled,
    eutroph
  ) %>%
  prcomp(.)

# Change names to aid interpretation
dimnames(chemPCA$x)[[2]] <- c("chemicalApp", "lessPesticide")

# Print summary
chemPCA; summary(chemPCA)

# Join to invData
invData <- cbind(invData, chemPCA$x)

### CORRELATION PLOTS WITH AGGREAGTED VARIABLES -------------------------------------------------------

# Create correlation data frame
corr_df <- invData %>%
  as_tibble(.) %>%
  select(pesticideShannon_scaled,
         chemicalApp,
         lessPesticide,
         residential,
         woodland,
         cattle_scaled,
         pigs_scaled,
         sheep_scaled,
         poultry_scaled,
         EDF_MEAN_scaled,
         HS_HMS_RSB_SubScore_scaled,,
         HS_HQA_scaled,
         totalArea_scaled,
         PC1,
         PC2,
        #  PC3,
         # PC4,
         # PC5,
         # PC6,
         # PC7,
         # PC8,
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

# Convert to raster
england_R <- rasterize(england, lcm2015,
                       touches = TRUE,
                       cover = TRUE)

# Restrict to cells with more than a small fraction (10% of land cover)
england_R <- ifel(england_R > 0.1,
             yes = 1,
             no = NA)

# Convert to sf object
england_sf <- as.polygons(england_R) %>% # Convert to polygon 
  st_as_sf(.) # Convert to sf object for smoothing

# Smooth
# N.B. This is used to created mesh and functions as modelling boundary (domain)
englandSmooth <- england_sf %>% 
  smoothr::smooth(., method = "chaikin") %>%
  smoothr::fill_holes(., threshold = Inf)

### SAVE DATASETS ---------------------------------------------------------------

# Save file
saveRDS(invData, file = paste0(dataDir, "Processed/Biosys/invData_forModel.Rds"))
saveRDS(englandSmooth, file = paste0(dataDir, "Raw/Country_data/EnglandSmooth.Rds"))
