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

### LOAD DATA ------------------------------------------------------------------

# Load screen data
screenData <- readRDS(paste0(dataDir,
                             "/Processed/Screen/Screen_data.Rds"))

# Read pesticide data
pestData <- read.csv(paste0(dataDir,
                            "Raw/Pesticide_data/Pesticide_properties.csv"))

### PROCESS PESTICIDE DATA -----------------------------------------------------

# Create joint wieght
pestData$Exposure <- #pestData$Leachability..GUS_leaching_potential_index. %>%
  #ifelse(is.na(.), mean(., na.rm = TRUE), .) %>%
 # ifelse(. <0.01, 0.01, .)# %>%
  #`*` (
    ifelse(is.na(pestData$Persistence_.DT.._.typical..),
            mean(pestData$Persistence_.DT.._.typical.., na.rm = TRUE),
            pestData$Persistence_.DT.._.typical..)
 #  )
           
   screenData$ExposureWeight <- NA 
   
for(pesticide in pestData$Pesticide) {

screenData[which(screenData$PESTICIDE_MOD %in% pesticide), "ExposureWeight"] <- 
     pestData[which(pestData$Pesticide  %in% pesticide), "Exposure"]
   
   }

### FILTER SCREEN DATA ---------------------------------------------------------

# Filter missing values, and convert to tibble
screenData <- screenData %>%
  # Remove sites with no name matched modelled data
  filter(!is.na(PESTICIDE_MOD)) %>%
  # Remove sites with no Limit of Detection data
  filter(!is.na(LOD)) %>%
  # Remove sites with concentrations below detection level
  filter(less_than == FALSE) %>%
  as_tibble()

# Separate out into GCMS and LCMS datasets to test separately
gcmsData <- screenData[screenData$method == "GCMS",]
lcmsData <- screenData[screenData$method == "LCMS",]

# DESCRIPTIVE STATS ------------------------------------------------------------

# GCMS

# Plot relationship between concentration and application per area
ggplot(data = gcmsData,
       aes(x = log(Concentration), 
           y = log(APP_PER_AREA_MOD *ExposureWeight))) +
  geom_point()

# Check correlation# Check correlationExposure
cor(gcmsData$Concentration, gcmsData$APP_PER_AREA_MOD * gcmsData$ExposureWeight)

# LCMS

# Plot relationship between concentration and application per area
ggplot(data = lcmsData,
       aes(x = log(Concentration), 
           y = log(APP_PER_AREA_MOD * ExposureWeight))) +
  geom_point()

# Check correlation
cor(lcmsData$Concentration, lcmsData$APP_PER_AREA_MOD * lcmsData$ExposureWeight)

# TRUNCATED DATA MODEL ---------------------------------------------------------

# GCMS 

#
gcmsTruncated<- gcmsData %>%
  filter(APP_PER_AREA_MOD > 0 ) %>%
  mutate(truncatedConc = ifelse(Concentration> LOD,
                                Concentration,
                                LOD))

gcmsTruncModel <- crch::trch ( truncatedConc ~ log(APP_PER_AREA_MOD*ExposureWeight),
                     truncated = TRUE,
                     left =  gcmsTruncated$LOD,
                     link.scale = "log",
                     data = gcmsTruncated)


summary(gcmsTruncModel)

# LCMS

#
lcmsTruncated <- lcmsData %>%
  filter(APP_PER_AREA_MOD > 0 ) %>%
  mutate(truncatedConc = ifelse(Concentration > LOD,
                                Concentration,
                                LOD))

lcmsTruncModel <- crch::trch ( truncatedConc ~ log(APP_PER_AREA_MOD*ExposureWeight),
                               truncated = TRUE,
                               left =  lcmsTruncated$LOD,
                               link.scale = "log",
                               data = lcmsTruncated)


summary(lcmsTruncModel)


### MIXED-EFFECTS MODEL --------------------------------------------------------

# GCMS

# Fit mixed model with concentration as response
gcmsModMixed <- lme4::glmer(Concentration ~ 1 +
                              # Fixed log-transformed estimated application
                              log(APP_PER_AREA_MOD*ExposureWeight) + 
                              # Pesticide name random effect
                              (1| PESTICIDE_MOD) +
                              # # Site random effect
                               (1 | Sample_Site_ID) ,
                            family = Gamma(link = "log"),
                            data = gcmsData %>% filter(APP_PER_AREA_MOD>0)) 

# Checl R^2 values
MuMIn::r.squaredGLMM(gcmsModMixed)

# Check fixed effect
summary(gcmsModMixed)

# LCMS

# Fit mixed model with concentration as response
lcmsModMixed <- lme4::glmer(Concentration ~ 1 +
                              # Fixed log-transformed estimated application
                              log(APP_PER_AREA_MOD*ExposureWeight) + 
                              # Pesticide name random effect
                              (1| PESTICIDE_MOD) +
                              # Site random effect
                              (1 | Sample_Site_ID) ,
                            family = Gamma(link = "log"),
                            data = lcmsData) 


# Check R^2 values
MuMIn::r.squaredGLMM(lcmsModMixed)

# Check fixed effect
summary(lcmsModMixed)
