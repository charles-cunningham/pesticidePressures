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
           y = log(APPLICATION_MOD))) +
  geom_point()

# Check correlation
cor(gcmsData$Concentration, gcmsData$APPLICATION_MOD)

# LCMS

# Plot relationship between concentration and application per area
ggplot(data = lcmsData,
       aes(x = log(Concentration), 
           y = log(APPLICATION_MOD))) +
  geom_point()

# Check correlation
cor(lcmsData$Concentration, lcmsData$APPLICATION_MOD)

# TRUNCATED DATA MODEL ---------------------------------------------------------

# GCMS 

#
gcmsTruncated<- gcmsData %>%
  mutate(truncatedConc = ifelse(Concentration> LOD,
                                Concentration,
                                LOD))

gcmsTruncModel <- crch::trch ( truncatedConc ~ log(APPLICATION_MOD),
                     truncated = TRUE,
                     left =  gcmsTruncated$LOD,
                     link.scale = "log",
                     data = gcmsTruncated)


summary(gcmsTruncModel)

# LCMS

#
lcmsTruncated <- lcmsData %>%
  mutate(truncatedConc = ifelse(Concentration > LOD,
                                Concentration,
                                LOD))

lcmsTruncModel <- crch::trch ( truncatedConc ~ log(APPLICATION_MOD),
                               truncated = TRUE,
                               left =  lcmsTruncated$LOD,
                               link.scale = "log",
                               data = lcmsTruncated)


summary(lcmsTruncModel)


### MIXED-EFFECTS MODEL --------------------------------------------------------

# For every compound, aggregate all concentrations at every site by the median
gcmsSummary <- gcmsData %>%
  group_by(Sample_Site_ID, PESTICIDE_MOD, OPCAT_NAME, APPLICATION_MOD) %>%
  summarise(ConcentrationMedian = median(Concentration)) %>%
  ungroup()


# Fit mixed model with concentration as response
gcmsModMixed <- lme4::glmer(ConcentrationMedian ~ 
                              # Fixed log-transformed estimated application
                              log(APPLICATION_MOD) + 
                              # Pesticide name random effect
                              (1| PESTICIDE_MOD) +
                              # # Site random effect
                               (1 | Sample_Site_ID) ,
                            family = Gamma(link = "log"),
                            data = gcmsSummary) 

# Checl R^2 values
MuMIn::r.squaredGLMM(gcmsModMixed)

# Check fixed effect
summary(gcmsModMixed)

# LCMS

# For every compound, aggregate all concentrations at every site by the median
lcmsSummary <- lcmsData %>%
  group_by(Sample_Site_ID, PESTICIDE_MOD, OPCAT_NAME, APPLICATION_MOD) %>%
  summarise(ConcentrationMedian = median(Concentration)) %>%
  ungroup()


# Fit mixed model with concentration as response
lcmsModMixed <- lme4::glmer(ConcentrationMedian ~ 1 +
                              # Fixed log-transformed estimated application
                              log(APPLICATION_MOD) + 
                              # Pesticide name random effect
                              (1| PESTICIDE_MOD) +
                              # Site random effect
                              (1 | Sample_Site_ID) ,
                            family = Gamma(link = "log"),
                            data = lcmsSummary) 


# Check R^2 values
MuMIn::r.squaredGLMM(lcmsModMixed)

# Check fixed effect
summary(lcmsModMixed)
