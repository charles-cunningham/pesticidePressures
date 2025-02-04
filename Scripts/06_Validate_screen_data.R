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

### VALIDATE FLOW DATA WITH SCREEN DATA ----------------------------------------

# FILTER SCREEN DATA

# Filter missing values, and convert to tibble
screenData <- screenData %>%
  # Remove sites with no name matched modelled data
  filter(!is.na(PESTICIDE_MOD)) %>%
  # Remove sites with concentrations below detection level
  filter(less_than == FALSE) %>%
  as_tibble()

# Separate out into GCMS and LCMS datasets to test separately
gcmsData <- screenData[screenData$method == "GCMS",]
lcmsData <- screenData[screenData$method == "LCMS",]

# MIXED-EFFECTS GLM

# GCMS 

# For every compound, aggregate all concentrations at every site by the median
gcmsSummary <- lcmsData %>%
  group_by(Sample_Site_ID, PESTICIDE_MOD, OPCAT_NAME, APPLICATION_MOD) %>% 
  summarise(ConcentrationMedian = median(Concentration)) %>%
  ungroup()


# Plot relationship
ggplot(data = gcmsSummary,
       aes(x = log(ConcentrationMedian), 
           y = log(APPLICATION_MOD))) +
  geom_point()

cor(gcmsSummary$ConcentrationMedian, gcmsSummary$APPLICATION_MOD)

gcmsModMixed <- lme4::glmer(ConcentrationMedian ~ log(APPLICATION_MOD) + 
                              (1 | PESTICIDE_MOD) +
                              (1 | OPCAT_NAME) ,
                            family = Gamma(link = "log"),
                            data = gcmsSummary) 


MuMIn::r.squaredGLMM(gcmsModMixed)

summary(gcmsModMixed)

gcmsData$Compound_Name %>% unique() %>% sort()
gcmsData$PESTICIDE_MOD %>% unique() %>% sort()
gcmsData[gcmsData$Compound_Name == "2,4-Dichlorophenol",]


# LCMS


# names not matching well



# Analysis 2










summary(testModMixed)
hist(log(testData$MODELLED_APPLICATION))
hist(log(testData$MODELLED_CONCENTRATION))

testData[testData$Concentration == 3070, ]
testData$unit %>% unique()
