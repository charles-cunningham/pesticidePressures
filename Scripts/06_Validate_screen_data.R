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
gcmsSummary <- gcmsData %>%
  group_by(Sample_Site_ID, PESTICIDE_MOD, OPCAT_NAME, APP_PER_AREA_MOD) %>% 
  summarise(ConcentrationMedian = median(Concentration)) %>%
  ungroup()

# Plot relationship
ggplot(data = gcmsSummary,
       aes(x = log(ConcentrationMedian), 
           y = log(APP_PER_AREA_MOD))) +
  geom_point()

# Fit mixed model with concentration as response
gcmsModMixed <- lme4::glmer(ConcentrationMedian ~ 
                              # Fixed log-transformed estimated application
                              log(APP_PER_AREA_MOD) + 
                              # Pesticide name random effect
                              (1 | PESTICIDE_MOD) +
                              # Catchment random effect
                              (1 | OPCAT_NAME) ,
                            family = Gamma(link = "log"),
                            data = gcmsSummary) 

# Checl R^2 values
MuMIn::r.squaredGLMM(gcmsModMixed)

# Check fixed effect
summary(gcmsModMixed)

# LCMS

# For every compound, aggregate all concentrations at every site by the median
lcmsSummary <- lcmsData %>%
  group_by(Sample_Site_ID, PESTICIDE_MOD, OPCAT_NAME, APP_PER_AREA_MOD) %>% 
  summarise(ConcentrationMedian = median(Concentration)) %>%
  ungroup()

# Plot relationship
ggplot(data = lcmsSummary,
       aes(x = log(ConcentrationMedian), 
           y = log(APP_PER_AREA_MOD))) +
  geom_point()

# Fit mixed model with concentration as response
lcmsModMixed <- lme4::glmer(Concentration ~ 
                              # Fixed log-transformed estimated application
                              log(APP_PER_AREA_MOD) + 
                              # Pesticide name random effect
                              (1 | PESTICIDE_MOD) +
                              # Catchment random effect
                              (1 | OPCAT_NAME) ,
                            family = Gamma(link = "log"),
                            data = lcmsData) 

# Check R^2 values
MuMIn::r.squaredGLMM(lcmsModMixed)

# Check fixed effect
summary(lcmsModMixed)



# Analysis 2


# Filter all records to a single value for each visit
visitData <- gcmsSummary %>%
  complete(Sample_Site_ID,PESTICIDE_MOD,
           fill = list(ConcentrationMedian=0))


ggplot(data = visitData,
       aes(x = log(ConcentrationMedian), 
           y = log(APP_PER_AREA_MOD))) +
  geom_point()


test <- glmmTMB::glmmTMB(formula = ConcentrationMedian ~ log(APP_PER_AREA_MOD) + (1|PESTICIDE_MOD),
                 zi = ~ (1|PESTICIDE_MOD),
                 family = Gamma(link = "log"),
                 data = visitData)


# Check R^2 values
MuMIn::r.squaredGLMM(test)

# Check fixed effect
summary(test)




cor(lcmsData$Concentration, lcmsData$APPLICATION_MOD)




test <- lm(log(ConcentrationMedian) ~ log(APP_PER_AREA_MOD),
           data = visitData)
summary(test)

test2 <- glm(Concentration ~ log(APPLICATION_MOD),
             family = Gamma(link = "log"),
             data = lcmsData[lcmsData$Compound_Name == "Acetamiprid",])
summary(test2)




summary(testModMixed)
hist(log(testData$MODELLED_APPLICATION))
hist(log(testData$MODELLED_CONCENTRATION))

testData[testData$Concentration == 3070, ]
testData$unit %>% unique()
