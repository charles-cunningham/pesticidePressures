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
# Needed to run on LINUX machine
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
library(corrplot)
library(GGally)

### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

# Create processed Biosys data folder
lapply(paste0(dataDir, "Processed/Species"), function(x) {
  if (!file.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
})

### LOAD DATA ------------------------------------------------------------------

# Load Biosys data
tempStoreForTesting <- readRDS(paste0(dataDir, "Processed/Biosys/invDataSpatialAll.Rds"))
invData <- tempStoreForTesting

# PROCESS DATA STRUCTURE -------------------------------------------------------

# Convert to tibble as non-spatially explicit model
invData <- as_tibble(invData)

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

# FILTER VARIABLES ----------------------------------------------------

# Filter temporally ( filter years between 2011 and 2020 inclusive)
invData <- invData %>%
  filter(YEAR >= 2010 & YEAR < 2020)

# Rescale months and years
invData$MONTH_NUM <- invData$MONTH_NUM - (min(invData$MONTH_NUM) - 1)
invData$YEAR <- invData$YEAR - (min(invData$YEAR) - 1)

# Filter available  data
invData <- invData %>%
  # Remove rows with no upstream data
  filter(!(is.na(pesticideLoad))) %>%
  # Remove rows with no site data
  filter(!(is.na(EDF_MEAN))) %>%
  filter(!(is.na(HS_HMS))) %>%
  filter(!(is.na(BIO_DEPTH))) 

### PROCESS TAXONOMY -----------------------------------------------------------

# Change species names to be file friendly
invData$TAXON <- invData$TAXON %>%
  # Remove slashes
  gsub(" ", "_", .) %>%
  # Remove
  gsub("/", "-", .)

# Schedule 2 species list
invDataS2 <- invData %>%
  filter(GROUP == "Schedule 2") %>%
  distinct(TAXON) %>%
  .$TAXON

# INNS species list
invDataINNS <- invData %>%
  filter(GROUP == "INNS") %>%
  distinct(TAXON) %>%
  .$TAXON

### CORRELATION PLOTS ----------------------------------------------------------

corr_df <- invData %>%
  select(pesticideShannon,
         pesticideToxicLoad,
         fertiliser_n,
         fertiliser_p,
         fertiliser_k,
         Arable,
         Urban,
         Improved_grassland,
         Deciduous_woodland,
         Coniferous_woodland,
         cattle,
         pigs,
         sheep,
         poultry,
         EDF_MEAN,
         HS_HMS,
         HS_HQA,
         BIO_ALTITUDE,
         BIO_SLOPE, 
         BIO_DISTANCE_FROM_SOURCE,
         BIO_DISCHARGE,
         BIO_WIDTH,
         BIO_DEPTH, 
         BIO_BOULDERS_COBBLES, 
         BIO_PEBBLES_GRAVEL, 
         BIO_SAND,
         BIO_SILT_CLAY,
         YEAR,
         MONTH_NUM)

corr_df <- corr_df %>%
  mutate(NPK = fertiliser_n + 
           fertiliser_p + 
           fertiliser_k) %>%
  select(-c(fertiliser_n, 
            fertiliser_p,
            fertiliser_k))

corPredictors <- cor(corr_df) 

corrplot::corrplot(corPredictors,
                   type = "upper", order = "original", diag = FALSE,
                   method = "number", addCoef.col="white", tl.col = "black",
                   tl.srt = 45)

# Plot
# effectPairs <- corr_df %>%
#   # Select columns we want for pairs plot
#   select(-c(fertiliser_n_scaled,
#             fertiliser_p_scaled,
#             fertiliser_k_scaled)) %>%
#   # Call ggpairs
#   ggpairs(
#     aes( alpha = 0.5),
#     lower = list(continuous = wrap("cor", size = 3)),
#     diag = list(continuous = wrap("densityDiag"))#,
#     # upper = list(
#     #   combo = wrap("box_no_facet"),
#     #   continuous = wrap("points")
#   )

### AGGREGATE VARIABLES --------------------------------------------------------

# NPK
invData$NPK <- invData$fertiliser_k + invData$fertiliser_n + invData$fertiliser_p

# Woodland
invData$woodland <- invData$Deciduous_woodland + invData$Coniferous_woodland

# SCALE VARIABLES --------------------------------------------------------------

# List variables to be scaled
modelVariables <- c(
  # Upstream variables
  "NPK",
  "Arable",
  "Urban",
  "Improved_grassland",
  "woodland",
  "pesticideShannon",
  "pesticideLoad",
  "pesticideToxicLoad",
  "cattle",
  "pigs",
  "sheep",
  "poultry",
  # Site variables
  "EDF_MEAN",
  "HS_HMS",
  "HS_HQA",
  "BIO_ALTITUDE",
  "BIO_SLOPE", 
  "BIO_DISTANCE_FROM_SOURCE",
  "BIO_DISCHARGE",
  "BIO_WIDTH",
  "BIO_DEPTH", 
  "BIO_BOULDERS_COBBLES", 
  "BIO_PEBBLES_GRAVEL", 
  "BIO_SAND",
  "BIO_SILT_CLAY")

# Create additional scaled column for each modelVariables
for(variable in modelVariables) {
  
  # Create new scaled column name
  colName <- paste0(variable, "_scaled")
  
  # Assign scaled variable to new column
  invData[, colName] <- scale(invData[[variable]])[,1]
}

# Convert categorical variables for random effects to factors
invData$WATER_BODY <- as.factor(invData$WATER_BODY)
invData$CATCHMENT <- as.factor(invData$CATCHMENT)
invData$REPORTING_AREA <- as.factor(invData$REPORTING_AREA)

### CONVERT SITE VARIABLES TO PCA ----------------------------------------------

sitePCA <- invData %>% 
  select(BIO_ALTITUDE,
         BIO_SLOPE,
         BIO_DISTANCE_FROM_SOURCE,
         BIO_DISCHARGE,
         BIO_WIDTH,
         BIO_DEPTH,
         BIO_BOULDERS_COBBLES,
         BIO_PEBBLES_GRAVEL,
         BIO_SAND,
         BIO_SILT_CLAY) %>%
  prcomp(.)

summary(sitePCA)

invData <- cbind(invData, sitePCA$x)

### RUN SPECIES-LEVEL MODELS ---------------------------------------------------

# Start taxa here
# Loop through taxa then species to preserve ordering
for (iTaxa in unique(invData$TAXON_GROUP_NAME)) {
  
# Set minimum number of records to model
minRecords <- 100

# Find species within taxa
taxaSpecies <- invData %>%
  
  # Filter to TAXON_GROUP
  filter(TAXON_GROUP_NAME == iTaxa) %>%
  
  # Filter only rows of species with over 100 records
  group_by(TAXON) %>% 
  filter(n() > minRecords) %>%
  
  # Get unique TAXON names
  .[["TAXON"]] %>%
  unique()
  
  # Loop through species here
  for (iSpecies in taxaSpecies) {

    # PROCESS TO PRESENCE-ABSENCE FORMAT

    # Create iSpecies abundance column with 0s
    speciesData <- invData %>%
      mutate(speciesAbundance = ifelse(TAXON == iSpecies,
                                       TOTAL_ABUNDANCE,
                                       0)) %>%
      # Remove TOTAL_ABUNDANCE column as deprecated
      select(-TOTAL_ABUNDANCE)
    
    # Covert to unique column for each SAMPLE_ID
    speciesData <- speciesData %>%
      # For each SAMPLE_ID ...
      group_by(SAMPLE_ID) %>%
      # Extract max abundance value of iSpecies from speciesAbundance
      slice(which.max(speciesAbundance)) %>%
      # Ungroup
      ungroup()
    
    # SET MODEL PARAMETERS
    
    # Priors for fixed effects
    fixedHyper <- list( mean.intercept = 0,
                        prec.intercept = 1,
                        mean = 0,
                        prec = 1 )
    
    # Priors for random effects
    iidHyper <- list(prec = list(prior = "pc.prec",
                                  param = c(0.5, 0.01)))
    rwHyper <- list(prec = list(prior="pc.prec",
                                     param=c(0.5, 0.01)))

    # SET MODEL COMPONENTS

    comps <- speciesAbundance ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      pesticideToxicity(pesticideToxicLoad_scaled, model = "linear") +
      NPK(NPK_scaled, model = "linear") +
      cattle(cattle_scaled, model = "linear") +
      pigs(pigs_scaled, model = "linear") +
      sheep(sheep_scaled, model = "linear") +
      poultry(poultry_scaled, model = "linear") +
      wastewater(EDF_MEAN_scaled, model = "linear") +
      modification(HS_HMS_scaled, model = "linear") +
      quality(HS_HQA_scaled, model = "linear") +
      arable(Arable_scaled, model = "linear") +
      pasture(Improved_grassland_scaled, model = "linear") +
      urban(Urban_scaled, model = "linear") +
      woodland(woodland_scaled, model = "linear") +
      PC1(PC1, model = "linear") +
      PC2(PC2, model = "linear") +
      PC3(PC3, model = "linear") +
      PC4(PC4, model = "linear") +
      PC5(PC5, model = "linear") +
      PC6(PC6, model = "linear") +
      month(main = MONTH_NUM,
            model = "rw2",
            cyclic = TRUE,
            hyper = rwHyper,
            scale.model = TRUE) +
      year(YEAR, model = "rw1",
           hyper = rwHyper,
           scale.model = TRUE) +
      basin(REPORTING_AREA, model = "iid", hyper = iidHyper) +
      catchment(CATCHMENT, model = "iid", hyper = iidHyper) +
      wb(WATER_BODY, model = "iid", hyper = iidHyper) +
      Intercept(1)

    # RUN MODEL
    
    # Remove previous model
    rm(model)
    
    # Add escape if model does not converge(
    try(
    
    model <- bru(
      components = comps,
      family = "zeroinflatednbinomial1",
      data = speciesData,
      options = list(
        control.fixed = fixedHyper,
        control.compute = list(waic = TRUE,
                               dic = TRUE,
                               cpo = TRUE),
        verbose = TRUE)
    )
    )

    # SAVE OUTPUT
    
    # If model converges...
    if (exists("model")) {
    
    # Model summary
    modelSummary <- summary(model) ; modelSummary

    # Set folder
    
    # If iTaxa is Schedule 2
    if (iSpecies %in% invDataS2) {
      group <- "Schedule_2"
    } else if (iSpecies %in% invDataINNS) {
      group <- "INNS"
    }
    
    # Create direcotry for iTaxa
    speciesDir <- paste0(
      dataDir,
      "Processed/Species/",
      "Model_outputs/",
      group,
      "/",
      iTaxa)
    
    if (!file.exists(speciesDir)) {
      dir.create(speciesDir, recursive = TRUE)
    }
                  
    # Save model summary
    save(modelSummary,
         file = paste0(
           speciesDir,
           "/",
           iSpecies,
           ".Rds")
         )
  }
}
}

### PLOT

# LABELS

linearEffLabels <- c('pesticideDiv' = "Pesticide diversity",
                     #'pesticideLoad' = "Pesticide total application",
                     'pesticideToxicity' = "Pesticide combined toxicity",
                     'N' = "Nitrogen",
                     'P' = "Phosphorus",
                     'K' = "Potassium",
                     'upstream' = "Total upstream area",
                     'arable' = "Upstream arable area",
                     'urban' = "Upstream urban area")

randomEffLabels <- c('alitude' = "alitude",
                     'slope' = "slope",
                     'discharge' = "discharge",
                     'width' = "width",
                     'depth' = "depth" ,
                     'boulders' = "boulders",
                     'pebbles' = "pebbles",
                     'sand' = "sand",
                     'silt' = "silt" ,
                     'alkalinity' = "alkalinity",
                     'month' = "month")

# RANDOM

# Extract random effects from model, and exclude spatial
randomEff_df <- model$summary.random

# Add name of random effect to each dataframe in list
randomEff_df <- imap(randomEff_df, ~mutate(.x, randomEff = .y))

# Unlist, then rename and select quantile columns
randomEff_df <- do.call(rbind, randomEff_df)%>%
  rename("q0.025" = "0.025quant",
         "q0.5" = "0.5quant",
         "q0.975" = "0.975quant") %>%
  dplyr::select(ID, q0.025, q0.5, q0.975, randomEff)

### Plot

randomEffPlot <- ggplot(randomEff_df) +

  # Random effect size
  geom_line(aes(x = as.numeric(ID), y = q0.5)) +
  geom_line(aes(x = as.numeric(ID), y = q0.025), lty = 2, alpha = .5) +
  geom_line(aes(x = as.numeric(ID), y = q0.975), lty = 2, alpha = .5) +

  # Thematics
  facet_wrap(~ randomEff, scale = 'free_x') +
            # labeller = as_labeller(randomEffLabels))  +
  ggtitle("Non-linear random effects") +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10)) +
  xlab("") +
  ylab("Count")
randomEffPlot

# LINEAR

# Loop through covariates and extract estimates
for (i in names(linearEffLabels)) {

  # For covariate i, extract effect size
  effectSize <- modelSummary$inla$fixed[i,] %>%
    t %>% # Transpose
    data.frame

  # Add covariate
  effectSize$Covariate <- i

  # If first covariate
  if( i == names(linearEffLabels)[1]) {

    # Create a new data frame
    effectSizeAll <- effectSize

  }  else {

    # Join data frames together
    effectSizeAll <- rbind(effectSizeAll, effectSize)

  }
}

# Plot fixed effects
fixedEffPlot <- ggplot(effectSizeAll,
                       aes(y = X0.5quant, x = Covariate,
                           ymin = X0.025quant, ymax=X0.975quant,
                           col = Covariate, fill = Covariate)) +
  #specify position here
  geom_linerange(linewidth=4, colour = "lightblue") +
  ggtitle("Linear effects") +
  geom_hline(yintercept=0, lty=2) +
  geom_point(size=2, shape=21, colour="white", fill = "black", stroke = 0.1) +
  scale_x_discrete(name="",
                   limits = rev(names(linearEffLabels)),
                   labels = as_labeller(linearEffLabels)) +
  scale_y_continuous(name="Effect size") +
  coord_flip() +
  theme_minimal() +
  guides(colour = "none") +
  theme(axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, vjust = -0.5))

fixedEffPlot



#### TESTING
# Create presence/absence column for ZAP model
speciesData <- mutate(speciesData,
                      speciesPresent = (speciesAbundance > 0) * 1L)

comps <- ~
  pesticideDiv(pesticideShannon_scaled, model = "linear") +
  pesticideToxicity(pesticideToxicLoad_scaled, model = "linear") +
  NPK(fertiliser_n_scaled +
        fertiliser_p_scaled +
        fertiliser_k_scaled, model = "linear") +
  cattle(cattle_scaled, model = "linear") +
  pigs(pigs_scaled, model = "linear") +
  sheep(sheep_scaled, model = "linear") +
  poultry(poultry_scaled, model = "linear") +
  wastewater(EDF_MEAN_scaled, model = "linear") +
  modification(HS_HMS_scaled, model = "linear") +
  quality(HS_HQA_scaled, model = "linear") +
  arable(Arable_scaled, model = "linear") +
  urban(Urban_scaled, model = "linear") +
  altitude(BIO_ALTITUDE_scaled_grp,
           model = "rw2",
           scale.model = TRUE,
           hyper = riverPrior) +
  discharge(BIO_DISCHARGE_scaled,
            model = "rw2",
            scale.model = TRUE,
            hyper = riverPrior) +
  month(main = MONTH_NUM,
        model = "seasonal",
        season.length = 12,
        hyper = timePrior,
        scale.model = TRUE) +
  year(YEAR,
       model = "rw2",
       scale.model = TRUE,
       hyper = timePrior) +
  basin(REPORTING_AREA, model = "iid", hyper = locPrior) +
  catchment(CATCHMENT, model = "iid", hyper = locPrior) +
  wb(WATER_BODY, model = "iid", hyper = locPrior) +
  Intercept_present(1) +
  Intercept_count(1)


poisson_obs <- 
  bru_obs(
    family = "zeroinflatedpoisson1",
    data = speciesData,
    formula = speciesAbundance ~ pesticideDiv +
      pesticideToxicity +
      NPK +
      cattle +
      pigs +
      sheep + 
      poultry +
      wastewater +
      modification +
      quality +
      arable +
      urban +
      altitude +
      month +
      discharge +
      basin + 
      catchment +
      wb +
      Intercept_count,
    control.family = list(hyper = list(theta = list(initial = -20,
                                                    fixed = TRUE)
    )))

present_obs <- bru_obs(
  family = "binomial",
  data = speciesData,
  formula = speciesPresent ~ 
    year +
    Intercept_present
)

# RUN MODEL

model <- bru(
  comps,
  present_obs,
  truncated_poisson_obs,
  options = list(
    control.compute = list(
      waic = TRUE,
      dic = TRUE,
      cpo = TRUE
    ),
    verbose = TRUE
  ))

# Create summary
modelSummary <- summary(model) ; modelSummary


