# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Model abundance of all species (species random effect)
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
library(cowplot)

# Set inla options
inla.setOption(num.threads = 2)
inla.setOption(inla.timeout = 0)

### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
plotDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Plots/"

### LOAD DATA ------------------------------------------------------------------

# Load Biosys data
invData <- readRDS(paste0(dataDir, "Processed/Biosys/invData_forModel.Rds"))

# SET PARAMETERS ---------------------------------------------------------------

linearLabels_NoW <- c('pesticideDiv' = "Pesticide diversity",
                      'pesticideToxicity' = "Pesticide combined toxicity",
                     'eutroph' = "Average Nitrogen and Potassium input",
                     'cattle' = "Cattle",
                     'pigs' = "Pigs",
                     'sheep' = "Sheep",
                     'poultry' = "Poultry",
                     'residential' = "Residential",
                     'woodland' = "Woodland",
                     'modification' = "Stream modification",
                     'quality' = "Habitat quality")

linearLabels_W <- c(linearLabels_NoW,
                    'wastewater' = "Wastewater")

randomLabels <- c( 'month' = "Month",
                   'year' = "Year")

# Priors for random effects
iidHyper <- list(prec = list(prior = "pc.prec",
                             param = c(100, 0.05)))
rwHyper <- list(prec = list(prior="pc.prec",
                            param=c(100, 0.05)))

### FILTER DATA ----------------------------------------------------------------

# Only keep needed columns for memory
invData <- invData %>%
  select(REPORTING_AREA,
         SITE_ID,
         SAMPLE_ID,
         TOTAL_ABUNDANCE,
         YEAR,
         MONTH_NUM,
         WEEK,
         TAXON,
         eutroph_PerArea_scaled,
         residential_PerArea_scaled,       
         woodland_PerArea_scaled,
         pesticideShannon_scaled,
         pesticideToxicLoad_PerArea_scaled,
         cattle_PerArea_scaled,            
         pigs_PerArea_scaled,
         sheep_PerArea_scaled,
         poultry_PerArea_scaled,          
         EDF_MEAN_scaled,
         HS_HMS_RSB_SubScore_scaled,
         HS_HQA_scaled,                    
         PC1_scaled,
         PC2_scaled,
         PC3_scaled,                       
         PC4_scaled,
         BASIN_F,
         CATCHMENT_F,                 
         WATER_BODY_F,
         GROUP)

# Filter to Schedule 2 species
invData <- filter(invData, GROUP == "Schedule 2")

### PROCESS TO RICHNESS FORMAT -------------------------------------------------

invData_SR <- add_count(invData, SAMPLE_ID, name = "numSpecies")

### PROCESS TO PSUEDO-ABSENCE FORMAT -------------------------------------------

# Create empty pseudo-absence table
invData_Abun_wZeroes <- NULL

# Loop through species here
  for (iSpecies in unique(invData$TAXON)) {
    
    # Create species abundance column with 0s
    speciesData <- invData %>%
      mutate(Abundance = ifelse(TAXON == iSpecies,
                                TOTAL_ABUNDANCE,
                                0)) %>%
      # Remove TOTAL_ABUNDANCE column as deprecated
      select(-TOTAL_ABUNDANCE)
    
    # Covert to unique column for each SAMPLE_ID
    speciesData <- speciesData %>%
      # For each SAMPLE_ID ...
      group_by(SAMPLE_ID) %>%
      # Extract max abundance value of iSpecies from speciesAbundance
      slice(which.max(Abundance)) %>%
      # Ungroup
      ungroup()
    
    # Add to dataframe with all species
    invData_Abun_wZeroes <- rbind(invData_Abun_wZeroes, speciesData)
    
  }

# Clear memory
rm(invData, speciesData)
gc()

### RUN RICHNESS MODELS --------------------------------------------------------
    
# SET MODEL COMPONENTS

# Richness model with wastewater
compsWastewater_SR <- numSpecies ~
  pesticideDiv(pesticideShannon_scaled, model = "linear") +
  pesticideToxicity(pesticideToxicLoad_PerArea_scaled, model = "linear") +
  eutroph(eutroph_PerArea_scaled, model = "linear") +
  cattle(cattle_PerArea_scaled, model = "linear") +
  pigs(pigs_PerArea_scaled, model = "linear") +
  sheep(sheep_PerArea_scaled, model = "linear") +
  poultry(poultry_PerArea_scaled, model = "linear") +
  residential(residential_PerArea_scaled, model = "linear") +
  woodland(woodland_PerArea_scaled, model = "linear") +
  modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
  quality(HS_HQA_scaled, model = "linear") +
  wastewater(EDF_MEAN_scaled, model = "linear") +
  PC1(PC1_scaled, model = "linear") +
  PC2(PC2_scaled, model = "linear") +
  PC3(PC3_scaled, model = "linear") +
  PC4(PC4_scaled, model = "linear") +
  month(
    main = MONTH_NUM,
    model = "rw1",
    scale.model = TRUE,
    hyper = rwHyper
  ) +
  year(YEAR,
       model = "rw1",
       scale.model = TRUE,
       hyper = rwHyper) +
  basin(BASIN_F, model = "iid", hyper = iidHyper) +
  catchment(CATCHMENT_F, model = "iid", hyper = iidHyper) +
  #wb(WATER_BODY_F, model = "iid", hyper = iidHyper) +
  Intercept(1)

# Richness model without wastewater
compsNoWastewater_SR <- numSpecies ~
  pesticideDiv(pesticideShannon_scaled, model = "linear") +
  pesticideToxicity(pesticideToxicLoad_PerArea_scaled, model = "linear") +
  eutroph(eutroph_PerArea_scaled, model = "linear") +
  cattle(cattle_PerArea_scaled, model = "linear") +
  pigs(pigs_PerArea_scaled, model = "linear") +
  sheep(sheep_PerArea_scaled, model = "linear") +
  poultry(poultry_PerArea_scaled, model = "linear") +
  residential(residential_PerArea_scaled, model = "linear") +
  woodland(woodland_PerArea_scaled, model = "linear") +
  modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
  quality(HS_HQA_scaled, model = "linear") +
  #wastewater(EDF_MEAN_scaled, model = "linear") +
  PC1(PC1_scaled, model = "linear") +
  PC2(PC2_scaled, model = "linear") +
  PC3(PC3_scaled, model = "linear") +
  PC4(PC4_scaled, model = "linear") +
  month(
    main = MONTH_NUM,
    model = "rw1",
    scale.model = TRUE,
    hyper = rwHyper
  ) +
  year(YEAR,
       model = "rw1",
       scale.model = TRUE,
       hyper = rwHyper) +
  basin(BASIN_F, model = "iid", hyper = iidHyper) +
  catchment(CATCHMENT_F, model = "iid", hyper = iidHyper) +
  #wb(WATER_BODY_F, model = "iid", hyper = iidHyper) +
  Intercept(1)

# RUN RICHNESS MODEL WITHOUT WASTEWATER

# Run model
modelNoWastewater_SR <- bru(
  components = compsNoWastewater_SR,
  family = "poisson",
  data = invData_SR,
  options = list(
    control.fixed = list(prec.intercept = 0.01),
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE),
    verbose = TRUE
  )
)

# RUN RICHNESS MODEL WITH WASTEWATER

# Run model
modelWastewater_SR <- bru(
  components = compsWastewater_SR,
  family = "poisson",
  data = invData_SR %>% filter(., !(is.na(EDF_MEAN_scaled))),
  options = list(
    control.fixed = list(prec.intercept = 0.01),
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE),
    verbose = TRUE
  )
)

gc()

### RUN ABUNDANCE MODELS -------------------------------------------------------

# Abundance model with wastewater
compsWastewater_Ab <- Abundance ~
  pesticideDiv(pesticideShannon_scaled, model = "linear") +
  pesticideToxicity(pesticideToxicLoad_PerArea_scaled, model = "linear") +
  eutroph(eutroph_PerArea_scaled, model = "linear") +
  cattle(cattle_PerArea_scaled, model = "linear") +
  pigs(pigs_PerArea_scaled, model = "linear") +
  sheep(sheep_PerArea_scaled, model = "linear") +
  poultry(poultry_PerArea_scaled, model = "linear") +
  residential(residential_PerArea_scaled, model = "linear") +
  woodland(woodland_PerArea_scaled, model = "linear") +
  modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
  quality(HS_HQA_scaled, model = "linear") +
  wastewater(EDF_MEAN_scaled, model = "linear") +
  PC1(PC1_scaled, model = "linear") +
  PC2(PC2_scaled, model = "linear") +
  PC3(PC3_scaled, model = "linear") +
  PC4(PC4_scaled, model = "linear") +
  month(
    main = MONTH_NUM,
    model = "rw1",
    scale.model = TRUE,
    hyper = rwHyper
  ) +
  year(YEAR,
       model = "rw1",
       scale.model = TRUE,
       hyper = rwHyper) +
  basin(BASIN_F, model = "iid", hyper = iidHyper) +
  catchment(CATCHMENT_F, model = "iid", hyper = iidHyper) +
  #wb(WATER_BODY_F, model = "iid", hyper = iidHyper) +
  Intercept(1)

# Abundance model without wastewater
compsNoWastewater_Ab <- Abundance ~
  pesticideDiv(pesticideShannon_scaled, model = "linear") +
  pesticideToxicity(pesticideToxicLoad_PerArea_scaled, model = "linear") +
  eutroph(eutroph_PerArea_scaled, model = "linear") +
  cattle(cattle_PerArea_scaled, model = "linear") +
  pigs(pigs_PerArea_scaled, model = "linear") +
  sheep(sheep_PerArea_scaled, model = "linear") +
  poultry(poultry_PerArea_scaled, model = "linear") +
  residential(residential_PerArea_scaled, model = "linear") +
  woodland(woodland_PerArea_scaled, model = "linear") +
  modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
  quality(HS_HQA_scaled, model = "linear") +
  #wastewater(EDF_MEAN_scaled, model = "linear") +
  PC1(PC1_scaled, model = "linear") +
  PC2(PC2_scaled, model = "linear") +
  PC3(PC3_scaled, model = "linear") +
  PC4(PC4_scaled, model = "linear") +
  month(
    main = MONTH_NUM,
    model = "rw1",
    scale.model = TRUE,
    hyper = rwHyper
  ) +
  year(YEAR,
       model = "rw1",
       scale.model = TRUE,
       hyper = rwHyper) +
  basin(BASIN_F, model = "iid", hyper = iidHyper) +
  catchment(CATCHMENT_F, model = "iid", hyper = iidHyper) +
  #wb(WATER_BODY_F, model = "iid", hyper = iidHyper) +
  Intercept(1)
    
# RUN ABUNDANCE MODEL WITHOUT WASTEWATER
    
# Run model
modelNoWastewater_Ab <- bru(
  components = compsNoWastewater_Ab,
  family = "zeroinflatednbinomial1",
  data = invData_Abun_wZeroes,
  options = list(
    control.fixed = list(prec.intercept = 0.01),
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE),
    verbose = TRUE
  )
)
gc()

# RUN ABUNDANCE MODEL WITH WASTEWATER

# Run model
modelWastewater_Ab <- bru(
  components = compsWastewater_Ab,
  family = "zeroinflatednbinomial1",
  data = invData_Abun_wZeroes %>% filter(., !(is.na(EDF_MEAN_scaled))),
  options = list(
    control.fixed = list(prec.intercept = 0.01),
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE),
    verbose = TRUE
  )
)
gc()
      
# PLOTS ------------------------------------------------------------------

# Loop through both models
models <- list(modelWastewater_SR = modelWastewater_SR,
               modelNoWastewater_SR = modelNoWastewater_SR,
               modelWastewater_Ab = modelWastewater_Ab,
               modelNoWastewater_Ab = modelNoWastewater_Ab)
      
for (modelName in names(models)) {
  
  # Get model
  model <- models[[modelName]]
  
  # Get model summary
  modelSummary <- summary(model)
  
  # Get linear effect labels
  if (modelName %in% c("modelNoWastewater_SR", "modelNoWastewater_Ab")) {
    linearLabels <- linearLabels_NoW
  }
  if (modelName %in% c("modelWastewater_SR", "modelWastewater_Ab")) {
    linearLabels <- linearLabels_W
  }
  
  # FIXED EFFECTS
  
  # Loop through variables and extract estimates
  for (i in names(linearLabels)) {
    # For covariate i, extract effect size
    effectSize <- modelSummary$inla$fixed[i, ] %>%
      t %>% # Transpose
      data.frame
    
    # Add covariate
    effectSize$Covariate <- i
    
    # If first covariate
    if (i == names(linearLabels)[1]) {
      # Create a new data frame
      effectSizeAll <- effectSize
      
    }  else {
      # Join data frames together
      effectSizeAll <- rbind(effectSizeAll, effectSize)
      
    }
  }
  
  # Plot fixed effects
  fixedEffPlot <- ggplot(
    effectSizeAll,
    aes(
      y = X0.5quant,
      x = Covariate,
      ymin = X0.025quant,
      ymax = X0.975quant,
      col = Covariate,
      fill = Covariate
    )
  ) +
    # Specify position here
    geom_linerange(linewidth = 4, colour = "lightblue") +
    ggtitle("Linear effects") +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(
      size = 2,
      shape = 21,
      colour = "white",
      fill = "black",
      stroke = 0.1
    ) +
    scale_x_discrete(
      name = "",
      limits = rev(names(linearLabels)),
      labels = as_labeller(linearLabels)
    ) +
    scale_y_continuous(name = "Effect size") +
    coord_flip() +
    theme_minimal() +
    guides(colour = "none") +
    theme(
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      legend.text = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, vjust = -0.5)
    )
  
  # RANDOM EFFECTS
  
  # Extract random effects from model, and exclude spatial
  randomEff_df <- model$summary.random
  
  # Add name of random effect to each dataframe in list
  randomEff_df <- imap(randomEff_df, ~ mutate(.x, randomEff = .y))
  
  # Unlist, then rename and select quantile columns
  randomEff_df <- do.call(rbind, randomEff_df) %>%
    rename("q0.025" = "0.025quant",
           "q0.5" = "0.5quant",
           "q0.975" = "0.975quant") %>%
    filter(randomEff %in% c("year", "month")) %>%
    select(ID, q0.025, q0.5, q0.975, randomEff)
  
  ### Plot
  
  randomEffPlot <- ggplot(randomEff_df) +
    
    # Random effect size
    geom_line(aes(x = as.numeric(ID), y = q0.5)) +
    geom_line(aes(x = as.numeric(ID), y = q0.025),
              lty = 2,
              alpha = .5) +
    geom_line(aes(x = as.numeric(ID), y = q0.975),
              lty = 2,
              alpha = .5) +
    
    # Thematics
    facet_wrap(~ randomEff,
               scale = 'free_x',
               labeller = as_labeller(randomLabels))  +
    ggtitle("Non-linear random effects") +
    theme(plot.title = element_text(hjust = 0.5),
          strip.text.x = element_text(size = 10)) +
    xlab("") +
    ylab("Count")
  
  # COMBINE PLOTS
  evalPlot <- plot_grid(fixedEffPlot,
                        randomEffPlot,
                        nrow = 2,
                        ncol = 1)
  
  # SAVE OUTPUT ------------------------------------------------------------

 # Create directory string for iSpecies
  if (modelName %in% c("modelNoWastewater_SR", "modelNoWastewater_Ab")) {
    # Create directory string for iSpecies
    iSpeciesDir <- paste0(
      dataDir,
      "Processed/Species/Model_outputs/NoWastewater/Schedule_2/AllSpecies")
  } else {
    iSpeciesDir <- paste0(
      dataDir,
      "Processed/Species/Model_outputs/Wastewater/Schedule_2/AllSpecies")
  }
  
  # Create directory
  dir.create(iSpeciesDir, recursive = TRUE, showWarnings = FALSE)
  
  # Save model summaries
  save(modelSummary, file = paste0(iSpeciesDir, "/AllSpecies.Rds"))
  ggsave(
    paste0(iSpeciesDir, "/AllSpecies.png"),
    evalPlot,
    width = 3000,
    height = 3000,
    units = "px",
    dpi = 400,
    limitsize = FALSE
  )
}
