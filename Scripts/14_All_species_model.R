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

# England boundary
#englandSmooth <- readRDS(paste0(dataDir, "Raw/Country_data/EnglandSmooth.Rds"))

# SET PARAMETERS ---------------------------------------------------------------

linearLabels_NoW <- c(
  'pesticideDiv' = "Pesticide diversity",
  'chemApp' = "Pesticides and NP",
  'lessPest' = "More NP than pesticide",
  'cattle' = "Cattle",
  'pigs' = "Pigs",
  'sheep' = "Sheep",
  'poultry' = "Poultry",
  'residential' = "Residential",
  'woodland' = "Woodland",
  'modification' = "Stream modification",
  'quality' = "Habitat quality",
  'upstreamArea' = "Upstream area"
)

linearLabels_W <- c(linearLabels_NoW,
                    'wastewater' = "Wastewater")

randomLabels <- c( 'month' = "Month",
                   'year' = "Year")

# Priors for random effects
iidHyper_SR <- list(prec = list(prior = "pc.prec",
                             param = c(100, 0.05)))
rwHyper_SR <- list(prec = list(prior="pc.prec",
                            param=c(100, 0.05)))
iidHyper_Ab <- list(prec = list(prior = "pc.prec",
                                param = c(100, 0.05)))
rwHyper_Ab <- list(prec = list(prior="pc.prec",
                               param=c(100, 0.05)))

### Download BNG WKT string
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = paste0(dataDir, "bng.prj"))

bng <- sf::st_crs(paste0(dataDir, "bng.prj"))$wkt

### FILTER DATA ----------------------------------------------------------------

# Only keep needed columns for memory
invData <- invData %>%
  select(
    pesticideShannon_scaled,
    chemicalApp,
    lessPesticide,
    residential,
    woodland,
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
    MONTH_NUM,
    REPORTING_AREA,
    SITE_ID,
    SAMPLE_ID,
    TOTAL_ABUNDANCE,
    TAXON,
    GROUP,
    BASIN_F,
    CATCHMENT_F,
    WATER_BODY_F
    )

# Filter to Schedule 2 species
invData <- filter(invData, GROUP == "Schedule 2")

### PROCESS TO RICHNESS FORMAT -------------------------------------------------

invData_SR <- add_count(invData, SAMPLE_ID, name = "numSpecies")

### PROCESS TO PSUEDO-ABSENCE FORMAT -------------------------------------------

# Create empty pseudo-absence table
invData_Ab_wZeroes <- NULL

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
    invData_Ab_wZeroes <- rbind(invData_Ab_wZeroes, speciesData)
    
  }

# Clear memory
rm(invData, speciesData)
gc()

### CREATE MESH ----------------------------------------------------------------

# # Max edge is as a rule of thumb (range/3 to range/10)
# maxEdge <- 50
# 
# # Create mesh
# mesh <- inla.mesh.2d(boundary = englandSmooth,
#                      max.edge =  maxEdge,
#                      cutoff = maxEdge/2,
#                      crs = gsub( "units=m", "units=km", st_crs(bng)$proj4string ))
# 
# # Define spatial SPDE priors
# spaceHyper <- inla.spde2.pcmatern(
#   mesh,
#   prior.range = c(1 * maxEdge, 0.5),
#   prior.sigma = c(1, 0.5))

### RUN RICHNESS MODELS --------------------------------------------------------
    
# SET MODEL COMPONENTS

# Richness model with wastewater
compsWastewater_SR <- numSpecies ~
  pesticideDiv(pesticideShannon_scaled, model = "linear") +
  chemApp(chemicalApp, model = "linear") +
  lessPest(lessPesticide, model = "linear") +
  cattle(cattle_scaled, model = "linear") +
  pigs(pigs_scaled, model = "linear") +
  sheep(sheep_scaled, model = "linear") +
  poultry(poultry_scaled, model = "linear") +
  residential(residential, model = "linear") +
  woodland(woodland, model = "linear") +
  modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
  quality(HS_HQA_scaled, model = "linear") +
  wastewater(EDF_MEAN_scaled, model = "linear") +
  upstreamArea(totalArea_scaled, model = "linear") +
  PC1(PC1, model = "linear") +
  PC2(PC2, model = "linear") +
  PC3(PC3, model = "linear") +
  PC4(PC4, model = "linear") +
  month(
    MONTH_NUM,
    model = "rw1",
    scale.model = TRUE,
    hyper = rwHyper_SR) +
  year(YEAR,
       model = "rw1",
       scale.model = TRUE,
       hyper = rwHyper_SR) +
  basin(BASIN_F, model = "iid", hyper = iidHyper_SR) +
  catchment(CATCHMENT_F, model = "iid", hyper = iidHyper_SR) +
  #wb(WATER_BODY_F, model = "iid", hyper = iidHyper_SR) +
  # space(main = geometry,
  #       model = spaceHyper) +
  Intercept(1)

# Richness model without wastewater
compsNoWastewater_SR <- numSpecies ~
  pesticideDiv(pesticideShannon_scaled, model = "linear") +
  chemApp(chemicalApp, model = "linear") +
  lessPest(lessPesticide, model = "linear") +
  cattle(cattle_scaled, model = "linear") +
  pigs(pigs_scaled, model = "linear") +
  sheep(sheep_scaled, model = "linear") +
  poultry(poultry_scaled, model = "linear") +
  residential(residential, model = "linear") +
  woodland(woodland, model = "linear") +
  modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
  quality(HS_HQA_scaled, model = "linear") +
  upstreamArea(totalArea_scaled, model = "linear") +
  #wastewater(EDF_MEAN_scaled, model = "linear") +
  PC1(PC1, model = "linear") +
  PC2(PC2, model = "linear") +
  PC3(PC3, model = "linear") +
  PC4(PC4, model = "linear") +
  month(
    main = MONTH_NUM,
    model = "rw1",
    scale.model = TRUE,
    hyper = rwHyper_SR) +
  year(YEAR,
       model = "rw1",
       scale.model = TRUE,
       hyper = rwHyper_SR) +
  basin(BASIN_F, model = "iid", hyper = iidHyper_SR) +
  catchment(CATCHMENT_F, model = "iid", hyper = iidHyper_SR) +
  #wb(WATER_BODY_F, model = "iid", hyper = iidHyper_SR) +
  # space(main = geometry,
  #            model = spaceHyper) +
  Intercept(1)

# RUN RICHNESS MODEL WITHOUT WASTEWATER

# Run model
modelNoWastewater_SR <- bru(
  components = compsNoWastewater_SR,
  family = "poisson",
  data = invData_SR,
  options = list(
    control.fixed = list(prec.intercept = 0.01),
    control.inla=list(int.strategy = "eb"),
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
    control.inla=list(int.strategy = "eb"),
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE),
    verbose = TRUE
  )
)

gc()

### RUN ABUNDANCE MODELS -------------------------------------------------------

# Abundance model with wastewater
compsWastewater_Ab <- Abundance ~
  pesticideDiv(pesticideShannon_scaled, model = "linear") +
  chemApp(chemicalApp, model = "linear") +
  lessPest(lessPesticide, model = "linear") +
  cattle(cattle_scaled, model = "linear") +
  pigs(pigs_scaled, model = "linear") +
  sheep(sheep_scaled, model = "linear") +
  poultry(poultry_scaled, model = "linear") +
  residential(residential, model = "linear") +
  woodland(woodland, model = "linear") +
  modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
  quality(HS_HQA_scaled, model = "linear") +
  wastewater(EDF_MEAN_scaled, model = "linear") +
  upstreamArea(totalArea_scaled, model = "linear") +
  PC1(PC1, model = "linear") +
  PC2(PC2, model = "linear") +
  PC3(PC3, model = "linear") +
  PC4(PC4, model = "linear") +
  month(
    main = MONTH_NUM,
    model = "rw1",
    scale.model = TRUE,
    hyper = rwHyper_Ab) +
  year(YEAR,
       model = "rw1",
       scale.model = TRUE,
       hyper = rwHyper_Ab) +
  basin(BASIN_F, model = "iid", hyper = iidHyper_Ab) +
  catchment(CATCHMENT_F, model = "iid", hyper = iidHyper_Ab) +
  #wb(WATER_BODY_F, model = "iid", hyper = iidHyper_Ab) +
  # space(main = geometry,
   #       model = spaceHyper) +
  species(TAXON,  model = "iid", hyper = iidHyper_Ab) +
  Intercept(1)

# Abundance model without wastewater
compsNoWastewater_Ab <- Abundance ~
  pesticideDiv(pesticideShannon_scaled, model = "linear") +
  chemApp(chemicalApp, model = "linear") +
  lessPest(lessPesticide, model = "linear") +
  cattle(cattle_scaled, model = "linear") +
  pigs(pigs_scaled, model = "linear") +
  sheep(sheep_scaled, model = "linear") +
  poultry(poultry_scaled, model = "linear") +
  residential(residential, model = "linear") +
  woodland(woodland, model = "linear") +
  modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
  quality(HS_HQA_scaled, model = "linear") +
  upstreamArea(totalArea_scaled, model = "linear") +
  #wastewater(EDF_MEAN_scaled, model = "linear") +
  PC1(PC1, model = "linear") +
  PC2(PC2, model = "linear") +
  PC3(PC3, model = "linear") +
  PC4(PC4, model = "linear") +
  month(
    main = MONTH_NUM,
    model = "rw1",
    scale.model = TRUE,
    hyper = rwHyper_Ab) +
  year(YEAR,
       model = "rw1",
       scale.model = TRUE,
       hyper = rwHyper_Ab) +
  basin(BASIN_F, model = "iid", hyper = iidHyper_Ab) +
  catchment(CATCHMENT_F, model = "iid", hyper = iidHyper_Ab) +
  #wb(WATER_BODY_F, model = "iid", hyper = iidHyper_Ab) +
  # space(main = geometry,
  #       model = spaceHyper) +
  species(TAXON,  model = "iid", hyper = iidHyper_Ab) +
  Intercept(1)
    
# RUN ABUNDANCE MODEL WITHOUT WASTEWATER
    
# Run model
modelNoWastewater_Ab <- bru(
  components = compsNoWastewater_Ab,
  family = "zeroinflatednbinomial1",
  data = invData_Ab_wZeroes,
  options = list(
    control.fixed = list(prec.intercept = 0.01),
    control.inla=list(cmin=0,
                      int.strategy = "eb"),
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
  data = invData_Ab_wZeroes %>% filter(., !(is.na(EDF_MEAN_scaled))),
  options = list(
    control.fixed = list(prec.intercept = 0.01),
    control.inla=list(cmin=0,
                      int.strategy = "eb"),
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
    filter(!(randomEff %in% c("basin", "catchment", "wb", "species"))) %>%
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
      "Processed/Species/Model_outputs/NoWastewater/Schedule_2/AllSpecies/")
  } else {
    iSpeciesDir <- paste0(
      dataDir,
      "Processed/Species/Model_outputs/Wastewater/Schedule_2/AllSpecies/")
  }
  
  # Create directory
  dir.create(iSpeciesDir, recursive = TRUE, showWarnings = FALSE)
  
  # Save model summaries
  save(modelSummary, file = paste0(iSpeciesDir, modelName, ".Rds"))
  ggsave(
    paste0(iSpeciesDir, modelName, ".png"),
    evalPlot,
    width = 3000,
    height = 3000,
    units = "px",
    dpi = 400,
    limitsize = FALSE
  )
}
