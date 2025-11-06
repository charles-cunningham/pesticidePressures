# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Model abundance of common species  using inlabru
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
inla.setOption(num.threads = 32)
inla.setOption(inla.timeout = 300) # 5 minutes

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
englandSmooth <- readRDS(paste0(dataDir, "Raw/Country_data/EnglandSmooth.Rds"))

# SET PARAMETERS ---------------------------------------------------------------

linearLabels_NoW <- c('pesticideDiv' = "Pesticide diversity",
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
                      'upstreamArea' = "Upstream area")


linearLabels_W <- c(linearLabels_NoW,
                    'wastewater' = "Wastewater")

randomLabels <- c( 'month' = "Month",
                   'year' = "Year",
                   'altitude' = "Altitude",
                   'slope' = "Slope",
                   'discharge' = "Discharge",
                   'ph' = "Alkalinity")

# Set minimum number of records to model - only commonly recorded species
minRecords <- 1000

# Schedule 2 species list
invDataS2 <- invData %>%
  as_tibble() %>%
  filter(GROUP == "Schedule 2") %>%
  distinct(TAXON) %>%
  .$TAXON

# INNS species list
invDataINNS <- invData %>%
  as_tibble() %>%
  filter(GROUP == "INNS") %>%
  distinct(TAXON) %>%
  .$TAXON

### Download BNG WKT string
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = paste0(dataDir, "bng.prj"))

bng <- sf::st_crs(paste0(dataDir, "bng.prj"))$wkt

### CREATE MESH ----------------------------------------------------------------

# Max edge is as a rule of thumb (range/3 to range/10)
maxEdge <- 50

# Create mesh
mesh <- inla.mesh.2d(boundary = englandSmooth,
                     max.edge =  maxEdge,
                     cutoff = maxEdge/2,
                     crs = gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

# Define spatial SPDE priors
space <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(1 * maxEdge, 0.5),
  prior.sigma = c(1, 0.5))

### MODEL SET UP FOR INDIVIDUAL SPECIES ----------------------------------------
# Loop through taxa then species to preserve ordering

# Loop through groups
for (iTaxa in unique(invData$TAXON_GROUP_NAME)) {
  
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
    
    # RUN MODEL ----------------------------------------------------------------
    
    # SET MODEL PARAMETERS

    # Priors for random effects
    iidHyper <- list(prec = list(prior = "pc.prec",
                                 param = c(100, 0.05)))
    rwHyper <- list(prec = list(prior="pc.prec",
                                param=c(100, 0.05)))
    
    # SET MODEL COMPONENTS
    
    # Model with wastewater
    compsWastewater <- speciesAbundance ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      #pesticideToxicity(pesticideToxicLoad_scaled, model = "linear") +
      #eutroph(ChemicalApplication_scaled, model = "linear") +
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
      # PC3(PC3, model = "linear") +
      # PC4(PC4, model = "linear") +
      # PC5(PC5, model = "linear") +
      # PC6(PC6, model = "linear") +
      # PC7(PC7, model = "linear") +
      # PC8(PC8, model = "linear") +
      month(
        MONTH_NUM,
        model = "rw1",
        scale.model = TRUE,
        hyper = rwHyper) +
      year(YEAR,
           model = "rw1",
           scale.model = TRUE,
           hyper = rwHyper) +
      altitude(ALTITUDE_GRP,
               model = "rw2",
               scale.model = TRUE,
               hyper = rwHyper) +
      slope(SLOPE_GRP,
            model = "rw2",
            scale.model = TRUE,
            hyper = rwHyper) +
      discharge(DISCHARGE_GRP,
                model = "rw2",
                scale.model = TRUE,
                hyper = rwHyper) +
      ph(ALKALINITY_GRP,
         model = "rw2",
         scale.model = TRUE,
         hyper = rwHyper) +
      basin(BASIN_F, model = "iid", hyper = iidHyper) +
      catchment(CATCHMENT_F, model = "iid", hyper = iidHyper) +
      wb(WATER_BODY_F, model = "iid", hyper = iidHyper) +
      # space(main = geometry,
      #           model = space) +
      Intercept(1)
    
    # Model without wastewater
    compsNoWastewater <- speciesAbundance ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      #pesticideToxicity(pesticideToxicLoad_scaled, model = "linear") +
      #eutroph(ChemicalApplication_scaled, model = "linear") +
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
      # PC3(PC3, model = "linear") +
      # PC4(PC4, model = "linear") +
      # PC5(PC5, model = "linear") +
      # PC6(PC6, model = "linear") +
      # PC7(PC7, model = "linear") +
      # PC8(PC8, model = "linear") +
      month(
        MONTH_NUM,
        model = "rw1",
        scale.model = TRUE,
        hyper = rwHyper) +
      year(YEAR,
           model = "rw1",
           scale.model = TRUE,
           hyper = rwHyper) +
      altitude(ALTITUDE_GRP,
               model = "rw2",
               scale.model = TRUE,
               hyper = rwHyper) +
      slope(SLOPE_GRP,
            model = "rw2",
            scale.model = TRUE,
            hyper = rwHyper) +
      discharge(DISCHARGE_GRP,
                model = "rw2",
                scale.model = TRUE,
                hyper = rwHyper) +
      ph(ALKALINITY_GRP,
         model = "rw2",
         scale.model = TRUE,
         hyper = rwHyper) +
      basin(BASIN_F, model = "iid", hyper = iidHyper) +
      catchment(CATCHMENT_F, model = "iid", hyper = iidHyper) +
      wb(WATER_BODY_F, model = "iid", hyper = iidHyper) +
      # space(main = geometry,
      #           model = space) +
      Intercept(1)
    
    # RUN MODEL WITH WASTEWATER
    
    # Remove previous model
    if (exists("modelWastewater")) {rm(modelWastewater)}
    
    # Add escape if model does not converge(
    try(

      modelWastewater <- bru(
        components = compsWastewater,
        family = "zeroinflatednbinomial1",
        data = speciesData %>% filter(., !(is.na(EDF_MEAN_scaled))),
        options = list(
          control.fixed = list(prec.intercept = 0.01),
          control.compute = list(waic = TRUE,
                                 dic = TRUE,
                                 cpo = TRUE),
          verbose = TRUE)
      )
    )
    
    # RUN MODEL WITHOUT WASTEWATER
    
    # Remove previous model
    if (exists("modelNoWastewater")) {rm("modelNoWastewater")}
    
    # Add escape if model does not converge(
    try(
      
      modelNoWastewater <- bru(
        components = compsNoWastewater,
        family = "zeroinflatednbinomial1",
        data = speciesData,
        options = list(
          control.fixed = list(prec.intercept = 0.01),
          control.compute = list(waic = TRUE,
                                 dic = TRUE,
                                 cpo = TRUE),
          verbose = TRUE)
      )
    )
    
    # Only plot and save if both models converge
    if (!is.null(summary(modelWastewater)$inla) & 
        !is.null(summary(modelNoWastewater)$inla)) {
      
      # Loop through both models
      for (modelName in c("modelWastewater", "modelNoWastewater")) {
        
        # Loop through both models
        models <- list(modelWastewater = modelWastewater,
                       modelNoWastewater = modelNoWastewater)
        
        # Get model
        model <- models[[modelName]]
        
        # Get model summary
        modelSummary <- summary(model)
        
        # Get linear effect labels
        if (modelName == "modelNoWastewater") { linearLabels <- linearLabels_NoW}
        if (modelName == "modelWastewater") { linearLabels <- linearLabels_W}
        
        # PLOTS ------------------------------------------------------------------
        
        # FIXED EFFECTS
        
        # Loop through variables and extract estimates
        for (i in names(linearLabels)) {
          
          # For covariate i, extract effect size
          effectSize <- modelSummary$inla$fixed[i,] %>%
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
          aes(y = X0.5quant,
              x = Covariate,
              ymin = X0.025quant,
              ymax = X0.975quant,
              col = Covariate,
              fill = Covariate )) +
          # Specify position here
          geom_linerange(linewidth = 4, colour = "lightblue") +
          ggtitle("Linear effects") +
          geom_hline(yintercept = 0, lty = 2) +
          geom_point(size = 2,
                     shape = 21,
                     colour = "white",
                     fill = "black",
                     stroke = 0.1) +
          scale_x_discrete(name = "",
                           limits = rev(names(linearLabels)),
                           labels = as_labeller(linearLabels)) +
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
          filter(!(randomEff %in% c("basin", "catchment", "wb", "space"))) %>%
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
          facet_wrap( ~ randomEff,
                      scale = 'free_x',
                      labeller = as_labeller(randomLabels))  +
          ggtitle("Non-linear random effects") +
          theme(plot.title = element_text(hjust = 0.5),
                strip.text.x = element_text(size = 10)) +
          xlab("") +
          ylab("Count")
        
        # COBINE PLOTS
        evalPlot <- plot_grid(fixedEffPlot, randomEffPlot,
                              nrow = 2, ncol = 1)
        
        # SAVE OUTPUT ------------------------------------------------------------
        
        # Set folder
        
        # If iSpecies is Schedule 2
        if (iSpecies %in% invDataS2) {
          group <- "Schedule_2"
        } else if (iSpecies %in% invDataINNS) {
          group <- "INNS"
        }
        
        # Create directory string for iSpecies
        iSpeciesDir <- paste0(
          dataDir,
          "Processed/Species/",
          "Model_outputs/",
          gsub("model", "", modelName),
          "/",
          group,
          "/",
          iTaxa
        )
        
        # Create directories for iTaxa if they don't exist
        lapply(paste0(iSpeciesDir, c("/ModelSummary", "/ModelPlots")),
               function(x) {
                 dir.create(x, recursive = TRUE, showWarnings = FALSE)
               })
        
        # Save model summaries
        save(modelSummary,
             file = paste0(iSpeciesDir,
                           "/ModelSummary/",
                           iSpecies,
                           ".Rds"))
        ggsave(paste0(iSpeciesDir,
                      "/ModelPlots/", iSpecies, ".png"),
               evalPlot,
               width = 3000, height = 3000, 
               units = "px", dpi = 400,
               limitsize = FALSE)
      }
    }
  }
}
