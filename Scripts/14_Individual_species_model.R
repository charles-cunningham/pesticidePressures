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

inla.setOption(num.threads = "16:1")

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
                      'N' = "Nitrogen",
                      'P' = "Phosporus",
                      'K' = "Potassium",
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

# Set minimum number of records to model - only commonly recorded species
minRecords <- 1000

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
                                 param = c(10, 0.05)))
    rwHyper <- list(prec = list(prior="pc.prec",
                                param=c(10, 0.05)))
    
    # SET MODEL COMPONENTS
    
    # Model with wastewater
    compsWastewater <- speciesAbundance ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      pesticideToxicity(pesticideToxicLoad_PerArea_scaled, model = "linear") +
      N(fertiliser_n_PerArea_scaled, model = "linear") +
      P(fertiliser_p_PerArea_scaled, model = "linear") +
      K(fertiliser_k_PerArea_scaled, model = "linear") +
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
      month(main = MONTH_NUM,
            model = "rw2",
            scale.model = TRUE,
            cyclic = TRUE,
            hyper = rwHyper) +
      year(YEAR,
           model = "rw1",
           scale.model = TRUE,
           hyper = rwHyper) +
      basin(REPORTING_AREA_NESTED, model = "iid", constr = TRUE, hyper = iidHyper) +
      catchment(CATCHMENT_NESTED, model = "iid", constr = TRUE, hyper = iidHyper) +
      #wb(WATER_BODY_NESTED, model = "iid", constr = TRUE, hyper = iidHyper) +
      Intercept(1)
    
    # Model with wastewater
    compsNoWastewater <- speciesAbundance ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      pesticideToxicity(pesticideToxicLoad_PerArea_scaled, model = "linear") +
      N(fertiliser_n_PerArea_scaled, model = "linear") +
      P(fertiliser_p_PerArea_scaled, model = "linear") +
      K(fertiliser_k_PerArea_scaled, model = "linear") +
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
      month(main = MONTH_NUM,
            model = "rw2",
            scale.model = TRUE,
            cyclic = TRUE,
            hyper = rwHyper) +
      year(YEAR,
           model = "rw1",
           scale.model = TRUE,
           hyper = rwHyper) +
      basin(REPORTING_AREA_NESTED, model = "iid", constr = TRUE, hyper = iidHyper) +
      catchment(CATCHMENT_NESTED, model = "iid", constr = TRUE, hyper = iidHyper) +
      #wb(WATER_BODY_NESTED, model = "iid", constr = TRUE, hyper = iidHyper) +
      Intercept(1)
    
    # RUN MODEL WITH WASTEWATER
    
    # Remove previous model
    if (exists("modelWastewater")) {rm(modelWastewater)}
    
    # Add escape if model does not converge(
    try(
      
      modelWastewater <- bru(
        components = compsWastewater,
        family = "zeroinflatednbinomial1",
        data = speciesData %>% filter(., !(is.na(EDF_MEAN))),
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
