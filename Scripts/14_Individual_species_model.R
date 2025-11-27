# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Model abundance of common species using inlabru
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
inla.setOption(num.threads = 16)
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

linearLabels_NoW <- c(
  'pesticideDiv' = "Pesticide diversity",
  'pestTox' = "Pesticide toxicity",
  'eutroph' = "Mean N and P application",
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

# Set minimum number of records to model - remove rare species
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

# ### Download BNG WKT string
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = paste0(dataDir, "bng.prj"))

bng <- sf::st_crs(paste0(dataDir, "bng.prj"))$wkt

### CREATE MESH ----------------------------------------------------------------

# Max edge is as a rule of thumb (range/3 to range/10)
maxEdge <- 20

# Create mesh
mesh <- fm_mesh_2d_inla(boundary = englandSmooth,
                        max.edge =  maxEdge,
                        cutoff = maxEdge/5,
                        min.angle = 26,
                        crs =  gsub( "units=m", "units=km",
                                     st_crs(bng)$proj4string ))

# Define spatial SPDE priors
spaceHyper <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(1 * maxEdge, 0.5),
  prior.sigma = c(1, 0.5))

# Create mesh dataframe for examining spatial field
mesh_df <- fm_pixels(mesh,
                     mask = st_transform(englandSmooth,
                                         crs = gsub( "units=m", "units=km",
                                                     st_crs(bng)$proj4string )))

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
    
    speciesTrendData <- speciesData %>%
      group_by(SITE_ID) %>%
      filter(n()>=10) %>%
      filter(any(speciesAbundance > 0)) %>%
      mutate(AbTrend = lm(speciesAbundance~YEAR)$coefficients[["YEAR"]]) %>%
      slice(which.max(AbTrend)) %>%
      # Ungroup
      ungroup()

    # RUN ABUNDANCE MODEL ------------------------------------------------------
    
    # SET MODEL PARAMETERS

    # Priors for random effects
    iidHyper <- list(prec = list(prior = "pc.prec",
                                 param = c(1, 0.5)))
    rwHyper <- list(prec = list(prior="pc.prec",
                                param=c(1, 0.5)))
    
    # SET MODEL COMPONENTS
    
    # Model with wastewater
    abWastewaterComps <- speciesAbundance ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      pestTox(pesticideToxicLoad_scaled, model = "linear") +
      eutroph(eutroph_scaled, model = "linear") +
      cattle(cattle_scaled, model = "linear") +
      pigs(pigs_scaled, model = "linear") +
      sheep(sheep_scaled, model = "linear") +
      poultry(poultry_scaled, model = "linear") +
      residential(residential_scaled, model = "linear") +
      woodland(woodland_scaled, model = "linear") +
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
        model = "rw2",
        scale.model = TRUE,
        cyclic = TRUE,
        hyper = rwHyper) +
      year(YEAR,
           model = "rw2",
           scale.model = TRUE,
           hyper = rwHyper) +
      #basin(BASIN_F, model = "iid", hyper = iidHyper_SR) +
      #catchment(CATCHMENT_F, model = "iid", hyper = iidHyper_SR) +
      #wb(WATER_BODY_F, model = "iid", hyper = iidHyper_SR)# +
       space(main = geometry,
             model = spaceHyper)
      
    # Model without wastewater
    abNoWastewaterComps <- speciesAbundance ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      pestTox(pesticideToxicLoad_scaled, model = "linear") +
      eutroph(eutroph_scaled, model = "linear") +
      cattle(cattle_scaled, model = "linear") +
      pigs(pigs_scaled, model = "linear") +
      sheep(sheep_scaled, model = "linear") +
      poultry(poultry_scaled, model = "linear") +
      residential(residential_scaled, model = "linear") +
      woodland(woodland_scaled, model = "linear") +
      modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
      quality(HS_HQA_scaled, model = "linear") +
      upstreamArea(totalArea_scaled, model = "linear") +
      #wastewater(EDF_MEAN_scaled, model = "linear") +
      PC1(PC1, model = "linear") +
      PC2(PC2, model = "linear") +
      PC3(PC3, model = "linear") +
      PC4(PC4, model = "linear") +
      month(
        MONTH_NUM,
        model = "rw2",
        cyclic = TRUE,
        scale.model = TRUE,
        hyper = rwHyper) +
      year(YEAR,
           model = "rw2",
           scale.model = TRUE,
           hyper = rwHyper) +
      #basin(BASIN_F, model = "iid", hyper = iidHyper_SR) +
      #catchment(CATCHMENT_F, model = "iid", hyper = iidHyper_SR) +
      #wb(WATER_BODY_F, model = "iid", hyper = iidHyper_SR) #+
       space(main = geometry,
             model = spaceHyper)
    
    # RUN MODEL WITH WASTEWATER
    
    # Remove previous model
    if (exists("abWastewater")) {rm(abWastewater)}
    
    # Add escape if model does not converge(
    try(

      abWastewater <- bru(
        components = abWastewaterComps,
        family = "zeroinflatednbinomial1",
        data = speciesData %>% filter(., !(is.na(EDF_MEAN_scaled))),
        options = list(
          control.inla=list(int.strategy = "eb"),
          control.compute = list(waic = TRUE,
                                 dic = TRUE,
                                 cpo = TRUE),
          verbose = TRUE)
      )
    )
    
    # RUN MODEL WITHOUT WASTEWATER
    
    # Remove previous model
    if (exists("abNoWastewater")) {rm("abNoWastewater")}
    
    # Add escape if model does not converge(
    try(
      
      abNoWastewater <- bru(
        components = abNoWastewaterComps,
        family = "zeroinflatednbinomial1",
        data = speciesData,
        options = list(
          control.inla=list(int.strategy = "eb"),
          control.compute = list(waic = TRUE,
                                 dic = TRUE,
                                 cpo = TRUE),
          verbose = TRUE)
      )
    )
    
    # RUN TREND MODEL ----------------------------------------------------------
    
    # SET MODEL COMPONENTS
    
    # Model with wastewater
    trendWastewaterComps <- AbTrend ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      pestTox(pesticideToxicLoad_scaled, model = "linear") +
      eutroph(eutroph_scaled, model = "linear") +
      cattle(cattle_scaled, model = "linear") +
      pigs(pigs_scaled, model = "linear") +
      sheep(sheep_scaled, model = "linear") +
      poultry(poultry_scaled, model = "linear") +
      residential(residential_scaled, model = "linear") +
      woodland(woodland_scaled, model = "linear") +
      modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
      quality(HS_HQA_scaled, model = "linear") +
      wastewater(EDF_MEAN_scaled, model = "linear") +
      upstreamArea(totalArea_scaled, model = "linear") +
      PC1(PC1, model = "linear") +
      PC2(PC2, model = "linear") +
      PC3(PC3, model = "linear") +
      PC4(PC4, model = "linear") +
      #basin(BASIN_F, model = "iid", hyper = iidHyper_SR) +
      #catchment(CATCHMENT_F, model = "iid", hyper = iidHyper_SR) +
      #wb(WATER_BODY_F, model = "iid", hyper = iidHyper_SR)# +
      space(main = geometry,
            model = spaceHyper)
    
    # Model without wastewater
    trendNoWastewaterComps <- AbTrend ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      pestTox(pesticideToxicLoad_scaled, model = "linear") +
      eutroph(eutroph_scaled, model = "linear") +
      cattle(cattle_scaled, model = "linear") +
      pigs(pigs_scaled, model = "linear") +
      sheep(sheep_scaled, model = "linear") +
      poultry(poultry_scaled, model = "linear") +
      residential(residential_scaled, model = "linear") +
      woodland(woodland_scaled, model = "linear") +
      modification(HS_HMS_RSB_SubScore_scaled, model = "linear") +
      quality(HS_HQA_scaled, model = "linear") +
      upstreamArea(totalArea_scaled, model = "linear") +
      #wastewater(EDF_MEAN_scaled, model = "linear") +
      PC1(PC1, model = "linear") +
      PC2(PC2, model = "linear") +
      PC3(PC3, model = "linear") +
      PC4(PC4, model = "linear") +
      #basin(BASIN_F, model = "iid", hyper = iidHyper_SR) +
      #catchment(CATCHMENT_F, model = "iid", hyper = iidHyper_SR) +
      #wb(WATER_BODY_F, model = "iid", hyper = iidHyper_SR) #+
      space(main = geometry,
            model = spaceHyper)
    
    # RUN MODEL WITH WASTEWATER
    
    # Remove previous model
    if (exists("trendWastewater")) {rm(trendWastewater)}
    
    # Add escape if model does not converge(
    try(
      
      trendWastewater <- bru(
        components = trendWastewaterComps,
        data = speciesTrendData %>% filter(., !(is.na(EDF_MEAN_scaled))),
        options = list(
          control.inla=list(int.strategy = "eb"),
          control.compute = list(waic = TRUE,
                                 dic = TRUE,
                                 cpo = TRUE),
          verbose = TRUE)
      )
    )
    
    # RUN MODEL WITHOUT WASTEWATER
    
    # Remove previous model
    if (exists("trendNoWastewater")) {rm("trendNoWastewater")}
    
    # Add escape if model does not converge(
    try(
      
      trendNoWastewater <- bru(
        components = trendNoWastewaterComps,
        data = speciesTrendData,
        options = list(
          control.inla=list(int.strategy = "eb"),
          control.compute = list(waic = TRUE,
                                 dic = TRUE,
                                 cpo = TRUE),
          verbose = TRUE)
      )
    )
    
    
    # Only plot and save if both models converge
    if (!is.null(summary(abWastewater)$inla) & 
        !is.null(summary(abNoWastewater)$inla) &
        !is.null(summary(trendWastewater)$inla) &
        !is.null(summary(trendNoWastewater)$inla)) {

      # Loop through both models
      models <- list(abWastewater = abWastewater,
                     abNoWastewater = abNoWastewater,
                     trendWastewater = trendWastewater,
                     trendNoWastewater = trendNoWastewater) 
      
      
      # Loop through both models
      for (modelName in c("abWastewater", "abNoWastewater",
                          "trendWastewater", "trendNoWastewater")) {
        
        # Get model
        model <- models[[modelName]]
        
        # Get model summary
        modelSummary <- summary(model)
        
        # Get linear effect labels
        if (modelName %in% c("abNoWastewater", 
                             "trendNoWastewater" )) { linearLabels <- linearLabels_NoW}
        if (modelName %in% c("abWastewater", 
                             "trendWastewater" )) { linearLabels <- linearLabels_W}
        
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
        
        if (modelName %in% c("abWastewater", "abNoWastewater")) {
          
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
        
         }
        
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
        
        # SPATIAL FIELD
        
        # Predict spatial field over domain
        pred_df <- predict(model, mesh_df, ~list(space = space))

        # Plot spatial field
        spatialEffPlot <- ggplot() +
          gg(pred_df$space["mean"], geom = "tile") +
          gg(st_transform(englandSmooth,
                          crs = gsub( "units=m", "units=km",
                                      st_crs(bng)$proj4string)),
             alpha = 0,
             col = "black",
             size = 1.5) +
          theme_void() +
          theme(legend.position = "bottom") +
          scale_fill_distiller(palette = 'RdYlBu', direction = 1,
                               limits = c(-1,1)*max(abs(pred_df$space$mean))) +
          labs(fill = "Spatial Field   ")
        
        # COMBINE PLOTS
        if (modelName %in% c("abWastewater", "abNoWastewater")) {
        
        evalPlot <- plot_grid(fixedEffPlot, randomEffPlot, spatialEffPlot,
                              nrow = 1, ncol = 3)
        
        } else {
          
          evalPlot <- plot_grid(fixedEffPlot, spatialEffPlot,
                                nrow = 1, ncol = 2)
        }
        
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
               width = 8000, height = 2000, 
               units = "px", dpi = 400,
               limitsize = FALSE)
      }
    }
  }
}
