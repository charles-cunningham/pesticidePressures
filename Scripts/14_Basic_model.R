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
library(cowplot)

### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
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

# SET PARAMETERS ---------------------------------------------------------------

linearEffLabels <- c('pesticideDiv' = "Pesticide diversity",
                     'pesticideToxicity' = "Pesticide combined toxicity",
                     'NPK' = "Total chemical input (NPK)",
                     'cattle' = "Cattle",
                     'pigs' = "Pigs",
                     'sheep' = "Sheep",
                     'poultry' = "Poultry",
                     'arable' = "Arable",
                     'grass' = "Intensive grassland",
                     'residential' = "Residential",
                     'modification' = "Stream modification",
                     'quality' = "Habitat quality")

randomEffLabels <- c( 'month' = "Month",
                      'year' = "Year")

# Set minimum number of records to model
minRecords <- 100

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

# FILTER VARIABLES -------------------------------------------------------------

# Filter temporally ( filter years after 2010)
invData <- invData %>%
  filter(YEAR > 2010)

# Rescale months and years
invData$MONTH_NUM <- invData$MONTH_NUM - (min(invData$MONTH_NUM) - 1)
invData$YEAR <- invData$YEAR - (min(invData$YEAR) - 1)

# Filter available  data
invData <- invData %>%
  # Remove rows with no upstream data
  filter(!(is.na(pesticideLoad))) %>%
  # Remove rows with no site data
  filter(!(is.na(HS_HMS_RSB_SubScore))) %>%
  filter(!(is.na(BIO_DEPTH)))%>%
  filter(!(is.na(BIO_SAND))) %>%
  filter(!(is.na(BIO_SILT_CLAY)))

### PROCESS TAXONOMY -----------------------------------------------------------

# Change species names to be file friendly
invData$TAXON_GROUP_NAME <- invData$TAXON_GROUP_NAME %>%
  # Remove "insect - " prefix
  gsub("insect - ", "", .) %>%
  # Remove spaces
  gsub(" ", "_", .)

# Change taxanomic group names to be file friendly
invData$TAXON <- invData$TAXON %>%
  # Remove spaces
  gsub(" ", "_", .) %>%
  # Remove slashes
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

### AGGREGATE VARIABLES --------------------------------------------------------

# NPK
invData$NPK <- invData$fertiliser_k + invData$fertiliser_n + invData$fertiliser_p

# Woodland
invData$woodland <- invData$Deciduous_woodland + invData$Coniferous_woodland

# Residential
invData$residential <- invData$Urban + invData$Suburban

# MODIFY UPSTREAM VARIABLES TO PER AREA VALUES----------------------------------

# Divide upstream variables by area (excluding diversity)
for(variable in c("NPK",
                  "Arable",
                  "residential",
                  "pesticideLoad",
                  "pesticideToxicLoad",
                  "cattle",
                  "pigs",
                  "sheep",
                  "poultry",
                  "woodland",
                  "Improved_grassland")) {
  
  # Create new scaled column name
  colName <- paste0(variable, "_PerArea")
  
  # Assign scaled variable to new column
  invData[, colName] <- invData[, variable] / invData[, "totalArea"]
  
}

### CORRELATION PLOTS ----------------------------------------------------------

# Create correlation data frame
corr_df <- invData %>%
  select(pesticideShannon,
         pesticideToxicLoad_PerArea,
         NPK_PerArea,
         Arable_PerArea,
         residential_PerArea,
         Improved_grassland_PerArea,
         woodland_PerArea,
         cattle_PerArea,
         pigs_PerArea,
         sheep_PerArea,
         poultry_PerArea,
         EDF_MEAN,
         HS_HMS_RSB_SubScore,
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

# WITH WASTEWATER

# Filter out missing wastewater values and create correlation object
corPredictors <- filter(corr_df, !(is.na(EDF_MEAN))) %>%
  cor(.) 

corrplot::corrplot(corPredictors,
                   type = "upper", order = "original", diag = FALSE,
                   method = "number", addCoef.col="white", tl.col = "black",
                   tl.srt = 45, tl.cex = 0.6)

# WITHOUT WASTEWATER

# Unselect wasterwater column and create correlation object
corPredictors <- select(corr_df, !(EDF_MEAN)) %>%
  cor(.) 

corrplot::corrplot(corPredictors,
                   type = "upper", order = "original", diag = FALSE,
                   method = "number", addCoef.col="white", tl.col = "black",
                   tl.srt = 45, tl.cex = 0.5)

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

# SCALE VARIABLES --------------------------------------------------------------

# List variables to be scaled
modelVariables <- c(
  # Upstream variables
  "NPK_PerArea",
  "Arable_PerArea",
  "residential_PerArea",
  "Improved_grassland_PerArea",
  "woodland_PerArea",
  "pesticideShannon",
  "pesticideLoad_PerArea",
  "pesticideToxicLoad_PerArea",
  "cattle_PerArea",
  "pigs_PerArea",
  "sheep_PerArea",
  "poultry_PerArea",
  # Site variables
  "EDF_MEAN",
  "HS_HMS_RSB_SubScore",
  "HS_HQA",
  "PC1",
  "PC2",
  "PC3",
  "PC4")

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
    
    # Priors for fixed effects
    fixedHyper <- list( mean = 0,
                        prec = 1 )
    
    # Priors for random effects
    iidHyper <- list(prec = list(prior = "pc.prec",
                                 param = c(0.5, 0.01)))
    rwHyper <- list(prec = list(prior="pc.prec",
                                param=c(0.5, 0.01)))
    
    # SET MODEL COMPONENTS
    
    # Model with wastewater
    compsWastewater <- speciesAbundance ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      pesticideToxicity(pesticideToxicLoad_PerArea_scaled, model = "linear") +
      NPK(NPK_PerArea_scaled, model = "linear") +
      cattle(cattle_PerArea_scaled, model = "linear") +
      pigs(pigs_PerArea_scaled, model = "linear") +
      sheep(sheep_PerArea_scaled, model = "linear") +
      poultry(poultry_PerArea_scaled, model = "linear") +
      arable(Arable_PerArea_scaled, model = "linear") +
      grass(Improved_grassland_PerArea_scaled, model = "linear") +
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
            cyclic = TRUE,
            hyper = rwHyper,
            scale.model = TRUE) +
      year(YEAR,
           model = "rw1",
           hyper = rwHyper,
           scale.model = TRUE) +
      basin(REPORTING_AREA, model = "iid", hyper = iidHyper) +
      catchment(CATCHMENT, model = "iid", hyper = iidHyper) +
      #wb(WATER_BODY, model = "iid", hyper = iidHyper) +
      Intercept(1)
    
    # Model with wastewater
    compsNoWastewater <- speciesAbundance ~
      pesticideDiv(pesticideShannon_scaled, model = "linear") +
      pesticideToxicity(pesticideToxicLoad_PerArea_scaled, model = "linear") +
      NPK(NPK_PerArea_scaled, model = "linear") +
      cattle(cattle_PerArea_scaled, model = "linear") +
      pigs(pigs_PerArea_scaled, model = "linear") +
      sheep(sheep_PerArea_scaled, model = "linear") +
      poultry(poultry_PerArea_scaled, model = "linear") +
      arable(Arable_PerArea_scaled, model = "linear") +
      grass(Improved_grassland_PerArea_scaled, model = "linear") +
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
            cyclic = TRUE,
            hyper = rwHyper,
            scale.model = TRUE) +
      year(YEAR,
           model = "rw1",
           hyper = rwHyper,
           scale.model = TRUE) +
      basin(REPORTING_AREA, model = "iid", hyper = iidHyper) +
      catchment(CATCHMENT, model = "iid", hyper = iidHyper) +
      #wb(WATER_BODY, model = "iid", hyper = iidHyper) +
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
          control.fixed = fixedHyper,
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
          control.fixed = fixedHyper,
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
        
        # Get model
        model <- get(modelName)
        
        # Model summary
        modelSummary <- summary(model)
        
        # PLOTS ------------------------------------------------------------------
        
        # FIXED EFFECTS
        
        # Loop through variables and extract estimates
        for (i in names(linearEffLabels)) {
          
          # For covariate i, extract effect size
          effectSize <- modelSummary$inla$fixed[i,] %>%
            t %>% # Transpose
            data.frame
          
          # Add covariate
          effectSize$Covariate <- i
          
          # If first covariate
          if (i == names(linearEffLabels)[1]) {
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
                           limits = rev(names(linearEffLabels)),
                           labels = as_labeller(linearEffLabels)) +
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
                      labeller = as_labeller(randomEffLabels))  +
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
