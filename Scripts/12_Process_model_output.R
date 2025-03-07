# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Process model output
#
# Script Description: 


### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(INLA)
library(inlabru)


### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"
# If working locally: "../Data/Processed/Species/
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"


# SET PARAMETERS ------------------------------------

# Set taxa groups to analyse
taxaGroups <- list.files(paste0(dataDir, "Model_outputs"))

# PROCESS OUTPUT ---------------------------------------------------------------

  # Loop through each taxa group here
  for (iTaxa in taxaGroups) {
    
    # Set taxa group
    taxaGroup <- iTaxa; print(taxaGroup)
    
    # Set taxa group folder
    taxaDir <- paste0(paste0(dataDir, "Model_outputs/", iTaxa))
    
    # List summary files
    taxaOutput <- list.files(taxaDir,
                                  full.names = TRUE,
                                  recursive = TRUE)
    
    # IMPORT FIXED EFFECTS ---------------------------------
    
    # Create empty list to populate with model fixed effect summaries
    speciesEffects <- list()
    
    # For every converged species model in taxa group...
    for (i in 1:length(taxaOutput)) {

      # Create list just for that species (will become sub-list)
      iSpeciesEffects <- list()
      
      # Assign species name from file path to list
      # (Remove file name and directory path)
      iSpeciesEffects$Species <-
        sub(".Rds",
            "",
            taxaOutput[i]) %>%
        sub(paste0('.*', taxaGroup, '/'), 
            '',
            .)
      
      # Load species 'i' summary file into global environment
      load(taxaOutput[i], envir = .GlobalEnv)
      
      # Assign fixed effects to list
      iSpeciesEffects$FixedEffects <- data.frame(modelSummary$inla$fixed)
      
      # Save list of species i as a sub-list of speciesEffects
      speciesEffects[[i]] <- iSpeciesEffects
      
    }
    
    # Add in species name to fixed effects data frames,
    # and effect from row names
    for (i in 1:length(speciesEffects)) {
      speciesEffects[[i]][[2]]$species <- speciesEffects[[i]][[1]]
      speciesEffects[[i]][[2]]$effect <- rownames(speciesEffects[[i]][[2]])
    }
    
    # Extract fixed effects data frames and bind together
    speciesEffects_df <- map(speciesEffects, 2, .default = NA) %>% 
      bind_rows
    
    # Drop redundant row names
    rownames(speciesEffects_df) <- NULL
    
    # Add taxa column
    speciesEffects_df$taxa <- taxaGroup
    
    # Add to allSpEff_df dataframe
    allSpEff_df <- bind_rows(allSpEff_df, speciesEffects_df)
    
  }
  
  # CHECK DUPLCIATES -----------------------------------------------
  
  # Find all duplicate rows
  duplicates <- allSpEff_df %>% 
    group_by(species) %>% 
    filter(n()>
             ( NROW(speciesEffects_df) / length(allSummaryFiles) ) 
    )# 8 effect sizes per species
  
  # Print if any duplicates - should be 0
  print(duplicates)
  
  # RESTRUCTURE DATA FRAME ------------------------------
  
  # Drop unneeded columns, then spread dataframe so that each effect mean, sd, and
  # quantile is a separate column
  SI_df <-
    dplyr::select( allSpEff_df, -c("mode", "kld", "X0.5quant", )) %>%
    pivot_wider(
      names_from = effect,
      values_from = c( "mean", "sd", "X0.025quant", "X0.975quant" ))
  
  # Remove any species with 0s for all mean and sd values (modelling error)
  SI_df <- SI_df %>%
    filter(rowSums(abs(SI_df[-c(1:2)])) > 0)
  

  # SAVE DATA FRAME -----------------------------------------------
  
  save(SI_df,
       file = paste0("../Data/Species_data/SDM_fixed_effect_summaries_",
                     SI_dir,
                     ".RData"))
}
