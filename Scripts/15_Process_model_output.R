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

# SET PARAMETERS ---------------------------------------------------------------

# Set taxa groups to analyse
taxaGroups <- list.files(paste0(dataDir, "Model_outputs/abWastewater/Schedule_2"))

# LOOP THROUGH TYPE AND GROUP --------------------------------------------------

# Loop through models that include and exclude wastewater
for (type in c("abWastewater", "abNoWastewater",
               "trendWastewater", "trendNoWastewater")) {
  
  # Loop through the two species groups
  for (group in c("Schedule_2")) {
    
    # Create output directory
    dir.create(paste0(dataDir, "Species_effects/", type), recursive = TRUE)
    
    # PROCESS OUTPUT -----------------------------------------------------------
    
    # Create empty data frame to populate with all species level effects for 
    # type and group
    effects_df <- data.frame()
    
    # Loop through each taxa group here
    for (iTaxa in taxaGroups) {
      
      # List all summary files for iTaxa
      taxaOutput <- paste0(paste0(dataDir,
                                  "Model_outputs/",
                                  type,
                                  "/",
                                  group,
                                  "/",
                                  iTaxa,
                                  "/ModelSummary")) %>%
        list.files(.,
                   full.names = TRUE,
                   recursive = TRUE)
  
      # Print taxa group as progress update
      paste0(type, "-", group, "-", iTaxa) %>% print()
  
      # IMPORT FIXED EFFECTS ---------------------------------------------------
  
      # If any species models within iTaxa...
      if (length(taxaOutput) > 0) {
      
        # For every species model in taxa group...
        for (i in 1:length(taxaOutput)) {
          
          # Load species 'i' summary file into global environment
          load(taxaOutput[i], envir = .GlobalEnv)
          
          # Assign fixed effects to dataframe
          iSpeciesEffects_df <- data.frame(modelSummary$inla$fixed)
        
          # Add fixed effects name column
          iSpeciesEffects_df$effect <- rownames(iSpeciesEffects_df)
          
          # Add in species name to fixed effects data frames
          iSpeciesEffects_df$species <- basename(taxaOutput[i]) %>%
            sub(".Rds",
                "",
                .)
          # Add taxa column
          iSpeciesEffects_df$taxa <- iTaxa
          
          # Drop redundant row names
          rownames(iSpeciesEffects_df) <- NULL
          
          # Bind iSpeciesEffects_df to the aggregated datadrame of all species
          effects_df <-
            bind_rows(effects_df, iSpeciesEffects_df)
          
        }
      }
    }
  
  # CHECK SPECIES DUPLCIATES ---------------------------------------------------
  
  # Find all duplicate species
  duplicates <- effects_df %>%
    group_by(species, effect) %>%
    filter(n() > 1)
  
  # Print any duplicates - should be 0
  print(duplicates)
  
  # RESTRUCTURE DATA FRAME -----------------------------------------------------
  
  # Drop unneeded columns, then spread dataframe so that each effect mean, sd, 
  # and quantile is a separate column
  effects_wide <-
    dplyr::select(effects_df,-c("mode", "kld", "X0.5quant",)) %>%
    pivot_wider(
      names_from = effect,
      values_from = c("mean", "sd", "X0.025quant", "X0.975quant")
    )
  
  # Remove any species with 0s for all values (modelling error)
  effects_wide <- effects_wide %>%
    filter(if_any(where(is.numeric), ~ .x != 0))
  
  # CHEMICAL AND ABUNDANCE RELATIONSHIPS --------------------------------------- 
  
  # Create discrete association categories
  effects_wide <- mutate( effects_wide,
    pestToxSig = case_when(
      X0.025quant_pestTox > 0 & X0.975quant_pestTox > 0 ~ "Pos",
      X0.025quant_pestTox < 0 &  X0.975quant_pestTox < 0 ~ "Neg",
      TRUE ~ "NS"),
    eutrophSig = case_when(
      X0.025quant_eutroph > 0 & X0.975quant_eutroph > 0 ~ "Pos",
      X0.025quant_eutroph < 0 & X0.975quant_eutroph < 0 ~ "Neg",
      TRUE ~ "NS"),
    pesticideDivSig = case_when(
      X0.025quant_pesticideDiv > 0 & X0.975quant_pesticideDiv > 0 ~ "Pos",
      X0.025quant_pesticideDiv < 0 & X0.975quant_pesticideDiv < 0 ~ "Neg",
      TRUE ~ "NS")
  )

  # SAVE DATA FRAME ------------------------------------------------------------
  
  # Save
  save(
    effects_df,
    effects_wide,
    file = paste0(dataDir, "Species_effects/", type, "/", group, ".Rdata")
  )
  
  }
}
