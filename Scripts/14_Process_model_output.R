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
taxaGroups <- list.files(paste0(dataDir, "Model_outputs"))

# Create empty data frame to populate with all species level effects
effects_df <- data.frame()

# PROCESS OUTPUT ---------------------------------------------------------------

# Loop through each taxa group here
for (iTaxa in taxaGroups) {
  
  # List all summary files for iTaxa
  taxaOutput <- paste0(paste0(dataDir, "Model_outputs/", iTaxa)) %>%
  list.files(.,
             full.names = TRUE,
             recursive = TRUE)
  
  # Print taxa group as progress update
  print(iTaxa)
  
  # IMPORT FIXED EFFECTS -------------------------------------------------------
  
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

# CHECK SPECIES DUPLCIATES -----------------------------------------------------

# Find all duplicate species
duplicates <- effects_df %>%
  group_by(species, effect) %>%
  filter(n() > 1)

# Print any duplicates - should be 0
print(duplicates)

# RESTRUCTURE DATA FRAME -------------------------------------------------------

# Drop unneeded columns, then spread dataframe so that each effect mean, sd, and
# quantile is a separate column
effects_wide <-
  dplyr::select(effects_df, 
                -c("mode", "kld", "X0.5quant", )) %>%
  pivot_wider(names_from = effect,
              values_from = c("mean", "sd", "X0.025quant", "X0.975quant"))

# Remove any species with 0s for all values (modelling error)
effects_wide <- effects_wide %>%
  filter(if_any(where(is.numeric), ~ .x != 0))

# SAVE DATA FRAME --------------------------------------------------------------

# Save
save(effects_df, effects_wide,
     file = paste0(dataDir, "Species_effects.Rdata"))
