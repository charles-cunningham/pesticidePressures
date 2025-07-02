# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Plot meta analysis
#
# Script Description:

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(brms)
library(tidyverse)
library(tidybayes)
library(ggridges)
library(rphylopic)
library(wesanderson)

### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"
# If working locally: "../Data/Processed/Species/
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"
plotDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Plots"

### RELOAD OBJECTS ---------------------------------------------------

# List of brms objects
brmsList <- c("pesticideDiv_brms", "pesticideToxicity_brms", "NPK_brms",
              "cattle_brms", "sheep_brms","pigs_brms", "poultry_brms",
              "wastewater_brms", "modification_brms", "quality_brms",
              "arable_brms", "urban_brms", "length_brms")

# Load brms objects
load(paste0(dataDir, "Species_brms.Rdata"))

# Load SDM fixed effect summaries
load(paste0(dataDir, "Species_effects.Rdata"))

# SUMMARISE EFFECT SIZES -------------------------------------------------------

# Summarise effect sizes by  woodland association and taxa
summary_taxa_df <- effects_wide %>%
  # Group by taxa group and effect category...
  group_by(taxa) %>%
  # ... and count total number
  summarise(nuSpecies = length(species)) %>%
  # Ungroup for ...
  ungroup %>%
  # ... proportion
  mutate(freq = prop.table(nuSpecies),
         taxa)

### TAXA SUMMARIES

# Loop through brms object list
for (i in brmsList) {
  
  # Assign brms object
  iBrms <- get(i)
  
  ### PROCESS BAYESIAN DRAWS 
  
  # Extract draws grouping by taxa
  studyDraws <- spread_draws(iBrms,
                             r_taxa[taxa, ],
                             b_Intercept) %>%
    mutate(taxa_mean = r_taxa + b_Intercept)
  
  # Extract draws (all pooled)
  pooledEffectDraws <- spread_draws(iBrms, b_Intercept) %>%
    rename(taxa_mean = b_Intercept) %>%
    mutate(taxa = "Pooled species")
  
  # Bind taxa draws and pooled draws together, and reorder
  forestData <- bind_rows(studyDraws,
                          pooledEffectDraws) %>%
    
    ungroup() %>%
    mutate(taxa = reorder(taxa, taxa_mean)) %>%
    mutate(taxa = relevel(taxa,
                          "Pooled species",
                          after = Inf))
  
  # Group by taxa and summarize by intercept
  forestDataSummary <- group_by(forestData, taxa) %>%
    mean_qi(taxa_mean)
  
  # Print taxa summaries, and assign out for plotting
  print( i )
  print( forestDataSummary )
  assign(paste0(i, "_draws"), forestData)
  assign(paste0(i, "_drawsSummary"), forestDataSummary)
}

# PLOT TAU ---------------------------------------------------------------------


# Loop through different meta analyses
for (i in brmsList) {
  
  # Extract draws, subset to sd, and process for plot
  tauPlot <- get(i) %>%
    as_draws_df %>%
    data.frame %>%
    select(starts_with("sd")) %>%
    gather(key, tau) %>%
    mutate(key = str_remove(key, "sd_") %>%
             str_remove(., "__Intercept")) %>%
    
    # Plot
    ggplot(aes(x = tau, fill = key)) +
    geom_density(color = "transparent", alpha = 2 / 3) +
    scale_fill_viridis_d(NULL, end = .85) +
    scale_y_continuous(NULL, breaks = NULL) +
    xlab(expression(tau)) +
    ggtitle(i) +
    theme(panel.grid = element_blank())
  
  # Save
  ggsave(filename = paste0(plotDir, "/Meta_analysis/Tau/Meta_", i, "tau.png"),
         tauPlot,
         dpi = 600,
         units = "px", width = 6000, height = 5000)
}

# PLOT POOLED ESTIMATES ------------------------------------------

# Loop through meta analysis subsets
for (i in brmsList) {
  
  # Assign brms objects
  iBrms <- get(i)
  iBrmsDraws <- get(paste0(i, "_draws"))
  iBrmsSummary <- get(paste0(i, "_drawsSummary"))
  
  ### PLOT 
  
  taxaSummaries <- ggplot(data = iBrmsDraws,
                          aes(taxa_mean,
                              taxa,
                              fill = taxa )) +
    
    # Add densities
    geom_density_ridges(scale = 0.95,
                        rel_min_height = 0.01) +
    geom_pointinterval(data = iBrmsSummary,
                       linewidth = 2,
                       aes(xmin = .lower, xmax = .upper)) +
    
    # Change colours and labels
    #scale_y_discrete(labels = taxaGroupLabels) +
    
    # Add vertical lines for pooled effect mean and CI, and 0
    geom_vline(xintercept = fixef(iBrms)[1, 1],
               color = "grey",
               linewidth = 1) +
    geom_vline( xintercept = fixef(iBrms)[1, 3:4],
                color = "grey",
                linetype = 2) +
    geom_vline(xintercept = 0,
               color = "black",
               linewidth = 1) +
    
    # Add taxon silouettes
    # geom_phylopic(data = data.frame(taxa = levels(iBrmsDraws$taxa )) %>%
    #                 left_join(., phylopicImages , by = "taxa"),
    #               inherit.aes = FALSE,
    #               aes(x = min(iBrmsDraws$taxa_mean) +
    #                     (max(iBrmsDraws$taxa_mean) - 
    #                        min(iBrmsDraws$taxa_mean)) / 100,
    #                   y = taxa,
    #                   img = svg ),
    #               size = 0.7,
    #               na.rm = TRUE) +
    
    # Add labels
    labs(x = "Effect size", # summary measure
         y = element_blank()) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16))
  
  # Save
  ggsave(filename = paste0(plotDir, "/Meta_analysis/Meta_",
                           i, "_taxaSummaries.png"),
         taxaSummaries,
         dpi = 600,
         units = "px",
         width = 5000,
         height = 5000)
}

