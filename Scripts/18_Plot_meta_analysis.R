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

### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"
# If working locally: "../Data/Processed/Species/
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"
plotDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Plots/"

# SET PARAMETERS ------------------------------------

# Set taxa groups to analyse
taxaGroups <- c("annelid",                "crustacean",
                "flatworm_(Turbellaria)", "alderfly_(Megaloptera)",
                "beetle_(Coleoptera)",    "caddis_fly_(Trichoptera)",
                "dragonfly_(Odonata)",    "mayfly_(Ephemeroptera)",
                "stonefly_(Plecoptera)",  "true_bug_(Hemiptera)",
                "true_fly_(Diptera)",     "mollusc")

# Set taxa group labels
taxaGroupLabels <- c( "annelid" = "Leeches",
                      "crustacean" = "Crustaceans",
                      "flatworm_(Turbellaria)" = "Flatworms",
                      "alderfly_(Megaloptera)" = "Alderflies",
                      "beetle_(Coleoptera)" = "Beetles",
                      "caddis_fly_(Trichoptera)" = "Caddisflies",
                      "dragonfly_(Odonata)" = "Dragonflies",
                      "mayfly_(Ephemeroptera)" = "Mayflies",
                      "stonefly_(Plecoptera)" = "Stoneflies",
                      "true_bug_(Hemiptera)" = "True bugs",
                      "true_fly_(Diptera)" = "True flies",
                      "mollusc" = "Molluscs")


# Set phylopic images (choose uuid manually)
phylopicImages <- data.frame(taxa = taxaGroups,
                             uuid = c(
                               "91973387-c0fb-4193-a28d-17fe4284c4aa",
                               "639a8581-cf15-4ab0-a81a-d4319f7078d7",
                               "22f8d558-6763-426d-a651-f4087090fc41",
                               "3c5af8a0-ea99-4df9-a344-88585844b24e",
                               "74e85bf0-e705-4512-9e53-be48e39251eb",
                               "e1307c88-3e8f-4ba8-9f93-751df3deb739",
                               "92a2b163-34dc-4fef-9aef-f348d8b3af6c",
                               "777bbd68-7924-4a4d-a5dd-452111d2a823",
                               "a0cbfbcb-71e3-4b21-bbfd-95a6ce1f246c",
                               "2d5699f4-f8e9-48d0-8668-eedff9a87343",
                               "558127d5-4e77-41b3-8d4b-6eb7d7deeb9c",
                               "c5835123-e2d3-4c20-9e7a-f7b6528bbf8e"
                             )) %>%
  mutate(svg = lapply(uuid, get_phylopic)) # Use uuids to get image objects

# Rotate images
phylopicImages$svg[[2]] <- rotate_phylopic(img = phylopicImages$svg[[2]], angle = 180)
phylopicImages$svg[[3]] <- rotate_phylopic(img = phylopicImages$svg[[3]], angle = 90)
phylopicImages$svg[[4]] <- flip_phylopic(img = phylopicImages$svg[[4]])
phylopicImages$svg[[5]] <- rotate_phylopic(img = phylopicImages$svg[[5]], angle = 90)
phylopicImages$svg[[9]] <- flip_phylopic(img = phylopicImages$svg[[9]])
phylopicImages$svg[[10]] <- rotate_phylopic(img = phylopicImages$svg[[10]], angle = 90)
phylopicImages$svg[[11]] <- rotate_phylopic(img = phylopicImages$svg[[11]], angle = 90)
phylopicImages$svg[[12]] <- rotate_phylopic(img = phylopicImages$svg[[12]], angle = 90)

# Get attribution using get_attribution():
#get_attribution(uuid = "5aeaf558-3c48-4173-83b4-dbf2846f8d75")

### LOOP THROUGH WASTEWATER AND NO WASTEWATER ----------------------------------

for (type in c("Wastewater", "NoWastewater")) {

### RELOAD OBJECTS ---------------------------------------------------

# List of brms objects
brmsList <- c("pesticideDiv_brms", "pesticideToxicity_brms", "NPK_brms",
              "cattle_brms", "sheep_brms","pigs_brms", "poultry_brms",
              "modification_brms", "quality_brms", "arable_brms",
              "residential_brms","grass_brms", "woodland_brms")

if (type == "Wastewater" ) {
  
  brmsList <- c(brmsList, "wastewater_brms")
  
}

# Load brms objects
load(paste0(dataDir, type, "_species_brms.Rdata"))

# Load SDM fixed effect summaries
load(paste0(dataDir, "Species_effects/", type, "/Schedule_2.Rdata"))

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
  allDraws <- bind_rows(studyDraws,
                          pooledEffectDraws) %>%
    
    ungroup() %>%
    mutate(taxa = reorder(taxa, taxa_mean)) %>%
    mutate(taxa = relevel(taxa,
                          "Pooled species",
                          after = Inf))
  
  # Group by taxa and summarize by intercept
  allDrawsSummary <- group_by(allDraws, taxa) %>%
    mean_qi(taxa_mean)
  
  # Print taxa summaries, and assign out for plotting
  print( i )
  print( allDraws )
  assign(paste0(i, "_draws"), allDraws)
  assign(paste0(i, "_drawsSummary"), allDrawsSummary)
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
  ggsave(filename = paste0(plotDir, "Meta_analysis/",
                           type, "/Tau/Meta_", i, "tau.png"),
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
    geom_density_ridges(scale = 0.9,
                        rel_min_height = 0.1) +
    geom_pointinterval(data = iBrmsSummary,
                       linewidth = 2,
                       aes(xmin = .lower, xmax = .upper)) +
    
    # Change colours and labels
    scale_y_discrete(labels = taxaGroupLabels) +
    
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
    geom_phylopic(data = data.frame(taxa = levels(iBrmsDraws$taxa )) %>%
                    left_join(., phylopicImages , by = "taxa"),
                  inherit.aes = FALSE,
                  aes(x = min(iBrmsDraws$taxa_mean) +
                        (max(iBrmsDraws$taxa_mean) -
                           min(iBrmsDraws$taxa_mean)) / 80,
                      y = taxa,
                      img = svg ),
                  width = (max(iBrmsDraws$taxa_mean) -
                             min(iBrmsDraws$taxa_mean)) / 10,
                  na.rm = TRUE) +

    # Add labels
    labs(x = "Effect size", # summary measure
         y = element_blank()) +
    ggtitle(gsub("_brms", "", i)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          plot.title = element_text(size = 24))
  
  # Save
  ggsave(filename = paste0(plotDir, "Meta_analysis/", type, "/Effects/Meta_",
                           i, "_taxaSummaries.png"),
         taxaSummaries,
         dpi = 600,
         units = "px",
         width = 5000,
         height = 5000)
}
}
