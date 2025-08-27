# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Descriptive statistics
#
# Script Description:

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)
library(GGally)

### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"
# If working locally: "../Data/Processed/Species/
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"
plotDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Plots/"

# SET PARAMETERS ---------------------------------------------------------------

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

# DATA FILES -------------------------------------------------------------------

# Loop through models that include and exclude wastewater
for (type in c("Wastewater", "NoWastewater")) {
  # Loop through the two species groups
  for (group in c("Schedule_2", "INNS")) {
    
    # Load SDM fixed effect summaries
    load(paste0(dataDir, "Species_effects/", type, "/", group, ".Rdata"))
    
    # assign(paste0(type, "_", group, "_effects_df" ), effects_df)
    # assign(paste0(type, "_", group, "_effects_wide" ), effects_wide)
    

# NoWastewater_INNS_effects_df
# Wastewater_INNS_effects_df
# NoWastewater_Schedule_2_effects_df
# Wastewater_Schedule_2_effects_df   
# 
# NoWastewater_INNS_effects_wide
# Wastewater_INNS_effects_wide
# NoWastewater_Schedule_2_effects_wide
# Wastewater_Schedule_2_effects_wide

effects_df <- Wastewater_Schedule_2_effects_df  
effects_wide <- Wastewater_Schedule_2_effects_wide

# BOX PLOT  --------------------------------------------------------------------

# Box plot of connectivity model effect sizes

# Plot
rawEffects <- effects_df[, ] %>%

  ggplot() +
  geom_boxplot(aes(x = effect, y = mean   )) +
   facet_wrap(~taxa) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #facet_wrap(~taxa) +
  #ylim(c(-8,8)) +
  coord_flip()

# Save
ggsave(filename = paste0(plotDir,
                         "Descriptive_analysis/",
                         type, "-", group,
                         "_raw_effect_sizes.png"),
       rawEffects,
       dpi = 600,
       units = "px", width = 6000, height = 5000)


# EFFECTS CORRELATION MATRIX ---------------------------------------------------

# Change taxa names to standardised names
# (as tricky to use labeller within ggpairs() as you normally would with ggplot())
effects_wide$labels <- str_replace_all(string = effects_wide$taxa,
                                       pattern = taxaGroupLabels)

# Plot
effectPairs <- effects_wide %>%
  # Select columns we want for pairs plot
  select(c(labels, mean_pesticideToxicity , mean_pesticideDiv , mean_NPK, mean_arable)) %>%
  # Call ggpairs
  ggpairs(
    aes(colour = labels, alpha = 0.5),
    #lower = list(continuous = wrap("cor", size = 3)),
    diag = list(continuous = wrap("densityDiag")),
    upper = list(
      combo = wrap("box_no_facet"),
      continuous = wrap("points")
    ),
    columnLabels = c(
      "Group",
      "Tox",
      "Div",
      "NPK",
      "Arable"
    ))

# Remove 'Taxa' histograms as redundant information and confusing

# Create plots list
plots = list()

# Iteratively collect the sub-plots we want
for (i in 1:4) {
  plots <- c(plots, lapply(2:effectPairs$ncol, function(j)
    getPlot(
      effectPairs, i = i, j = j
    )))
}

# Use ggmatrix to plot the subsetted sub-plots
effectPairsTrim <- ggmatrix(
  plots,
  nrow = 4,
  ncol = effectPairs$ncol - 1,
  xAxisLabels = effectPairs$xAxisLabels[2:effectPairs$ncol],
  yAxisLabels = effectPairs$yAxisLabels
)

# Save
ggsave(
  filename = paste0("../Writing/Plots/", "rawPairsPlot.png"),
  effectPairsTrim,
  dpi = 600,
  units = "px",
  width = 6000,
  height = 6000
)
