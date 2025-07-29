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

### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"
# If working locally: "../Data/Processed/Species/
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"

# SET PARAMETERS ---------------------------------------------------------------

# DATA FILES -------------------------------------------------------------------

# Load SDM fixed effect summaries
load(paste0(dataDir, "Species_effects.Rdata"))


# BOX PLOT  --------------------------------------------------------------------


# Box plot of connectivity model effect sizes

# Plot
 effects_df[, ] %>%

  ggplot() +
  geom_boxplot(aes(x = effect, y = mean   )) +
   #facet_wrap(~taxa) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #facet_wrap(~taxa) +
  #ylim(c(-8,8)) +
  coord_flip()

# Save
ggsave(filename = paste0("../Writing/Plots/", "Raw_connectivity_effect_sizes.png"),
       connEffPlot,
       dpi = 600,
       units = "px", width = 6000, height = 5000)

# DATA SUMMARIES --------------------------------------

# Proportion of ...

e.g. NROW(subset(meta_df,broadleafAssociation == "Y")) / NROW(meta_df) * 100


# EFFECTS CORRELATION MATRIX ---------------------------------------------------

# Change taxa names to standardised names
# (as tricky to use labeller within ggpairs() as you normally would with ggplot())
meta_df$labels <- str_replace_all(string = meta_df$taxa,
                                  pattern = taxaGroupLabels)

# Plot
effectPairs <- 
  # Subset to broadleaf species
  subset(meta_df, broadleafAssociation == "Y") %>%
  # Select columns we want for pairs plot
  select(c(labels, mean_coverBF, mean_connectivity, mean_BFconnINT)) %>%
  # Call ggpairs
  ggpairs(
    aes(colour = labels, alpha = 0.5),
    lower = list(continuous = wrap("cor", size = 3)),
    diag = list(continuous = wrap("densityDiag")),
    upper = list(
      combo = wrap("box_no_facet"),
      continuous = wrap("points")
    ),
    columnLabels = c(
      "Recording scheme",
      "Broadleaf cover effect",
      "Connectivity effect",
      "Interaction effect"
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
