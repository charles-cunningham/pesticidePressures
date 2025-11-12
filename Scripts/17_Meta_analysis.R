# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Meta analysis
#
# Script Description:

### LOAD LIBRARIES -------------------------------------------------------------

# Run once for R 4.4.2 to get the brms package working (belt and braces)
# Install RTools
#remove.packages(c("StanHeaders", "rstan", "brms"))
# if (file.exists(".RData")) file.remove(".RData")
# # RESTART R
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# options(mc.cores = parallel::detectCores())
# example(stan_model, package = "rstan", run.dontrun = TRUE) # This checks rstan and the C++ compiler are correctly installed
# # RESTART R
# install.packages("brms")
# rstan_options(auto_write = TRUE)

# Load packages
library(brms)
library(tidyverse)
library(tidybayes)
library(ggridges)

### DIRECTORY MANAGEMENT -------------------------------------------------------
# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"
# If working locally: "../Data/Processed/Species/
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/Processed/Species/"

# SET PARAMETERS ---------------------------------------------------------------

# DATA FILES -------------------------------------------------------------------

# Loop through models that include and exclude wastewater
for (type in c("Wastewater", "NoWastewater")) {

  # Load SDM fixed effect summaries
  load(paste0(dataDir, "Species_effects/", type, "/Schedule_2.Rdata"))

### META ANALYSIS --------------------------------------
# In this stage we want to check overall data patterns,
# partitioning of variance between taxa and species levels etc.
# Run all brm functions first in a single section. It takes a while
# so can save the outputs so only need to run once.

# RUN BAYESIAN MODELS ---------------------------------

# Fit Bayesian meta analysis model
# N.B. Species effect nested within taxa effect
# (i.e. taxa + taxa:species random effects)

pesticideDiv_brms <- brm(data = effects_wide,
                 family = gaussian,
                 mean_pesticideDiv | se(sd_pesticideDiv) ~
                   1 + (1 | taxa) + (1 | taxa:species),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(cauchy(0, 1), class = sd)),
                 iter = 10000,
                 warmup = 5000,
                 # control=list(adapt_delta = 0.99,
                 #              stepsize = 0.01,
                 #             max_treedepth = 15),
                 cores = 4,
                 chains = 4)

chemApp_brms <- brm(data = effects_wide,
                         family = gaussian,
                         mean_chemApp | se(sd_chemApp) ~
                           1 + (1 | taxa) + (1 | taxa:species),
                         prior = c(prior(normal(0, 1), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                         iter = 10000,
                         warmup = 5000,
                         # control=list(adapt_delta = 0.99,
                         #              stepsize = 0.01,
                         #              max_treedepth = 15),
                         cores = 4,
                         chains = 4)

lessPest_brms <- brm(data = effects_wide,
                family = gaussian,
                mean_lessPest | se(sd_lessPest) ~
                  1 + (1 | taxa) + (1 | taxa:species),
                prior = c(prior(normal(0, 1), class = Intercept),
                          prior(cauchy(0, 1), class = sd)),                
                iter = 10000,
                warmup = 5000,
                # control=list(adapt_delta = 0.99,
                #              stepsize = 0.01,
                #              max_treedepth = 15),
                         cores = 4,
                         chains = 4)

cattle_brms <- brm(data = effects_wide,
                         family = gaussian,
                         mean_cattle | se(sd_cattle) ~
                           1 + (1 | taxa) + (1 | taxa:species),
                         prior = c(prior(normal(0, 1), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                   iter = 10000,
                   warmup = 5000,
                   # control=list(adapt_delta = 0.99,
                   #              stepsize = 0.01,
                   #              max_treedepth = 15),
                         cores = 4,
                         chains = 4)

pigs_brms <- brm(data = effects_wide,
                         family = gaussian,
                         mean_pigs | se(sd_pigs) ~
                           1 + (1 | taxa) + (1 | taxa:species),
                         prior = c(prior(normal(0, 1), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                 iter = 10000,
                 warmup = 5000,
                 # control=list(adapt_delta = 0.99,
                 #              stepsize = 0.01,
                 #              max_treedepth = 15),
                         cores = 4,
                         chains = 4)

sheep_brms <- brm(data = effects_wide,
                         family = gaussian,
                         mean_sheep | se(sd_sheep) ~
                           1 + (1 | taxa) + (1 | taxa:species),
                         prior = c(prior(normal(0, 1), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                  iter = 10000,
                  warmup = 5000,
                  # control=list(adapt_delta = 0.99,
                  #              stepsize = 0.01,
                  #              max_treedepth = 15),
                         cores = 4,
                         chains = 4)

poultry_brms <- brm(data = effects_wide,
                         family = gaussian,
                         mean_poultry | se(sd_poultry) ~
                           1 + (1 | taxa) + (1 | taxa:species),
                         prior = c(prior(normal(0, 1), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                    iter = 10000,
                    warmup = 5000,
                    # control=list(adapt_delta = 0.99,
                    #              stepsize = 0.01,
                    #              max_treedepth = 15),
                         cores = 4,
                         chains = 4)

modification_brms <- brm(data = effects_wide,
                         family = gaussian,
                         mean_modification | se(sd_modification) ~
                           1 + (1 | taxa) + (1 | taxa:species),
                         prior = c(prior(normal(0, 1), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                         iter = 10000,
                         warmup = 5000,
                         # control=list(adapt_delta = 0.99,
                         #              stepsize = 0.01,
                         #              max_treedepth = 15),
                         cores = 4,
                         chains = 4)

quality_brms <- brm(data = effects_wide,
                         family = gaussian,
                         mean_quality | se(sd_quality) ~
                           1 + (1 | taxa) + (1 | taxa:species),
                         prior = c(prior(normal(0, 1), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                    iter = 10000,
                    warmup = 5000,
                    # control=list(adapt_delta = 0.99,
                    #              stepsize = 0.01,
                    #              max_treedepth = 15),
                         cores = 4,
                         chains = 4)

residential_brms <- brm(data = effects_wide,
                         family = gaussian,
                         mean_residential | se(sd_residential) ~
                           1 + (1 | taxa) + (1 | taxa:species),
                         prior = c(prior(normal(0, 1), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                  iter = 10000,
                  warmup = 5000,
                  # control=list(adapt_delta = 0.99,
                  #              stepsize = 0.01,
                  #              max_treedepth = 15),
                         cores = 4,
                         chains = 4)

woodland_brms <- brm(data = effects_wide,
                  family = gaussian,
                  mean_woodland | se(sd_woodland) ~
                    1 + (1 | taxa) + (1 | taxa:species),
                  prior = c(prior(normal(0, 1), class = Intercept),
                            prior(cauchy(0, 1), class = sd)),
                  iter = 10000,
                  warmup = 5000,
                  # control=list(adapt_delta = 0.99,
                  #              stepsize = 0.01,
                  #              max_treedepth = 15),
                  cores = 4,
                  chains = 4)

upstreamArea_brms <- brm(data = effects_wide,
                     family = gaussian,
                     mean_upstreamArea | se(sd_upstreamArea) ~
                       1 + (1 | taxa) + (1 | taxa:species),
                     prior = c(prior(normal(0, 1), class = Intercept),
                               prior(cauchy(0, 1), class = sd)),
                     iter = 10000,
                     warmup = 5000,
                     # control=list(adapt_delta = 0.99,
                     #              stepsize = 0.01,
                     #              max_treedepth = 15),
                     cores = 4,
                     chains = 4)

if(type == "Wastewater") {
  
  wastewater_brms <- brm(data = effects_wide,
                       family = gaussian,
                       mean_wastewater | se(sd_wastewater) ~
                         1 + (1 | taxa) + (1 | taxa:species),
                       prior = c(prior(normal(0, 1), class = Intercept),
                                 prior(cauchy(0, 1), class = sd)),
                       iter = 10000,
                       warmup = 5000,
                       # control=list(adapt_delta = 0.99,
                       #              stepsize = 0.01,
                       #              max_treedepth = 15),
                       cores = 4,
                       chains = 4)
}



### Save

# List of brms objects
brmsList <- c("pesticideDiv_brms", "chemApp_brms", "lessPest_brms",
              "cattle_brms", "sheep_brms", "pigs_brms", "poultry_brms",
              "wastewater_brms", "modification_brms", "quality_brms",
              "residential_brms", "upstreamArea_brms", "woodland_brms")

# Save brms objects
save(list = brmsList,
     file = paste0(dataDir, type, "_species_brms.Rdata"))

}
