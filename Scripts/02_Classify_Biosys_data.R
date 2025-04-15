# HEADER -----------------------------------------------------------------------
#
# Author: Charles Cunningham
# 
# Script Name: Create schedule 2 and INNS groups for biosys data
#
# Script Description: Identify and filter invData from script 01 to only include
# species which we are interested in and have sufficienct records to model, i.e.
# Schedule 2 and INNS species. Also apply other filtering steps on sample types
# Many thanks to Dr Martin Wilkes for code snippets for this script.

### LOAD LIBRARIES -------------------------------------------------------------

# Load packages
library(tidyverse)

### DATA MANAGEMENT ------------------------------------------------------------

# Set data directory
# If working on Databricks: "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"
# If working locally: "../Data/"
dataDir <- "/dbfs/mnt/lab/unrestricted/charles.cunningham@defra.gov.uk/Pesticides/Data/"

# Read Biosys data from script 01
invData <- readRDS(file = paste0(dataDir, "Raw/Biosys/invData_all.Rds"))

### PRODUCE SPECIES LISTS ------------------------------------------------------

# SCHEDULE 2

# Information obtained from The Environmental Targets (Biodiversity) 
# (England) Regulations 2023, available here:
# https://www.legislation.gov.uk/uksi/2023/91/made

# Pre-processing requires copying the scientific names of 'Freshwater 
# invertebrates' from the 'SCHEDULE 2 Species for the targets relating to the
# abundance of species' list into a .csv file.
# The .csv should be set up such that there is a single column called "Species",
# each row contains the scientific name of a single species (235 rows), and 
# called ' Schedule_2_fw_list.csv' and saved in 'Raw/Biosys/' directory

# Read in schedule 2 species
s2 <- paste0(dataDir, "Raw/Biosys/Schedule_2_fw_list.csv") %>%
  read.csv()

# INVASIVE NON-NATIVE SPECIES

# Information obtained from the UK Biodiversity Indicators 2024. Pressure from 
# invasive species 2024, available here:
# https://hub.jncc.gov.uk/assets/647caed5-93d0-4dc0-92bf-13d231a37dda

# Pre-processing requires downloading the 'UK-BDI-2024-invasive-species.xlsx'
# file using the 'UKBI 2024. Pressure from invasive species (Datafile, Excel)' 
# link within the above page.
# The sheet within this file named 'Worksheet 2: List of species' should then 
# be saved as a .csv file and named 'INNS_species_list.csv'

# Read in INNS species
inns <- paste0(dataDir, "Raw/Biosys/INNS_species_list.csv") %>%
  read.csv(., skip = 4)

### FILTER DATA ----------------------------------------------------------------

# Filter on sample type
invData <- invData %>%
  # ...to remove species which are not invertebrates
  filter(TAXON_TYPE  == "Other Macroinvertebrates") %>%
  #  ...to rivers only
  filter(WATERBODY_TYPE == "WBRV") %>%
  # ...to standard 3-minute kick samples only
  filter(SAMPLE_METHOD == "S3PO") %>%
  # ...to only include main abundance estimate methods
  filter(ANALYSIS_METHOD %in% c("ANAA", "ANLA", "ANLE"))

# Remove problematic species
invData <- invData %>%
  # Remove as the target shouldn't include this single INNS
  filter(PREFERRED_TAXON_NAME != "Musculium transversum") %>%
  # Not Capnia any more
  filter(PREFERRED_TAXON_NAME != "Capnia bifron")

# Rename for consistency
invData$PARENT_TAXON_NAME[which(invData$TAXON_NAME=="Baetis niger")] <- 
  "Nigrobaetis" #To keep consistent with Schedule 2 taxonomy

# ADD SPECIES NAME AND GROUP COLUMNS

# Add columms to populate species group name and type - 'schedule 2' or 'INNS'
invData$TAXON <- NA 
invData$GROUP <- NA

### IDENTIFY SCHEDULE 2 --------------------------------------------------------

# ADD TAXON LEVEL

# Add taxon level column; populate with "species"
s2$level <- "species"

# For species aggregates, change level to "genus"
s2$level[which(str_detect(s2$Species, " spp."))] <- "genus"

# Remove ".spp" from  species name
s2$Species <- str_replace(s2$Species, " spp.", "")

# Change name column to "taxon"
colnames(s2)[1] <- "taxon"

# CREATE SYNONYM LIST TO DEAL WITH SYNONYMS/GROUPS/AGGREGATES

# Add additional columns
s2$taxon.alt1 <- s2$taxon.alt2 <- s2$taxon.alt3 <- s2$taxon.alt4 <- NA

# Populate additional columns with synonyms
s2$taxon.alt1[s2$taxon=="Alboglossiphonia heteroclita"] <- "Glossiphonia heteroclita"
s2$taxon.alt1[s2$taxon=="Ampullaceana balthica"] <- "Radix balthica"
s2$taxon.alt1[s2$taxon=="Baetis muticus"] <- "Alainites muticus"
s2$taxon.alt1[s2$taxon=="Baetis scambus group"] <- "Baetis scambus"
s2$taxon.alt2[s2$taxon=="Baetis scambus group"] <- "Baetis scambus/fuscatus"
s2$taxon.alt3[s2$taxon=="Baetis scambus group"] <- "Baetis fuscatus"
s2$taxon.alt1[s2$taxon=="Baetis vernus"] <- "Baetis tenax"
s2$taxon.alt1[s2$taxon=="Bithynia leachii"] <- "Bithynia leachi"
s2$taxon.alt1[s2$taxon=="Caenis luctuosa group"] <- "Caenis luctuosa"
s2$taxon.alt2[s2$taxon=="Caenis luctuosa group"] <- "Caenis luctuosa/macrura"
s2$taxon.alt3[s2$taxon=="Caenis luctuosa group"] <- "Caenis macrura"
s2$taxon.alt4[s2$taxon=="Caenis luctuosa group"] <- "Caenis moesta"
s2$taxon.alt1[s2$taxon=="Dugesia polychroa group"] <- "Dugesia polychroa"
s2$taxon.alt2[s2$taxon=="Dugesia polychroa group"] <- "Dugesia lugubris or polychroa"
s2$taxon.alt3[s2$taxon=="Dugesia polychroa group"] <- "Dugesia lugubris"
s2$taxon.alt1[s2$taxon=="Gammarus pulex group"] <- "Gammarus pulex"
s2$taxon.alt2[s2$taxon=="Gammarus pulex group"] <- "Gammarus pulex/fossarum agg."
s2$taxon.alt3[s2$taxon=="Gammarus pulex group"] <- "Gammarus fossarum"
s2$taxon.alt1[s2$taxon=="Haliplus ruficollis group"] <- "Haliplus ruficollis"
s2$taxon.alt1[s2$taxon=="Lepidostoma basale"] <- "Lasiocephala basalis"
s2$taxon.alt1[s2$taxon=="Nebrioporus depressus group"] <- "Nebrioporus depressus"
s2$taxon.alt2[s2$taxon=="Nebrioporus depressus group"] <- "Nebrioporus elegans"
s2$taxon.alt1[s2$taxon=="Nemurella pictetii"] <- "Nemurella picteti"
s2$taxon.alt1[s2$taxon=="Perlodes mortoni"] <- "Perlodes microcephala"
s2$taxon.alt1[s2$taxon=="Simulium angustitarse group"] <- "Simulium angustitarse"
s2$taxon.alt1[s2$taxon=="Simulium argyreatum group"] <- "Simulium argyreatum"
s2$taxon.alt2[s2$taxon=="Simulium argyreatum group"] <- "Simulium argyreatum/variegatum"
s2$taxon.alt3[s2$taxon=="Simulium argyreatum group"] <- "Simulium variegatum"
s2$taxon.alt1[s2$taxon=="Simulium aureum group"] <- "Simulium aureum"
s2$taxon.alt1[s2$taxon=="Simulium cryophilum-vernum group"] <- "Simulium cryophilum"
s2$taxon.alt2[s2$taxon=="Simulium cryophilum-vernum group"] <- "Simulium vernum"
s2$taxon.alt1[s2$taxon=="Simulium ornatum group"] <- "Simulium ornatum"
s2$taxon.alt2[s2$taxon=="Simulium ornatum group"] <- "Simulium ornatum/intermedium/trifasciatum"
s2$taxon.alt3[s2$taxon=="Simulium ornatum group"] <- "Simulium intermedium"
s2$taxon.alt4[s2$taxon=="Simulium ornatum group"] <- "Simulium trifasciatum"
s2$taxon.alt1[s2$taxon=="Simulium tuberosum complex"] <- "Simulium tuberosum"
s2$taxon.alt1[s2$taxon=="Siphonoperla torrentium"] <- "Chloroperla torrentium"

# IDENTIFY SPECIES

# For each species in the schedule 2 list...
for (i in 1:NROW(s2)) {
  
  # Create taxon synonyms list for i
  taxaNames <- s2[i, ] %>% # For row i
    select(-level) %>% # Drop 'level' column
    select_if( ~ !(is.na(.))) # Drop NA columns
  
  # If taxon level is species ...
  if (s2$level[i] == "species") {
    
    # Identify invDatarows where PREFERRED_TAXON_NAME matches taxaNames
    iRows <- which(invData$PREFERRED_TAXON_NAME %in% taxaNames)
    
    # Else if taxon level is genus...
  } else {
    
    # Identify invData rows where PREFERRED_TAXON_NAME or
    # PARENT_TAXON_NAME matches taxaNames
    iRows <- which(invData$PREFERRED_TAXON_NAME %in% taxaNames |
                     invData$PARENT_TAXON_NAME %in% taxaNames)
  }
  
  # Assign invData taxon name and group based on iRows (match with s2)
  invData[iRows, "TAXON"] <- s2$taxon[i]
  invData[iRows, "GROUP"] <- "Schedule 2"
}

### IDENTIFY INNS --------------------------------------------------------------

# For each species in the inns list...
for (i in 1:NROW(inns)) {
  
  # Identify invData rows where PREFERRED_TAXON_NAME matches Scientific.name for i
  iRows <-
    which(invData$PREFERRED_TAXON_NAME %in% inns$Scientific.name[i])
  
  # Assign invData taxon name and group based on iRows (match with inns)
  invData[iRows, "TAXON"] <- inns$Scientific.name[i]
  invData[iRows, "GROUP"] <- "INNS"
}

# FILTER TO SCHEDULE 2 AND INNS SPECIES ----------------------------------------

# Remove species which do not have a group, e.g. are not schedule 2 or INNS
invData <- invData %>%
  filter(!is.na(GROUP))

### SAVE BIOSYS DATA -----------------------------------------------------------

# Save file
saveRDS(invData, file = paste0(dataDir, "Processed/Biosys/invData.Rds"))
