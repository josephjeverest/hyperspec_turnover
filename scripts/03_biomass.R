# 03 - Saddle Biomass Data
# Joseph Everest
# February 2023


# LOAD PACKAGES ----

# Load packages
library(tidyverse)


# **[CHANGE]** - DECIDE WHETHER TO USE ALL HITS OR TOP HITS ONLY ----

# Decision
top.hits.only <- "No" # Default = "No"

# Generate output folder path
if (top.hits.only == "No"){ filepath.top.hits <- "" } else { filepath.top.hits <- "_top_hits_only" }


# IMPORT DATA (biomass.1) ----

# Load Niwot biomass data
biomass.1 <- read.csv("data/nwt_biomass.csv")


# DATA MANIPULATION (biomass.2) ----

# Trim dataframe to correct years columns
biomass.2 <- biomass.1 %>% 
  dplyr::select(year, subsample, grid_pt, NPP) %>%
  group_by(year, grid_pt) %>% 
  summarise(NPP = mean(NPP)) %>% # Average the two subsamples per year
  ungroup() %>% 
  rename(Year = year, Plot = grid_pt) %>% 
  filter(Year %in% c(2017, 2018, 2019, 2020))


# PLOT MATCHING (biomass.3) ----

# Determine plots retained in composition/trait/spectral dataframes
unique.plots <- unique(read.csv(paste0("outputs/output_saddle_cover", filepath.top.hits, ".csv"))$PLOT)

# Remove plots not in those dataframes
biomass.3 <- biomass.2 %>% 
  filter(Plot %in% unique.plots)
  

# OUTPUT DATAFRAME ----

# Write dataframe to .csv
write.csv(biomass.3, paste0("outputs/output_biomass", filepath.top.hits, ".csv"), row.names = FALSE)
