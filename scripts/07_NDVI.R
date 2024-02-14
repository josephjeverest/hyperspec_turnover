# 08 - Extract NDVI and compare to biomass
# Joseph Everest
# March 2023, modified April 2023, May 2023, September 2023


# LOAD PACKAGES, FUNCTIONS & THEMES ----

# Packages
library(tidyverse)
library(reshape2)
library(gridExtra)
library(vegan)
library(sp)
library(rgdal)
library(raster)
# install.packages("BiocManager")
# BiocManager::install("rhdf5")
library(rhdf5)
library(neonhs)
library(viridis)

# Functions and themes
source("scripts/07_NDVI_FUNCTION.R")
source("scripts/EX1_ggplot_themes.R")


# **[CHANGE]** - DECIDE ON PRE-PROCESSING DECISIONS ----

# ** [1] - Decisions
remove.37 <- "No" # Default = "No"
top.hits.only <- "No" # Default = "No"
buffer <- "1" # Default = "1"
exponentiate <- "No" # Default = "No"
spectral.metric = "Euclidean" # Default = "Euclidean"

# Generate output folder paths
if (remove.37 == "No"){ filepath.37 <- "" } else { filepath.37 <- "_removed_37" }
if (top.hits.only == "No"){ filepath.top.hits <- "" } else { filepath.top.hits <- "_top_hits_only" }
if (exponentiate == "No"){ filepath.exponentiate <- ""} else { filepath.exponentiate <- "_exponentiated" }


# **[CHANGE]** - DECIDE WHICH ANALYSES TO RUN ----

# Decide on which analyses to run
run.ndvi.extraction <- FALSE
run.ndvi.beta.spatial <- FALSE
run.ndvi.beta.temporal <- FALSE
run.ndvi.mantel.spatial <- FALSE
run.ndvi.mantel.temporal <- FALSE


# LOAD NDVI DATA ----

# Load in processed spectra (not PCA)
spectra <- read.csv(paste0("outputs/output_saddle_spectra_b", buffer, filepath.top.hits, ".csv"))


# EXTRACT NDVI FOR SADDLE LOCATIONS ----

# Extract NDVI from the full hyperspectra
if (run.ndvi.extraction == TRUE){extract.ndvi(spectra)}

# Load in NDVI output
saddle.ndvi <- read.csv(paste0("outputs/output_saddle_ndvi_b", buffer, filepath.top.hits, ".csv"))


# CALCULATE BETA STATISTICS ON EXTRACTED NDVI VALUES - SPATIAL ----

# Create vector of years in data
composition.years <- c(2017, 2018, 2019, 2020)

# Run NDVI spatial beta calculations
if (run.ndvi.beta.spatial == TRUE){calc.beta.ndvi.spatial(composition.years)}

# Import output as dataframe
beta.ndvi.spatial <- read.csv(paste0("outputs/output_beta_ndvi_spatial_b", buffer, "_", spectral.metric,
                                     filepath.exponentiate, filepath.top.hits, ".csv"))


# CALCULATE BETA STATISTICS ON EXTRACTED NDVI VALUES - TEMPORAL ----

# Load combined composition and trait data
composition.traits <- read.csv(paste0("outputs/output_saddle_composition_traits.csv"))

# Create vector of plots in the dataset
composition.plots <- sort(unique(composition.traits$PLOT))

# Remove composition traits object
rm(composition.traits)

# Create vector of the paired years in the dataset
composition.year.pairs <- c("2017_2018", "2017_2019", "2017_2020", "2018_2019", "2018_2020", "2019_2020")

# Run NDVI temporal beta calculations
if (run.ndvi.beta.temporal == TRUE){calc.beta.ndvi.temporal(composition.year.pairs, composition.plots)}

# Import output as dataframe
beta.ndvi.temporal <- read.csv(paste0("outputs/output_beta_ndvi_temporal_b", buffer, "_", spectral.metric,
                                      filepath.exponentiate, filepath.top.hits, ".csv"))


# LOAD EXISTING BETA METRICS - SPATIAL & TEMPORAL ----

# Load in the complete datasets of dissimilarities
beta.spatial <- read.csv(paste0("outputs/output_beta_FULL_spatial_b", buffer, filepath.37, filepath.top.hits, ".csv"))
beta.temporal <- read.csv(paste0("outputs/output_beta_FULL_temporal_b", buffer, filepath.37, filepath.top.hits, ".csv"))

# Join the ndvi information to these datasets
beta.spatial.full <- beta.ndvi.spatial %>% 
  rename(NDVI_Dis = Dissimilarity) %>% 
  dplyr::select(Year, PLOT_1, PLOT_2, NDVI_Dis) %>% 
  left_join(beta.spatial, ., by = c("Year" = "Year", "PLOT_1" = "PLOT_1", "PLOT_2" = "PLOT_2"))

beta.temporal.full <- beta.ndvi.temporal %>% 
  rename(NDVI_Dis = Dissimilarity) %>% 
  dplyr::select(PLOT, Years, NDVI_Dis) %>% 
  left_join(beta.temporal, ., by = c("PLOT" = "PLOT", "Years" = "Years"))


# RUN STATISTICAL TESTS BETWEEN NDVI AND OTHER BETA METRICS ----

# Run mantel tests for NDVI spatial comparisons
run.mantel.spatial.ndvi(composition.years)

# Run mantel tests for NDVI temporal comparisons
run.mantel.temporal.ndvi(composition.plots, composition.years)
