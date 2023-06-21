# 08 - Extract NDVI and compare to biomass
# Joseph Everest
# March 2023


# LOAD PACKAGES, FUNCTIONS & THEMES ----

# Packages
library(tidyverse)
library(viridis)
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

# Functions and themes
source("scripts/08_NDVI_FUNCTION.R")
source("scripts/EX1_ggplot_themes.R")


# LOAD NDVI ----

# Load in processed spectra
spectra <- read.csv("outputs/output_saddle_spectra_b1.csv")

# Extract NDVI for saddle locations
extract.ndvi.buffer(spectra)

# Load in NDVI output
saddle.ndvi <- read.csv("outputs/output_saddle_ndvi.csv")


# CALCULATE NDVI BETA MATRICES - SPATIAL ----

# Create vector of years in data
composition.years <- c(2017, 2018, 2019, 2020)

# Run NDVI spatial beta calculations
calc.beta.ndvi.spatial(composition.years)

# Import output as dataframe
beta.ndvi.spatial <- read.csv("outputs/output_beta_ndvi_spatial.csv")


# CALCULATE NDVI BETA MATRICES - TEMPORAL ----

# Load combined composition and trait data
composition.traits <- read.csv("outputs/output_saddle_composition_traits.csv")

# Create vector of plots in the dataset
composition.plots <- sort(unique(composition.traits$PLOT))

# Remove composition traits object
rm(composition.traits)

# Create vector of the paired years in the dataset
composition.year.pairs <- c("2017_2018", "2017_2019", "2017_2020", "2018_2019", "2018_2020", "2019_2020")

# Run NDVI temporal beta calculations
calc.beta.ndvi.temporal(composition.year.pairs, composition.plots)

# Import output as dataframe
beta.ndvi.temporal <- read.csv("outputs/output_beta_ndvi_temporal.csv")


# LOAD EXISTING BETA METRICS (SPATIAL AND TEMPORAL) ----

# Load in the complete datasets of dissimilarities
beta.spatial <- read.csv("outputs/output_beta_full_spatial.csv")
beta.temporal <- read.csv("outputs/output_beta_full_temporal.csv")

# Join the ndvi information to these datasets
beta.spatial.full <- beta.ndvi.spatial %>% 
  rename(NDVI_Dis = Dissimilarity) %>% 
  dplyr::select(Year, PLOT_1, PLOT_2, NDVI_Dis) %>% 
  left_join(beta.spatial, ., by = c("Year" = "Year", "PLOT_1" = "PLOT_1", "PLOT_2" = "PLOT_2"))

beta.temporal.full <- beta.ndvi.temporal %>% 
  rename(NDVI_Dis = Dissimilarity) %>% 
  dplyr::select(Plot, Years, NDVI_Dis) %>% 
  left_join(beta.temporal, ., by = c("Plot" = "Plot", "Years" = "Years"))


# PRODUCE COMPARISON VISUALISATIONS OVER SPACE AND TIME ----

# Generate spatial plots
beta.visualisations.spatial.ndvi(beta.spatial.full)

# Generate temporal plots
beta.visualisations.temporal.ndvi(beta.temporal.full)


# RUN STATISTICAL TESTS BETWEEN NDVI AND OTHER BETA METRICS ----

# Run mantel tests for NDVI spatial comparisons
run.mantel.spatial.ndvi(composition.years)

# Run mantel tests for NDVI temporal comparisons
run.mantel.temporal.ndvi(composition.plots, composition.years)

