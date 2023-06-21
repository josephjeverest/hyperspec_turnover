# 07b - Statistically analysing the output beta statistics
# Joseph Everest
# February 2023, adpated April 2023


# LOAD PACKAGES, THEMES & FUNCTIONS ----

# Load packages
library(tidyverse)
library(viridis)
library(gridExtra)
library(vegan)

# Load themes and functions
source("scripts/EX1_ggplot_themes.R")
source("scripts/07_beta_analysis_FUNCTION.R")


# **[CHANGE]** - DECIDE WHETHER TO RETAIN PLOT 37 OR NOT ----

# Decision
retain.37 <- "Yes" # Default = "Yes"

# Generate output folder path
if (retain.37 == "Yes"){ filepath.37 <- "" } else { filepath.37 <- "_removed_37" }


# **[CHANGE]** - DECIDE WHETHER TO RETAIN PLOT 37 OR NOT ----

# Decision
top.hits.only <- "No" # Default = "No"

# Generate output folder path
if (top.hits.only == "No"){ filepath.top.hits <- "" } else { filepath.top.hits <- "_top_hits_only" }


# MANTEL TESTS: SPATIAL COMPARISON ----

# Create a vector of the years for which data is available
saddle.years <- c("2017", "2018", "2019", "2020")

# Run function to calculate mantel statistics and plots
# run.mantel.spatial(saddle.years, spectral.metric = "Euclidean")


# MANTEL TESTS: TEMPORAL COMPARISON ----

# Load in the temmporal beta outputs
beta.temporal <- read.csv(paste0("outputs/output_beta_FULL_temporal", filepath.37, filepath.top.hits, ".csv"))

# Create vector of plots
saddle.plots <- sort(unique(beta.temporal$Plot))

# Run function to calculate mantel statistics and plots
# run.mantel.temporal(saddle.plots, spectral.metric = "Euclidean")
