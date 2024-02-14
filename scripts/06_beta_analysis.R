# 07b - Statistically analysing the output beta statistics
# Joseph Everest
# February 2023, modified April 2023, May 2023, September 2023


# LOAD PACKAGES & FUNCTIONS ----

# Load packages
library(tidyverse)
library(viridis)
library(gridExtra)

# Load functions & themes
source("scripts/06_beta_analysis_FUNCTION.R")


# **[CHANGE]** - DECIDE ON PRE-PROCESSING DECISIONS ----

# ** [1] - Decisions
brightness <- "Yes" # Default = "No"
smoothing <- "No" # Default = "No"
PCA <- "No" # Default = "No"
remove.37 <- "No" # Default = "No"
top.hits.only <- "No" # Default = "No"
buffer <- "1" # Default = "1"
spectral.metric = "Euclidean" # Default = "Euclidean"

# Generate output folder paths
if (brightness == "No"){ filepath.brightness <- "" } else { filepath.brightness <- "_brightness_normalized" }
if (smoothing == "No"){ filepath.smoothing <- "" } else {filepath.smoothing <- "_smoothed"}
if (PCA == "No"){ filepath.PCA <- "" } else {filepath.PCA <- "_PCA"}
if (remove.37 == "No"){ filepath.37 <- "" } else { filepath.37 <- "_removed_37" }
if (top.hits.only == "No"){ filepath.top.hits <- "" } else { filepath.top.hits <- "_top_hits_only" }


# **[CHANGE]** - DECIDE WHICH MANTEL ANALYSES TO RUN ----

# Decide on which analyses to run
run.mantel.test.spatial <- TRUE
run.mantel.test.temporal <- TRUE


# MANTEL TESTS: SPATIAL COMPARISON ----

# Create a vector of the years for which data is available
saddle.years <- c("2017", "2018", "2019", "2020")

# Run function to calculate mantel statistics and plots
if (run.mantel.test.spatial == TRUE){run.mantel.spatial(saddle.years)}


# MANTEL TESTS: TEMPORAL COMPARISON ----

# Load in the temporal beta outputs
beta.temporal <- read.csv(paste0("outputs/output_beta_FULL_temporal_b", buffer, filepath.brightness,
                                 filepath.smoothing, filepath.37, filepath.top.hits, ".csv"))

# Create vector of plots
saddle.plots <- sort(unique(beta.temporal$PLOT))

# Run function to calculate mantel statistics and plots
if (run.mantel.test.temporal == TRUE){run.mantel.temporal(saddle.plots)}
