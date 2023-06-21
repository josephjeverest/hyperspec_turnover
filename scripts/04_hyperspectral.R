# 04b - Processing NEON hyperspectral
# Joseph Everest, with Sarah Elmendorf
# February 2023, adapted April 2023


# LOAD PACKAGES & FUNCTIONS ----

# Load packages
library(tidyverse)
library(sp)
library(rgdal)
library(viridis)
library(gridExtra)
library(raster)
# install.packages("BiocManager")
# BiocManager::install("rhdf5")
library(rhdf5)
library(neonhs)


# Load functions for working with NEON data and plotting
source("scripts/EX1_ggplot_themes.R")
source("scripts/04_hyperspectral_FUNCTION.R")


# **[CHANGE]** - DECIDE WHETHER TO USE ALL HITS OR TOP HITS ONLY ----

# Decision
top.hits.only <- "No" # Default = "No"

# Generate output folder path
if (top.hits.only == "No"){ filepath.top.hits <- "" } else { filepath.top.hits <- "_top_hits_only" }


# USE MOST-RECENT COORDINATES PROVIDED TO FIND PLOT LOCATIONS ----

# Load in saddle data to get plot locations
saddle.1 <- read.csv(paste0("outputs/output_saddle_composition_traits", filepath.top.hits, ".csv"))

# Determine list of plots retained in composition data (78)
saddle.plots <- sort(unique(saddle.1$PLOT))

# Load in plot locations
saddle.locations.raw <- read.csv("data/nwt_lter_locations_utm_II_v10.csv")

# Clean up dataframe to retain only required information
saddle.locations <- saddle.locations.raw %>% 
  dplyr::select(ALT_SITECODE, UTM_E, UTM_N) %>% 
  filter(str_detect(ALT_SITECODE, pattern = "Pt Quad")) %>% 
  rename(PLOT = ALT_SITECODE,
         EASTING = UTM_E,
         NORTHING = UTM_N) %>% 
  mutate(PLOT = str_remove(PLOT, pattern = "Saddle Pt Quad "),
         PLOT = str_remove(PLOT, pattern = "^0+"),
         PLOT = str_replace(PLOT, pattern = "A", replacement = "1"),
         PLOT = as.integer(PLOT)) %>% 
  filter(PLOT %in% saddle.plots) %>% 
  arrange(PLOT) %>% 
  filter(PLOT != 3) # Remove plots removed in composition dataframe also

# Export list of saddle locations
write.csv(saddle.locations, file = paste0("outputs/output_saddle_locations", filepath.top.hits, ".csv"), row.names = FALSE)


# EXTRACT TILE METADATA ----

# Create a vector of all the HDF5 files to work with
tiles.folderpath <- "data/NEON/tiles/"
tiles.filename <- list.files(path = tiles.folderpath, pattern = "*.h5")
tiles.filepath <- paste0(tiles.folderpath, tiles.filename)

# Extract metadata using function
# extract.tile.metadata(tiles.filepath)

# Import output of function
tile.information <- read.csv(paste0("outputs/output_NIWOT_tile_metadata", filepath.top.hits, ".csv"))

# Remove intermediate objects
rm(tiles.folderpath, tiles.filename)


# WORKING WITH SADDLE PLOTS ----

# Create a vector of the years in which hyperspectral data was available
saddle.years <- c("2017", "2018", "2019", "2020")

# Set buffer to run with
buffer <- 1

# Run pixel extraction function
# extract.pixel.spectra(tiles.filepath, saddle.years, buffer)

# Import output of function
spectra.output <- read.csv(paste0("outputs/output_saddle_spectra_b", buffer, filepath.top.hits, ".csv"))

# Remove intermediate objects
rm(saddle.years, tiles.filepath)


# BAND PER YEAR TESTS ----

# Create checks dataframe to see all years/plots have same number of bands/spectra
spectra.checks <- spectra.output %>%
  group_by(Year) %>%
  mutate(num_plots = length(unique(Plot))) %>%
  ungroup() %>%
  group_by(Year, Plot) %>%
  mutate(num_bands = length(unique(Band))) %>%
  ungroup() %>%
  dplyr::select(Year, Plot, num_plots, num_bands) %>%
  distinct()


# GENERATE OUTPUT PLOTS BETWEEN BUFFERS ----

# Create vector of three buffers to compare (in metres, must have already run above)
buffers <- c(0, 1, 3)

# Run the function to generate output plots
# extract.comparison.plots(buffers)


# RUN PCA DIMENSIONALITY REDUCTION ----

# Restructure the spectral outputs for running PCAs
spectra.restructure <- spectra.output %>%
  mutate(PlotYear = paste0(Plot, ":", Year)) %>%
  dplyr::select(PlotYear, Band, smooth_Reflectance) %>%
  pivot_wider(names_from = "Band", values_from = "smooth_Reflectance") %>%
  column_to_rownames(var = "PlotYear")

# Run PCA on the restructured spectral data
spectra.PCA <- prcomp(spectra.restructure, scale = TRUE)

# Reverse sign of results as eigenvectors negative in R by default
spectra.PCA$x <- -1*spectra.PCA$x

# Get results of the spectra PCA
spectra.PCA.results <- data.frame(spectra.PCA$x)

# Calculate variance explained by each principle component
spectra.variance <- spectra.PCA$sdev^2 / sum(spectra.PCA$sdev^2)

# Retain PCAs that incorporate 99% of the variance in the data (5)
spectra.variance.99 <- sum(spectra.variance[1:5])

# Retain only the required PCs
spectra.PCA.selected <- spectra.PCA.results[, 1:5]

# Restructure data back to the original format for mantel test calculations
spectra.PCA.restructure <- spectra.PCA.selected %>%
  rownames_to_column(.) %>%
  rename(PlotYear = rowname) %>%
  pivot_longer(names_to = "PC", values_to = "PC_value", cols = 2:6) %>%
  separate(PlotYear, into = c("Plot", "Year"), sep = ":") %>%
  mutate(Plot = as.integer(Plot),
         Year = as.integer(Year))

# Create NIR and NDVI key for rejoining back to PCA data
spectra.NIR.NDVI <- spectra.output %>%
  dplyr::select(Plot, Year, mean_NIR, NDVI_broad, NDVI_narrow) %>%
  distinct()

# Join NIR and NDVI values to PCA data
spectra.PCA.output <- spectra.PCA.restructure %>%
  left_join(., spectra.NIR.NDVI, by = c("Plot" = "Plot", "Year" = "Year"))

# Export PCA output
write.csv(spectra.PCA.output, file = paste0("outputs/output_saddle_spectra_b", buffer, "_PCA", filepath.top.hits, ".csv"),
          row.names = FALSE)
