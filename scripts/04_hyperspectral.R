# 04b - Processing NEON hyperspectral
# Joseph Everest, with Sarah Elmendorf
# February 2023, modified April 2023, May 2023, September 2023


# LOAD PACKAGES & FUNCTIONS ----

# Load packages
library(tidyverse)
library(neonUtilities)
library(sp)
library(rgdal)


# Load functions for working with NEON data and plotting
source("scripts/EX1_ggplot_themes.R")
source("scripts/04_hyperspectral_FUNCTION.R")


# **[CHANGE]** - DECIDE WHICH ANALYSES TO RUN ----

download.NEON.tiles <- FALSE
run.tile.metadata <- FALSE
run.spectra.extraction <- FALSE
run.buffer.comparison <- FALSE


# **[CHANGE]** - DECIDE ON PRE-PROCESSING DECISIONS & WHETHER TO USE TOP HITS ONLY OR NOT ----

# Decisions
brightness <- "Yes" # Default = "Yes"
smoothing <- "No" # Default = "No"
top.hits.only <- "No" # Default = "No"
buffer <- 1

# Generate output folder paths
if (brightness == "No"){ filepath.brightness <- "" } else { filepath.brightness <- "_brightness_normalized" }
if (smoothing == "No"){ filepath.smoothing <- "" } else {filepath.smoothing <- "_smoothed"}
if (PCA == "No"){ filepath.PCA <- "" } else {filepath.PCA <- "_PCA"}
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

# Generate vectors of unique eastings & northings for which to download spectral data
saddle.eastings <- c(unique(saddle.locations$EASTING))
saddle.northings <- c(unique(saddle.locations$NORTHING))

# Remove intermediate objects
rm(saddle.plots, saddle.locations.raw)


# DOWNLOAD HYPERSPECTRAL DATA (only run with top.hits.only == "No") ----

if (download.NEON.tiles == TRUE){
  
  
  # Create directory to move files to
  dir.create("data/NEON/tiles")
  
  # Create counter of Niwot revisit for determining file renaming (2017 was revisit 1, 2018 was 2...)
  niwot.revisit <- 0
  
  
  # Run for loop for each year
  for (year in c("2017", "2018", "2019", "2020")){
    
    
    # Incresae Niwot revisit counter by 1
    niwot.revisit <- niwot.revisit + 1
    
    # Run download function (has to be run on separate years, can't vector them)
    byTileAOP("DP3.30006.001",
              site = "NIWO",
              year = year, # Data available from 2017 to 2020
              easting = saddle.eastings,
              northing = saddle.northings,
              buffer = 1,
              check.size = TRUE,
              savepath = "data/NEON/BOB",
              token = NA_character_)
    
    
    # Run loop to rename (move) files to easier location with year in filename
    for (tile.extent in c("4433", "4434")){
      
      
      file.copy(paste0("data/NEON/DP3.30006.001/neon-aop-products/", year, "/FullSite/D13/",
                       year, "_NIWO_", niwot.revisit, "/L3/Spectrometer/Reflectance/NEON_D13_NIWO_DP3_",
                       "449000_", tile.extent, "000_reflectance.h5"),
                paste0("data/NEON/tiles/NEON_D13_NIWO_DP3_", year, "_449000_", tile.extent,
                       "000_reflectance.h5"))
      
      
    } # End of extent loop
    
    
  } # End of year loop
  
  
} # End of if statement


# Remove intermediate objects
rm(saddle.northings, saddle.eastings)


# EXTRACT TILE METADATA ----

# Create a vector of all the HDF5 files to work with
tiles.folderpath <- "data/NEON/tiles/"
tiles.filename <- list.files(path = tiles.folderpath, pattern = "*.h5")
tiles.filepath <- paste0(tiles.folderpath, tiles.filename)

# Extract metadata using function
if (run.tile.metadata == TRUE){extract.tile.metadata(tiles.filepath)}

# Import output of function
tile.information <- read.csv(paste0("outputs/output_NIWOT_tile_metadata", filepath.top.hits, ".csv"))

# Remove intermediate objects
rm(tiles.folderpath, tiles.filename)


# WORKING WITH SADDLE PLOTS ----

# Create a vector of the years in which hyperspectral data was available
saddle.years <- c("2017", "2018", "2019", "2020")

# Run pixel extraction function
if (run.spectra.extraction == TRUE){extract.pixel.spectra(tiles.filepath, saddle.years, buffer)}

# Import output of function
spectra.output <- read.csv(paste0("outputs/output_saddle_spectra_b", buffer,
                                  filepath.brightness, filepath.smoothing, filepath.top.hits, ".csv"))

# Remove intermediate objects
rm(saddle.years, tiles.filepath)


# BAND PER YEAR TESTS ----

# Create checks dataframe to see all years/plots have same number of bands/spectra
spectra.checks <- spectra.output %>%
  group_by(Year) %>%
  mutate(num_plots = length(unique(PLOT))) %>%
  ungroup() %>%
  group_by(Year, PLOT) %>%
  mutate(num_bands = length(unique(Band))) %>%
  ungroup() %>%
  dplyr::select(Year, PLOT, num_plots, num_bands) %>%
  distinct()


# GENERATE OUTPUT PLOTS BETWEEN BUFFERS ----

# Create vector of three buffers to compare (in metres, must have already run above)
buffers <- c(0, 1, 3)

# Run the function to generate output plots
if(run.buffer.comparison == TRUE){extract.comparison.plots(buffers)}


# RUN PCA DIMENSIONALITY REDUCTION ----

  
# Restructure the spectral outputs for running PCAs
spectra.restructure <- spectra.output %>%
  mutate(PlotYear = paste0(PLOT, ":", Year)) %>%
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
  separate(PlotYear, into = c("PLOT", "Year"), sep = ":") %>%
  mutate(PLOT = as.integer(PLOT),
         Year = as.integer(Year))

# Create NIR and NDVI key for rejoining back to PCA data
spectra.NIR.NDVI <- spectra.output %>%
  dplyr::select(PLOT, Year, mean_NIR, NDVI_broad, NDVI_narrow) %>%
  distinct()

# Join NIR and NDVI values to PCA data
spectra.PCA.output <- spectra.PCA.restructure %>%
  left_join(., spectra.NIR.NDVI, by = c("PLOT" = "PLOT", "Year" = "Year"))

# Export PCA output
write.csv(spectra.PCA.output, file = paste0("outputs/output_saddle_spectra_b", buffer,
                                            filepath.brightness, filepath.smoothing, "_PCA",
                                            filepath.top.hits, ".csv"), row.names = FALSE)
