# 06 - Comparing Beta Diversity - Taxonomic, Functional & Spectral
# Joseph Everest
# February 2023, adapted March 2023


# LOAD PACKAGES, THEMES & FUNCTIONS ----

# Load packages
library(tidyverse)
library(viridis)
library(gridExtra)
library(broom)

# Load themes and functions
source("scripts/EX1_ggplot_themes.R")
source("scripts/06_beta_visualisation_FUNCTION.R")


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


# LOAD DATA ----

# Load in the complete datasets of dissimilarities
beta.spatial <- read.csv(paste0("outputs/output_beta_FULL_spatial", filepath.37, filepath.top.hits, ".csv"))
beta.temporal <- read.csv(paste0("outputs/output_beta_FULL_temporal", filepath.37, filepath.top.hits, ".csv"))


# EARLY VISUALISATIONS: SPACE ----

# Run spatial visualisation function [OPTIONS: "Euclidean", "Manhattan", "SAM"]
# beta.visualisations.spatial(beta.spatial, spectral.metric = "Euclidean")


# EARLY VISUALISATIONS: TEMPORAL ----

# Run spatial visualisation function [OPTIONS: "Euclidean", "Manhattan", "SAM"]
# beta.visualisations.temporal(beta.temporal, spectral.metric = "Euclidean")


# BOXPLOT: SPECTRAL DISTANCE BY YEAR ----

# Create a boxplot of spectral distance by year pairs
(boxplot.spectral <- ggplot(data = beta.temporal, aes(x = as.factor(Years), y = Spectral_PCA_Euclidean_Dis, fill = Years)) +
   geom_boxplot() +
   scale_fill_viridis(option = "viridis", begin = 0.3, end = 1, direction = 1, discrete = TRUE) +
   scale_y_continuous(limits = c(0, 35), expand = c(0,0)) +
   labs(title = "Spectral PCA Distance",
        subtitle = "By Year Pair",
        x = "Year Pair",
        y = "Spectral (Euclidean) Distance") +
   theme_1() +
   theme(legend.position = "none"))

# Export plot
ggsave(boxplot.spectral, filename = paste0("outputs/figures/boxplot_spectral_distance_year_pairs", filepath.37, filepath.top.hits, ".png"), width = 8, height = 5)

# Create dataframes for each year
beta.no.2017 <- filter(beta.temporal, !str_detect(Years, pattern = "2017")) %>% 
  mutate(YearNo = 2017) %>% 
  dplyr::select(YearNo, Spectral_PCA_Euclidean_Dis)
beta.no.2018 <- filter(beta.temporal, !str_detect(Years, pattern = "2018")) %>% 
  mutate(YearNo = 2018) %>% 
  dplyr::select(YearNo, Spectral_PCA_Euclidean_Dis)
beta.no.2019 <- filter(beta.temporal, !str_detect(Years, pattern = "2019")) %>% 
  mutate(YearNo = 2019) %>% 
  dplyr::select(YearNo, Spectral_PCA_Euclidean_Dis)
beta.no.2020 <- filter(beta.temporal, !str_detect(Years, pattern = "2020")) %>% 
  mutate(YearNo = 2020) %>% 
  dplyr::select(YearNo, Spectral_PCA_Euclidean_Dis)

# Join into one dataframe for plotting
beta.no.year <- rbind(beta.no.2017, beta.no.2018, beta.no.2019, beta.no.2020)

# Create a boxplot of spectral distance by missing years
(boxplot.spectral.missing <- ggplot(data = beta.no.year, aes(x = as.factor(YearNo), y = Spectral_PCA_Euclidean_Dis, fill = as.factor(YearNo))) +
    geom_boxplot() +
    scale_fill_viridis(option = "viridis", begin = 0.3, end = 1, direction = 1, discrete = TRUE) +
    scale_y_continuous(limits = c(0, 35), expand = c(0,0)) +
    labs(title = "Spectral PCA Distance",
         subtitle = "By Year NOT Included",
         x = "Year NOT Included",
         y = "Spectral (Euclidean) Distance") +
    theme_1() +
    theme(legend.position = "none"))

# Export plot
ggsave(boxplot.spectral.missing, filename = paste0("outputs/figures/boxplot_spectral_distance_missing_year.", filepath.37, filepath.top.hits, "png"), width = 8, height = 5)
