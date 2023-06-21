# 05b - Calculating Beta Diversity - Taxonomic, Functional & Spectral
# Joseph Everest
# February 2023, adapted April 2023, May 2023


# LOAD PACKAGES, THEMES & FUNCTIONS ----

# Load packages
library(tidyverse)
library(viridis)
library(gridExtra)
library(reshape2)
library(vegan)
library(biotools)

# Load themes and functions
source("scripts/EX1_ggplot_themes.R")
source("scripts/05_beta_diversity_FUNCTION.R")


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


# **[CHANGE]** - GENERATE INPUT PARAMETERS ---

# ** [1] - combined composition and trait data
composition.traits <- read.csv(paste0("outputs/output_saddle_composition_traits", filepath.top.hits, ".csv")) %>% 
  mutate(REMOVE.37 = ifelse(retain.37 == "No" & PLOT == 37, TRUE, FALSE)) %>% 
  filter(REMOVE.37 != TRUE) %>% 
  dplyr::select(-REMOVE.37)

# ** [2] - Spectral output data
spectra <- read.csv(paste0("outputs/output_saddle_spectra_b1", filepath.top.hits, ".csv")) %>% 
  mutate(REMOVE.37 = ifelse(retain.37 == "No" & Plot == 37, TRUE, FALSE)) %>% 
  filter(REMOVE.37 != TRUE) %>% 
  dplyr::select(-REMOVE.37)

# ** [3] - Spectral output data
spectra.PCA <- read.csv(paste0("outputs/output_saddle_spectra_b1_PCA", filepath.top.hits, ".csv")) %>% 
  mutate(REMOVE.37 = ifelse(retain.37 == "No" & Plot == 37, TRUE, FALSE)) %>% 
  filter(REMOVE.37 != TRUE) %>% 
  dplyr::select(-REMOVE.37)

# ** [4] - biomass output data
biomass <- read.csv(paste0("outputs/output_biomass", filepath.top.hits, ".csv")) %>% 
  mutate(REMOVE.37 = ifelse(retain.37 == "No" & Plot == 37, TRUE, FALSE)) %>% 
  filter(REMOVE.37 != TRUE) %>% 
  dplyr::select(-REMOVE.37)

# ** [5] - Create vector of the years in the dataset
composition.years <- c(2017, 2018, 2019, 2020)

# ** [6] - Create vector of the paired years in the dataset
composition.year.pairs <- c("2017_2018", "2017_2019", "2017_2020", "2018_2019", "2018_2020", "2019_2020")

# ** [7] - Create vector of plots in the dataset
composition.plots <- sort(unique(composition.traits$PLOT))


# CALCULATE BRAY-CURTIS DISSIMILARITY MATRIX - SPATIAL (beta.taxonomic.spatial) ----

# Run Bray-Curtis calculations for each year of dataset (~ 2-3 seconds)
# calc.beta.taxonomic.spatial(composition.years)

# Import output as dataframe
beta.taxonomic.spatial <- read.csv(paste0("outputs/output_beta_taxonomic_spatial", filepath.37, filepath.top.hits, ".csv"))


# CALCULATE FUNCTIONAL DISSIMILARITY MATRIX - SPATIAL (beta.functional.spatial) ----

# Run dissimilarity calculations for each year of dataset (~ 2-3 minutes)
# calc.beta.functional.spatial(composition.years)

# Import output as dataframe
beta.functional.spatial <- read.csv(paste0("outputs/output_beta_functional_spatial", filepath.37, filepath.top.hits, ".csv"))


# CALCULATE BIOMASS DISSIMILARITY MATRIX - SPATIAL (beta.biomass.spatial) ----

# Run biomass distance calculations for each year of dataset (~ 2-3 seconds)
# calc.beta.biomass.spatial(composition.years)

# Import output as dataframe
beta.biomass.spatial <- read.csv(paste0("outputs/output_beta_biomass_spatial", filepath.37, filepath.top.hits, ".csv"))


# CALCULATE SPECTRAL DISSIMILARITY MATRIX - SPATIAL (beta.spectral.spatial) ----

# Run spectral dissimilarity calculations for each year of dataset (~ 2-3 seconds)
# calc.beta.spectral.spatial(composition.years)

# Import output as dataframe
beta.spectral.spatial <- read.csv(paste0("outputs/output_beta_spectral_spatial", filepath.37, filepath.top.hits, ".csv"))


# CALCULATE SPECTRAL DISSIMILARITY MATRIX BY PCA - SPATIAL (beta.spectral.pca.spatial) ----

# Run spectral dissimilarity calculations for each year of dataset (~ 2-3 seconds)
# calc.beta.spectral.pca.spatial(composition.years)

# Import output as dataframe
beta.spectral.pca.spatial <- read.csv(paste0("outputs/output_beta_spectral_PCA_spatial", filepath.37, filepath.top.hits, ".csv"))


#  COMBINE DISSIMILARITY OUTPUTS - SPATIAL (beta.full.spatial) ----

# Bind the datasets together and convert to wide format
beta.full.spatial <- rbind(beta.taxonomic.spatial, beta.functional.spatial,
                           beta.biomass.spatial, beta.spectral.spatial,
                           beta.spectral.pca.spatial) %>% 
  mutate(Type = ifelse(str_detect(Type, pattern = "Spectral"), paste0(Type, "_", Method), Type),
         Type = paste0(Type, "_Dis")) %>% 
  dplyr::select(-Method) %>% 
  pivot_wider(names_from = "Type", values_from = "Dissimilarity")

# Export full wide dataframe
write.csv(beta.full.spatial, file = paste0("outputs/output_beta_FULL_spatial", filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)


# CALCULATE BRAY-CURTIS DISSIMILARITY MATRIX - TEMPORAL (beta.taxonomic.temporal) ----

# Run Bray-Curtis calculations for each year of dataset (~ 2-3 seconds)
# calc.beta.taxonomic.temporal(composition.year.pairs, composition.plots)

# Import output as dataframe
beta.taxonomic.temporal <- read.csv(paste0("outputs/output_beta_taxonomic_temporal", filepath.37, filepath.top.hits, ".csv"))


# CALCULATE FUNCTIONAL DISSIMILARITY MATRIX - TEMPORAL (beta.functional.temporal) ----

# Run distance calculations for each year of dataset (~ 2-3 minutes)
# calc.beta.functional.temporal(composition.year.pairs, composition.plots)

# Import output as dataframe
beta.functional.temporal <- read.csv(paste0("outputs/output_beta_functional_temporal", filepath.37, filepath.top.hits, ".csv"))


# CALCULATE BIOMASS DISSIMILARITY MATRIX - TEMPORAL (beta.biomass.temporal) ----

# Run distance calculations for each year of dataset (~ 2-3 seconds)
# calc.beta.biomass.temporal(composition.year.pairs, composition.plots)

# Import output as dataframe
beta.biomass.temporal <- read.csv(paste0("outputs/output_beta_biomass_temporal", filepath.37, filepath.top.hits, ".csv"))


# CALCULATE SPECTRAL DISSIMILARITY MATRIX - TEMPORAL (beta.spectral.temporal) ----

# Run distance calculations for each year of dataset (~ 2-3 seconds)
# calc.beta.spectral.temporal(composition.year.pairs, composition.plots)

# Import output as dataframe
beta.spectral.temporal <- read.csv(paste0("outputs/output_beta_spectral_temporal", filepath.37, filepath.top.hits, ".csv"))


# CALCULATE SPECTRAL DISSIMILARITY MATRIX BY PCA - TEMPORAL (beta.spectral.pca.temporal) ----

# Run spectral dissimilarity calculations for each year of dataset (~ 2-3 seconds)
# calc.beta.spectral.pca.temporal(composition.year.pairs, composition.plots)

# Import output as dataframe
beta.spectral.pca.temporal <- read.csv(paste0("outputs/output_beta_spectral_PCA_temporal", filepath.37, filepath.top.hits, ".csv"))


#  COMBINE DISSIMILARITY OUTPUTS - TEMPORAL (beta.full.temporal) ----

# Bind the datasets together and convert to wide format
beta.full.temporal <- rbind(beta.taxonomic.temporal, beta.functional.temporal,
                           beta.biomass.temporal, beta.spectral.temporal,
                           beta.spectral.pca.temporal) %>% 
  mutate(Type = ifelse(str_detect(Type, pattern = "Spectral"), paste0(Type, "_", Method), Type),
         Type = paste0(Type, "_Dis")) %>% 
  dplyr::select(-Method) %>% 
  pivot_wider(names_from = "Type", values_from = "Dissimilarity")

# Export full wide dataframe
write.csv(beta.full.temporal, file = paste0("outputs/output_beta_FULL_temporal", filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)


# PLOT SPECTRAL VS SPECTRAL PCA DISTANCES AGAINST ONE ANOTHER ----

# Comparison plot - spatial
(spectral.comparison.spatial <- ggplot(data = beta.full.spatial, aes(x = Spectral_Euclidean_Dis,
                                                                     y = Spectral_PCA_Euclidean_Dis,
                                                                     fill = as.factor(Year))) +
   geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
   stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
   
   ggpmisc::stat_fit_glance(method = 'lm', # Don't load this package (breaks tidyverse)
                            method.args = list(formula = y ~ x),  geom = 'text',
                            aes(label = paste0("\nR2 = ", signif(..r.squared.., digits = 3)))) +
   
   scale_fill_viridis(option = "magma", begin = 0, end = 1, direction = 1, discrete = TRUE) +
   scale_x_continuous(limits = c(min(beta.full.spatial$Spectral_Euclidean_Dis),
                                 max(beta.full.spatial$Spectral_Euclidean_Dis)), expand = c(0,0)) +
   scale_y_continuous(limits = c(min(beta.full.spatial$Spectral_PCA_Euclidean_Dis),
                                 max(beta.full.spatial$Spectral_PCA_Euclidean_Dis)), expand = c(0,0)) +
   labs(# title = "Spectral (Euclidean) Distance: Raw vs PCA",
        # subtitle = "SPATIAL: All Years",
        x = "\nSpectra (Raw)",
        y = "Spectra (PCA)\n",
        fill = "Year") +
   theme_1() +
   theme(legend.position = "right"))

# Comparison plot - temporal
(spectral.comparison.temporal <- ggplot(data = beta.full.temporal, aes(x = Spectral_Euclidean_Dis,
                                                                       y = Spectral_PCA_Euclidean_Dis,
                                                                       fill = Years)) +
    geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
    stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
    scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
    scale_x_continuous(limits = c(min(beta.full.temporal$Spectral_Euclidean_Dis),
                                  max(beta.full.temporal$Spectral_Euclidean_Dis)), expand = c(0,0)) +
    scale_y_continuous(limits = c(min(beta.full.temporal$Spectral_PCA_Euclidean_Dis),
                                  max(beta.full.temporal$Spectral_PCA_Euclidean_Dis)), expand = c(0,0)) +
    labs(title = "Spectral (Euclidean) Distance: Raw vs PCA",
         subtitle = "TEMPORAL: All Years",
         x = "Spectra (Raw)",
         y = "Spectra (PCA)",
         fill = "Year Pair") +
    theme_1() +
    theme(legend.position = "right"))

# Create panel of two comparison plots
spectral.comparison.panel <- grid.arrange(spectral.comparison.spatial, spectral.comparison.temporal, ncol = 2)

# Export plot
ggsave(spectral.comparison.panel, filename = paste0("outputs/figures/beta_spectral_vs_PCA_comparison", filepath.37, filepath.top.hits, ".png"),
       width = 16, height = 7.5)


# FURTHER EXPLORATION ----

# Import dataframe of plot vegetation classes
veg.class <- read.csv("data/nwt_biomass.csv") %>% 
  dplyr::select(grid_pt, veg_class) %>% 
  distinct()

# Produce wide format data to determine the plots with abnormally high functional distance
beta.explore.w <- beta.functional.spatial %>% 
  filter(Dissimilarity > 0.22) %>% 
  left_join(., veg.class, by = c("PLOT_1" = "grid_pt")) %>% 
  rename(VEG_1 = veg_class) %>% 
  left_join(., veg.class, by = c("PLOT_2" = "grid_pt")) %>% 
  rename(VEG_2 = veg_class) %>% 
  dplyr::select(Year, PLOT_1, PLOT_2, Dissimilarity, VEG_1, VEG_2)

# Determine whether every instance includes plot 37
beta.explore.37 <- beta.explore.w %>% 
  mutate(includes_37 = ifelse(PLOT_1 == 37 | PLOT_2 == 37, TRUE, FALSE))
    # Every record
      # PLOT 37 is an ST (shrub tundra) plot


# EXPORT EXAMPLE MATRIX - SPATIAL ----

# Load in the matrices
matrices.functional.spatial <- get(load(paste0("outputs/output_beta_functional_spatial.RData")))

# Trimmed matrices
matrix.spatial <- as.data.frame(as.matrix(matrices.functional.spatial[[1]]))[1:10, 1:10]

# Write matrix to .csv
write.csv(matrix.spatial, file = paste0("outputs/output_example_matrix_spatial.csv"))


# EXPORT EXAMPLE MATRIX - TEMPORAL ----

# Load in temporal beta outputs
beta.temporal <- read.csv(paste0("outputs/output_beta_FULL_temporal.csv"))

# Cut to correct years
beta.temporal.years <- beta.temporal %>% 
  separate(Years, into = c("Year_1", "Year_2"), sep = "_")

# Cut to single plot (2)
beta.temporal.cut <- beta.temporal.years %>% 
  filter(Plot == 2)

# Cut dataframe to required information for single matrix
beta.temporal.df.f <- dplyr::select(beta.temporal.cut, Year_1, Year_2, Functional_Dis)

# Create a vector of the years for which data is available
saddle.years <- c("2017", "2018", "2019", "2020")

# Produce matrix
beta.temporal.matrix.f <- with(beta.temporal.df.f,
                               structure(Functional_Dis, Size = length(saddle.years), Labels = saddle.years,
                                         Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))

# Trimmed matrices
matrix.temporal <- as.data.frame(as.matrix(beta.temporal.matrix.f))

# Write matrix to .csv
write.csv(matrix.temporal, file = paste0("outputs/output_example_matrix_temporal.csv"))
