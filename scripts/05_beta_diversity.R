# 05b - Calculating Beta Diversity - Taxonomic, Functional & Spectral
# Joseph Everest
# February 2023, modified April 2023, May 2023, September 2023


# LOAD PACKAGES & FUNCTIONS ----

# Load packages
library(tidyverse)

# Load functions
source("scripts/05_beta_diversity_FUNCTION.R")


# **[CHANGE]** - DECIDE ON PRE-PROCESSING DECISIONS ----

# ** [1] - Decisions
brightness <- "Yes" # Default = "No"
smoothing <- "No" # Default = "No"
PCA <- "No" # Default = "No"
remove.37 <- "No" # Default = "No"
top.hits.only <- "No" # Default = "No"
buffer <- "1" # Default = "1"

# Generate output folder paths
if (brightness == "No"){ filepath.brightness <- "" } else { filepath.brightness <- "_brightness_normalized" }
if (smoothing == "No"){ filepath.smoothing <- "" } else {filepath.smoothing <- "_smoothed"}
if (remove.37 == "No"){ filepath.37 <- "" } else { filepath.37 <- "_removed_37" }
if (top.hits.only == "No"){ filepath.top.hits <- "" } else { filepath.top.hits <- "_top_hits_only" }


# **[CHANGE]** - DECIDE WHAT ANALYSES TO RUN ----

# ** [2] - Decide which spatial analyses to run
run.spatial.taxonomic <- TRUE
run.spatial.functional <- TRUE
run.spatial.biomass <- TRUE
run.spatial.spectral <- TRUE

# ** [3] - Decide which temporal analyses to run
run.temporal.taxonomic <- TRUE
run.temporal.functional <- TRUE
run.temporal.biomass <- TRUE
run.temporal.spectral <- TRUE


# GENERATE INPUT PARAMETERS ---

# Combined composition and trait data
composition.traits <- read.csv(paste0("outputs/output_saddle_composition_traits", filepath.top.hits, ".csv")) %>% 
  mutate(REMOVE.37 = ifelse(remove.37 == "Yes" & PLOT == 37, TRUE, FALSE)) %>% 
  filter(REMOVE.37 != TRUE) %>% 
  dplyr::select(-REMOVE.37)

# Biomass output data
biomass <- read.csv(paste0("outputs/output_biomass", filepath.top.hits, ".csv")) %>% 
  mutate(REMOVE.37 = ifelse(remove.37 == "Yes" & PLOT == 37, TRUE, FALSE)) %>% 
  filter(REMOVE.37 != TRUE) %>% 
  dplyr::select(-REMOVE.37)
  
# Spectral output data
spectra <- read.csv(paste0("outputs/output_saddle_spectra_b", buffer, filepath.brightness,
                           filepath.smoothing, filepath.top.hits, ".csv")) %>% 
  mutate(REMOVE.37 = ifelse(remove.37 == "Yes" & PLOT == 37, TRUE, FALSE)) %>% 
  filter(REMOVE.37 != TRUE) %>% 
  dplyr::select(-REMOVE.37)

# Spectral PCA output data
spectra.PCA <- read.csv(paste0("outputs/output_saddle_spectra_b", buffer, filepath.brightness,
                               filepath.smoothing, "_PCA", filepath.top.hits, ".csv")) %>% 
  mutate(REMOVE.37 = ifelse(remove.37 == "Yes" & PLOT == 37, TRUE, FALSE)) %>% 
  filter(REMOVE.37 != TRUE) %>% 
  dplyr::select(-REMOVE.37)

# Create vector of the years in the dataset
composition.years <- c(2017, 2018, 2019, 2020)

# Create vector of the paired years in the dataset
composition.year.pairs <- c("2017_2018", "2017_2019", "2017_2020", "2018_2019", "2018_2020", "2019_2020")

# Create vector of plots in the dataset
composition.plots <- sort(unique(composition.traits$PLOT))


# CALCULATE SPATIAL DISSIMILARITY MATRICES ----

# Run Bray-Curtis calculations for each year of dataset (~ 2-3 seconds)
if (run.spatial.taxonomic == TRUE){calc.beta.taxonomic.spatial(composition.years)}

# Load in results
beta.taxonomic.spatial <- read.csv(paste0("outputs/output_beta_taxonomic_spatial", filepath.brightness, filepath.37, filepath.top.hits, ".csv"))

# Run dissimilarity calculations for each year of dataset (~ 2-3 minutes)
if (run.spatial.functional == TRUE){calc.beta.functional.spatial(composition.years)}

# Load in results
beta.functional.spatial <- read.csv(paste0("outputs/output_beta_functional_spatial", filepath.brightness, filepath.37, filepath.top.hits, ".csv"))

# Run biomass distance calculations for each year of dataset (~ 2-3 seconds)
if (run.spatial.biomass == TRUE){calc.beta.biomass.spatial(composition.years)}

# Load in results
beta.biomass.spatial <- read.csv(paste0("outputs/output_beta_biomass_spatial", filepath.brightness, filepath.37, filepath.top.hits, ".csv"))

# Run spectral dissimilarity calculations for each year of dataset (~ 2-3 seconds)
if (run.spatial.spectral == TRUE){calc.beta.spectral.spatial(composition.years)}

# Load in results
beta.spectral.spatial <- read.csv(paste0("outputs/output_beta_spectral_spatial_b", buffer, filepath.brightness,
                                         filepath.smoothing, filepath.37, filepath.top.hits, ".csv"))

# Run spectral (PCA) dissimilarity calculations for each year of dataset (~ 2-3 seconds)
if (run.spatial.spectral == TRUE){calc.beta.spectral.pca.spatial(composition.years)}

# Load in results
beta.spectral.pca.spatial <- read.csv(paste0("outputs/output_beta_spectral_spatial_b", buffer, filepath.brightness,
                                             filepath.smoothing, "_PCA", filepath.37, filepath.top.hits, ".csv"))
  

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
write.csv(beta.full.spatial, file = paste0("outputs/output_beta_FULL_spatial_b", buffer, filepath.brightness,
                                           filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)


# CALCULATE TEMPORAL DISSIMILARITY MATRICES ----

# Run Bray-Curtis calculations for each plot in dataset (~ 2-3 seconds)
if (run.temporal.taxonomic == TRUE){calc.beta.taxonomic.temporal(composition.year.pairs, composition.plots)}

# Load in results
beta.taxonomic.temporal <- read.csv(paste0("outputs/output_beta_taxonomic_temporal", filepath.brightness, filepath.37, filepath.top.hits, ".csv"))

# Run dissimilarity calculations for each plot in dataset (~ 2-3 minutes)
if (run.temporal.functional == TRUE){calc.beta.functional.temporal(composition.year.pairs, composition.plots)}

# Load in results
beta.functional.temporal <- read.csv(paste0("outputs/output_beta_functional_temporal", filepath.brightness, filepath.37, filepath.top.hits, ".csv"))

# Run biomass distance calculations for each plot in dataset (~ 2-3 seconds)
if (run.temporal.biomass == TRUE){calc.beta.biomass.temporal(composition.year.pairs, composition.plots)}

# Load in results
beta.biomass.temporal <- read.csv(paste0("outputs/output_beta_biomass_temporal", filepath.brightness, filepath.37, filepath.top.hits, ".csv"))

# Run spectral dissimilarity calculations for each plot in dataset (~ 2-3 seconds)
if (run.temporal.spectral == TRUE){calc.beta.spectral.temporal(composition.year.pairs, composition.plots)}
  
# Load in results
beta.spectral.temporal <- read.csv(paste0("outputs/output_beta_spectral_temporal_b", buffer, filepath.brightness,
                                          filepath.smoothing, filepath.37, filepath.top.hits, ".csv"))


# Run spectral (PCA) dissimilarity calculations for each plot in dataset (~ 2-3 seconds)
if (run.temporal.spectral == TRUE){calc.beta.spectral.pca.temporal(composition.year.pairs, composition.plots)}
  
# Load in results
beta.spectral.pca.temporal <- read.csv(paste0("outputs/output_beta_spectral_temporal_b", buffer, filepath.brightness,
                                              filepath.smoothing, "_PCA", filepath.37, filepath.top.hits, ".csv"))


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
write.csv(beta.full.temporal, file = paste0("outputs/output_beta_FULL_temporal_b", buffer, filepath.brightness,
                                            filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)


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
ggsave(spectral.comparison.panel, filename = paste0("outputs/figures/beta_spectral_vs_PCA_comparison_b", buffer, filepath.brightness,
                                                    filepath.smoothing, filepath.37, filepath.top.hits, ".png"), width = 16, height = 7.5)


# FURTHER EXPLORATION ----

# Import dataframe of plot vegetation classes
veg.class <- read.csv("hyperspectral/data/saddgrid_npp.hh.data.csv") %>% 
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
matrices.functional.spatial <- get(load(paste0("outputs/output_beta_functional_spatial", filepath.37, filepath.top.hits, ".RData")))

# Trimmed matrices
matrix.spatial <- as.data.frame(as.matrix(matrices.functional.spatial[[1]]))[1:10, 1:10]

# Write matrix to .csv
write.csv(matrix.spatial, file = paste0("outputs/output_example_matrix_spatial", filepath.37, filepath.top.hits, ".csv"))


# EXPORT EXAMPLE MATRIX - TEMPORAL ----

# Load in temporal beta outputs
beta.temporal <- read.csv(paste0("outputs/output_beta_FULL_temporal_b", buffer, filepath.brightness,
                                 filepath.smoothing, filepath.37, filepath.top.hits, ".csv"))

# Cut to correct years
beta.temporal.years <- beta.temporal %>% 
  separate(Years, into = c("Year_1", "Year_2"), sep = "_")

# Cut to single plot (2)
beta.temporal.cut <- beta.temporal.years %>% 
  filter(PLOT == 2)

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
write.csv(matrix.temporal, file = paste0("outputs/output_example_matrix_temporal", filepath.brightness,
                                         filepath.smoothing, filepath.37, filepath.top.hits, ".csv"))

