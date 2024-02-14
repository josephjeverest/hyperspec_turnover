# 09b - Creating statistics and figures for manuscript 
# Joseph Everest
# March - May 2023


# LOAD PACKAGES & THEMES ----

# Load packages
library(tidyverse)
library(raster)
library(gridExtra)
library(ggpattern)
library(viridis)

# Load themes
source("scripts/EX1_ggplot_themes.R")
source("scripts/08_manuscript_FUNCTION.R")


# ** [1] - Decisions
brightness <- "Yes" # Default = "Yes"
smoothing <- "No" # Default = "No"
PCA <- "No" # Default = "No"
remove.37 <- "Yes" # Default = "No"
top.hits.only <- "No" # Default = "No"
buffer <- "1" # Default = "1"
exponentiate <- "No" # Default = "No"
spectral.metric = "Euclidean" # Default = "Euclidean"
single.flight.plot <- TRUE


# Generate output folder paths
if (brightness == "No"){ filepath.brightness <- "" } else { filepath.brightness <- "_brightness_normalized" }
if (smoothing == "No"){ filepath.smoothing <- "" } else {filepath.smoothing <- "_smoothed"}
if (PCA == "No"){ filepath.PCA <- "" } else {filepath.PCA <- "_PCA"}
if (remove.37 == "No"){ filepath.37 <- "" } else { filepath.37 <- "_removed_37" }
if (top.hits.only == "No"){ filepath.top.hits <- "" } else { filepath.top.hits <- "_top_hits_only" }
if (exponentiate == "No"){ filepath.exponentiate <- ""} else { filepath.exponentiate <- "_exponentiated" }


# **[CHANGE]** - DECIDE WHICH PLOTS TO GENERATE ----

# Decide on which plots to generate
run.spectral.plots <- FALSE
run.ndvi.plots <- FALSE


# [1] - DETERMINE SADDLE EXTENTS ----

# Import saddle locations
saddle.locations <- read.csv(paste0("outputs/output_saddle_locations.csv"))

# Determine range of eastings and northings
range.easting <- max(saddle.locations$EASTING) - min(saddle.locations$EASTING)
range.northing <- max(saddle.locations$NORTHING) - min(saddle.locations$NORTHING)


# [2] - DETERMINE SADDLE ELEVATIONS ----

# Re-import merged DTM as a raster
saddle.dtm <- raster("data/Niwot_DTM_2020_MERGED.tif")

# Extract elevation of each point
saddle.elevations <- raster::extract(saddle.dtm, dplyr::select(saddle.locations, EASTING, NORTHING))

# Determine minimum and maximum elevations
elevation.min <- min(saddle.elevations)
elevation.max <- max(saddle.elevations)


# [3] - SADDLE CLIMATE DATA ----

# Load in met data
saddle.metdata <- read.csv("hyperspectral/data/nwt_met-data.csv")

# Determine mean supper temp
saddle.temp <- saddle.metdata %>% 
  filter(year %in% c(2017, 2018, 2019, 2020),
         local_site == "SDL") %>%
  separate(date, into = c("y", "m", "d"), sep = "-") %>% 
  dplyr::select(y, m, airtemp_avg_homogenized) %>% 
  filter(m %in% c("06", "07", "08")) %>% 
  group_by(y) %>% 
  summarise(mean_temp = mean(airtemp_avg_homogenized)) %>% 
  ungroup()


# [4] - 'LIVING' & 'NON-VASCULAR' TOP HIT COUNTS ----

# Load in composition data
composition <- read.csv(paste0("outputs/output_saddle_composition.csv"))
  
# Convert all pins with only a single hit to 'top hits'
composition.2 <- composition %>%
  group_by(PlotYear, x, y) %>%
  mutate(num_hits = length(hit_type)) %>%
  ungroup() %>%
  mutate(hit_type = ifelse(num_hits == 1, "top", hit_type)) %>%
  dplyr::select(-num_hits)

# Remove barren plots and retain only top hits
composition.3 <- filter(composition.2, CLASS != "barren", hit_type == "top")

# Determine percentage of each plot that comprise living vascular material
composition.living.vasc <- composition.3 %>% 
  mutate(LivingVasc = ifelse(SPECIES == "NON-VASC. PLANT MATERIAL", FALSE, TRUE)) %>% 
  filter(LivingVasc == TRUE) %>% 
  group_by(PlotYear) %>% 
  summarise(percentage_LivingVasc = length(point)) %>% 
  ungroup()

# Overall mean percentage of living vasc. per plot
mean_percentage_LivingVasc <- mean(composition.living.vasc$percentage_LivingVasc)

# Determine what is listed as "NON-VASC. PLANT MATERIAL" as NOTE: not all non-vasc plant material
non.vasc.objects <- unique(filter(composition.3, SPECIES == "NON-VASC. PLANT MATERIAL")$USDA_name)

# Determine percentage of each plot that comprise living non-vascular material
composition.living.nonvasc <- composition.3 %>% 
  mutate(LivingNonVasc = ifelse(SPECIES == "NON-VASC. PLANT MATERIAL", TRUE, FALSE)) %>% 
  filter(LivingNonVasc == TRUE,
         USDA_name %in% c("Moss", "Lichen")) %>% 
  group_by(PlotYear) %>% 
  summarise(percentage_LivingNonVasc = length(point)) %>% 
  ungroup() %>% 
  left_join(composition.3, ., by = c("PlotYear" = "PlotYear")) %>% 
  dplyr::select(PlotYear, percentage_LivingNonVasc) %>% 
  distinct(PlotYear, percentage_LivingNonVasc) %>% 
  mutate(percentage_LivingNonVasc = ifelse(is.na(percentage_LivingNonVasc), 0, percentage_LivingNonVasc))

# Overall mean percentage of living non-vasc. per plot
mean_percentage_LivingNonVasc <- mean(composition.living.nonvasc$percentage_LivingNonVasc)


# [5] - LOW BIOMASS YEARS ----

# Import biomass data
biomass <- read.csv(paste0("outputs//output_biomass.csv"))

# Work out mean biomass
biomass.mean <- biomass %>% 
  group_by(Year) %>% 
  summarise(mean_NPP = mean(NPP)) %>% 
  ungroup()


# [6] - PLOT: FLIGHT AND COMPOSITION DATES ----

# Load in flights data
flight <- read.csv("data/niwot_flight_dates.csv")

# Modify flight data and generate stats
flight.stats <- flight %>% 
  rename(Year = y, Month = m, Day = d, DayOfYear = doy) %>% 
  mutate(Year = as.factor(Year)) %>% 
  group_by(Year) %>% 
  mutate(Period = max(DayOfYear) - min(DayOfYear),
         Max = max(DayOfYear),
         Min = min(DayOfYear)) %>% 
  mutate(date = as.POSIXct(date, format = "%d/%m/%Y")) %>% 
  ungroup()

# Trim flight dates to dates actually flown
if (single.flight.plot == TRUE){
  
  # Retain just the actual dates flown
  flight.stats <- flight.stats %>% 
    mutate(RETAIN = FALSE,
           RETAIN = ifelse(Year == 2017 & DayOfYear == 264, TRUE, RETAIN),
           RETAIN = ifelse(Year == 2018 & DayOfYear == 212, TRUE, RETAIN),
           RETAIN = ifelse(Year == 2019 & DayOfYear == 226, TRUE, RETAIN),
           RETAIN = ifelse(Year == 2019 & DayOfYear == 227, TRUE, RETAIN),
           RETAIN = ifelse(Year == 2020 & DayOfYear == 202, TRUE, RETAIN)) %>% 
    filter(RETAIN == TRUE) %>% 
    dplyr::select(-RETAIN) %>% 
    group_by(Year) %>% 
    mutate(Period = max(DayOfYear) - min(DayOfYear),
           Max = max(DayOfYear),
           Min = min(DayOfYear)) %>% 
    ungroup()
  
}


# Create dataframe of peak greenness stats (from phenocam dates)
evi.greenness.year <- c(2017, 2018, 2019, 2020)
evi.greenness.start <- c("04/07/2017", "20/06/2018", "14/07/2019", "28/06/2020")
evi.greenness.end <- c("05/08/2017", "17/07/2018", "07/08/2019", "27/07/2020")
evi.greenness.start.doy <- c(185, 171, 195, 180)
evi.greenness.end.doy <- c(217, 198, 219, 209)

# Create dataframe of evi greenness
peak.greenness.date <- data.frame(cbind(evi.greenness.year, evi.greenness.start, evi.greenness.end)) %>% 
  rename(Year = evi.greenness.year, Start = evi.greenness.start, End = evi.greenness.end) %>% 
  mutate(Start = as.POSIXct(Start, format = "%d/%m/%Y"),
         End = as.POSIXct(End, format = "%d/%m/%Y"),
         Year = as.numeric(Year)) %>% 
  pivot_longer(values_to = "Date", names_to = "Period", cols = 2:3)

peak.greenness.doy <- data.frame(cbind(evi.greenness.year, evi.greenness.start.doy, evi.greenness.end.doy)) %>% 
  rename(Year = evi.greenness.year, Start = evi.greenness.start.doy, End = evi.greenness.end.doy) %>% 
  pivot_longer(values_to = "DayOfYear", names_to = "Period", cols = 2:3)

peak.greenness <- left_join(peak.greenness.date, peak.greenness.doy, by = c("Year" = "Year", "Period" = "Period")) %>% 
  mutate(Window_Start = 164, Window_End = 208, Window_Middle = 186) # NEON Niwot Peak Greenness window stats

# Determine gap between end of peak greenness and start of acquisition
evi.greenness.gap.start <- c(
  filter(peak.greenness, Period == "End", Year == 2017)$DayOfYear,
  filter(peak.greenness, Period == "End", Year == 2018)$DayOfYear,
  filter(peak.greenness, Period == "End", Year == 2019)$DayOfYear,
  filter(peak.greenness, Period == "End", Year == 2020)$DayOfYear)

evi.greenness.gap.end <- c(
  unique(filter(flight.stats, Year == 2017)$Min),
  unique(filter(flight.stats, Year == 2018)$Min),
  unique(filter(flight.stats, Year == 2019)$Min),
  unique(filter(flight.stats, Year == 2020)$Min))

evi.greenness.gap <- c(
  unique(filter(flight.stats, Year == 2017)$Min) -
    filter(peak.greenness, Period == "End", Year == 2017)$DayOfYear,
  unique(filter(flight.stats, Year == 2018)$Min) -
    filter(peak.greenness, Period == "End", Year == 2018)$DayOfYear,
  unique(filter(flight.stats, Year == 2019)$Min) -
    filter(peak.greenness, Period == "End", Year == 2019)$DayOfYear,
  unique(filter(flight.stats, Year == 2020)$Min) -
    filter(peak.greenness, Period == "End", Year == 2020)$DayOfYear)

greenness.gap <- data.frame(cbind(evi.greenness.year, evi.greenness.gap.start,
                                  evi.greenness.gap.end, evi.greenness.gap)) %>% 
  rename(Year = evi.greenness.year, Start = evi.greenness.gap.start,
         End = evi.greenness.gap.end, Gap = evi.greenness.gap) %>% 
  pivot_longer(names_to = "Period", values_to = "DayOfYear", cols = 2:3) %>% 
  relocate(Gap, .after = DayOfYear) %>% 
  mutate(Year = as.factor(Year)) %>% 
  group_by(Year) %>% 
  mutate(Gap_Middle = min(DayOfYear) + (max(DayOfYear) - min(DayOfYear))/2) %>% 
  ungroup()


# Re-order year levels
flight.stats$Year <- factor(flight.stats$Year, levels = c(2020, 2019, 2018, 2017))
peak.greenness$Year <- factor(peak.greenness$Year, levels = c(2020, 2019, 2018, 2017))
greenness.gap$Year <- factor(greenness.gap$Year, levels = c(2020, 2019, 2018, 2017))

# Extract middle of window value
window.middle <- unique(peak.greenness$Window_Middle)

# Plot flight dates
(flight.plot <- ggplot() +
    
    # Rectangle for the MODIS peak greenness window
    annotate("rect", xmin = unique(peak.greenness$Window_Start),
             xmax =  unique(peak.greenness$Window_End), ymin = -Inf, ymax = Inf,
             fill = "#E4FAE4", alpha = 0.75) +
    geom_vline(aes(xintercept = unique(peak.greenness$Window_Start)), colour = "#000000", size = 0.5, linetype = "dashed") +
    geom_vline(aes(xintercept = unique(peak.greenness$Window_End)), colour = "#000000", size = 0.5, linetype = "dashed") +
    
    # Line for greenness gap between greenness window and flight line
    geom_line(data = greenness.gap, aes(x = DayOfYear, y = Year),
              size = 1, alpha = 0.75, colour = "#777777", linetype = "dashed") +
    
    # Line for phenocam peak greenness window
    geom_line(data = peak.greenness, aes(x = DayOfYear, y = Year),
              size = 1, alpha = 0.75, colour = "#72bf6a") +
    
    # Line for flight times
    geom_line(data = filter(flight.stats, DayOfYear %in% c(Max, Min)),
              aes(x = DayOfYear, y = Year, colour = Year),
              size = 1, alpha = 1) +
    
    # Triangles for phenocam peak greenness window limits
    geom_point(data = peak.greenness, aes(x = DayOfYear, y = Year),
               shape = 24, size = 4, alpha = 1, colour = "#000000", fill = "#72bf6a") +
    
    # Circles for limits of flight times
    geom_point(data = filter(flight.stats, DayOfYear %in% c(Max, Min)),
               aes(x = DayOfYear, y = Year, fill = Year),
               shape = 24, size = 5, alpha = 1, colour = "#000000") +
    
    # Labels of gap between peak greenness window and flight lines
    geom_label(data = greenness.gap, aes(x = Gap_Middle, y = Year, label = paste0(Gap)),
               label.padding = unit(0.2, "lines"),
               label.r = unit(0.15, "lines"),
               label.size = 0.25,
               size = 5) +
    
    scale_fill_viridis(option = "magma", begin = 0.5, end = 0.9, direction = 1, discrete = TRUE) +
    scale_colour_viridis(option = "magma", begin = 0.5, end = 0.9, direction = 1, discrete = TRUE) +
    labs(x = "\nDay of Year",
         y = "Acquisition Year\n") +
    theme_1() +
    theme(legend.position = "none",
          panel.border = element_rect(color = "black", fill = "transparent"),
          axis.line.x = element_blank(), 
          axis.line.y = element_blank()))

# Export plot
ggsave(flight.plot, filename = paste0("outputs/figures/manuscript/flight_dates.png"),
       width = 10, height = 5)


# [7] - LOAD IN DISTANCE MATRICES (BETA) AND MANTEL OUTPUTS (MANTEL) ----

# Load in the distance matrix values (beta)
beta.spatial <- read.csv(paste0("outputs/output_beta_FULL_spatial_b", buffer, filepath.brightness,
                                filepath.smoothing, filepath.37, filepath.top.hits, ".csv"))

beta.temporal <- read.csv(paste0("outputs/output_beta_FULL_temporal_b", buffer, filepath.brightness,
                                 filepath.smoothing, filepath.37, filepath.top.hits, ".csv"))

beta.ndvi.spatial <- read.csv(paste0("outputs/output_beta_ndvi_spatial_b", buffer, "_", spectral.metric,
                                     filepath.exponentiate, filepath.top.hits, ".csv"))

beta.ndvi.temporal <- read.csv(paste0("outputs/output_beta_ndvi_temporal_b", buffer, "_", spectral.metric,
                                      filepath.exponentiate, filepath.top.hits, ".csv"))

# Load in the mantel outputs (mantel)
mantel.spatial <- read.csv(paste0("outputs/output_statistics_spatial_b", buffer, "_", spectral.metric, filepath.brightness,
                                  filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, ".csv"))

mantel.temporal.full <- read.csv(paste0("outputs/output_statistics_temporal_FULL_b", buffer, "_", spectral.metric, filepath.brightness,
                                        filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, ".csv"))

mantel.temporal.summ <- read.csv(paste0("outputs/output_statistics_temporal_SUMMARIZED_b", buffer, "_", spectral.metric, filepath.brightness,
                                        filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, ".csv"))

mantel.spatial.ndvi <- read.csv(paste0("outputs/output_statistics_spatial_ndvi_b", buffer, "_", spectral.metric, filepath.exponentiate,
                                       filepath.37, filepath.top.hits, ".csv"))

mantel.temporal.full.ndvi <- read.csv(paste0("outputs/output_statistics_temporal_FULL_ndvi_b", buffer, "_", spectral.metric, filepath.exponentiate,
                                             filepath.37, filepath.top.hits, ".csv"))

mantel.temporal.summ.ndvi <- read.csv(paste0("outputs/output_statistics_temporal_SUMMARIZED_ndvi_b", buffer, "_", spectral.metric, filepath.exponentiate,
                                             filepath.37, filepath.top.hits, ".csv"))


# [8] - DETERMINE THE SPECTRAL METRIC BEING INVESTIGATED ----

# Determine what spectral metric we wish to analyse
selected.spectral.metric <- paste0("Spectral", filepath.PCA, "_", spectral.metric, "_Dis")


# [9] - COUNT OF PLOTS RETAINED AFTER NDVI AND NIR MASKS ----

# Spatial counts
count.spatial <- beta.spatial %>% 
  dplyr::select(Year, PLOT_1, PLOT_2, selected.spectral.metric) %>% 
  mutate(Simplified_Spectral_Euclidean_Dis = ifelse(is.na(.[,4]), TRUE, FALSE)) %>% 
  group_by(Simplified_Spectral_Euclidean_Dis) %>% 
  summarise(Count = length(Simplified_Spectral_Euclidean_Dis)) %>% 
  ungroup()

# Temporal counts
count.temporal <- beta.temporal %>% 
  dplyr::select(Years, PLOT, selected.spectral.metric) %>% 
  mutate(Simplified_Spectral_Euclidean_Dis = ifelse(is.na(.[,3]), TRUE, FALSE)) %>% 
  group_by(Simplified_Spectral_Euclidean_Dis) %>% 
  summarise(Count = length(Simplified_Spectral_Euclidean_Dis)) %>% 
  ungroup()


# [9] - PREPARE DATAFRAMES FOR PLOTTING ----

# Modify NDVI dataframe - spatial
plot.ndvi.spatial <- beta.ndvi.spatial %>% 
  dplyr::select(Year, PLOT_1, PLOT_2, Dissimilarity) %>% 
  rename(NDVI_Dis = Dissimilarity)

# Modify NDVI dataframe - temporal
plot.ndvi.temporal <- beta.ndvi.temporal %>% 
  dplyr::select(PLOT, Years, Dissimilarity) %>% 
  rename(NDVI_Dis = Dissimilarity)

# Create input dataframes for plotting
plot.spatial <- beta.spatial %>% 
  rename(Spectral_Dis = selected.spectral.metric) %>% 
  left_join(., plot.ndvi.spatial, by = c("Year" = "Year", "PLOT_1" = "PLOT_1", "PLOT_2" = "PLOT_2")) %>% 
  dplyr::select(Year, Taxonomic_Dis, Functional_Dis, Biomass_Dis, Spectral_Dis, NDVI_Dis) %>% 
  rename(Taxonomic = Taxonomic_Dis,
         Functional = Functional_Dis,
         Biomass = Biomass_Dis,
         Spectral = Spectral_Dis,
         NDVI = NDVI_Dis)

plot.temporal <- beta.temporal %>% 
  rename(Spectral_Dis = selected.spectral.metric) %>% 
  left_join(., plot.ndvi.temporal, by = c("Years" = "Years", "PLOT" = "PLOT")) %>% 
  dplyr::select(PLOT, Years, Taxonomic_Dis, Functional_Dis, Biomass_Dis, Spectral_Dis, NDVI_Dis) %>% 
  separate(Years, into = c("Year_1", "Year_2"), sep = "_") %>% 
  mutate(Year_1 = as.numeric(Year_1),
         Year_2 = as.numeric(Year_2),
         timeframe = as.character(Year_2 - Year_1)) %>% 
  rename(Taxonomic = Taxonomic_Dis,
         Functional = Functional_Dis,
         Biomass = Biomass_Dis,
         Spectral = Spectral_Dis,
         NDVI = NDVI_Dis)


# [10] - PLOT: COMPARISON OF SPATIAL HYPERSPECTRAL VS NDVI (R-VALUES) ----

# Create new dataframe with grouped variables
mantel.spatial.ALL <- mantel.spatial.ndvi %>% 
  rename(Matrix_1 = Matrix_2,
         Matrix_2 = Matrix_1) %>% 
  rbind(mantel.spatial, .) %>% 
  filter(Matrix_1 %in% c("Functional", "Taxonomic", "Biomass"),
         Matrix_2 %in% c("Spectral", "NDVI")) %>% 
  mutate(Matrix_2 = ifelse(Matrix_2 == "Spectral", "Spectral", Matrix_2)) %>% 
  arrange(Year) %>% 
  mutate(Year = as.factor(Year)) %>% 
  dplyr::select(-p_value) %>%
  pivot_wider(names_from = "Matrix_2", values_from = "R_value") %>% 
  mutate(R_diff = Spectral - NDVI,
         plot_NDVI = ifelse(R_diff >= 0, NDVI, NDVI - abs(R_diff)),
         plot_Spectral = ifelse(R_diff >= 0, Spectral, Spectral + abs(R_diff)),
         plot_Spectral_L = ifelse(R_diff >= 0, plot_Spectral, 0), # Green bars for when spectral bigger
         plot_Spectral_S = ifelse(R_diff < 0, plot_Spectral, 0)) # Red bars for when spectral smaller

# Re-order factor levels
mantel.spatial.ALL$Year <- factor(mantel.spatial.ALL$Year,
                                  levels = c(2020, 2019, 2018, 2017))
mantel.spatial.ALL$Matrix_1 <- factor(mantel.spatial.ALL$Matrix_1,
                                      levels = c("Biomass", "Taxonomic", "Functional"))


# Generate stacked bar plot
(plot.stats.all <- ggplot() +
    
    # NEW COLOURS
    geom_bar(position = "dodge", stat = "identity", # Green bars for when spectral bigger
             data = mantel.spatial.ALL,
             aes(x = Matrix_1, y = plot_Spectral_L, group = Year),
             colour = "#000000", fill = "#CD3278") +
    
    geom_bar(position = "dodge", stat = "identity", # Red bars for when spectral smaller
             data = mantel.spatial.ALL,
             aes(x = Matrix_1, y = plot_Spectral_S, group = Year),
             colour = "#000000", fill = "#F0CC68") +

    geom_bar(position = "dodge", stat = "identity", # Grey bars NDVI (except when spectral smaller, then spectral)
             data = mantel.spatial.ALL,
             aes(x = Matrix_1, y = plot_NDVI, group = Year, fill = Year),
             colour = "#000000") +
    
    scale_fill_manual(values = c("#949494", "#BEBEBE", "#D3D3D3", "#EBECF0"), # Filling grey bars
                      guide = guide_legend(reverse = TRUE)) +
    
    labs(x = expression(beta*" Metric"),
         y = "\nMantel R Statistic",
         fill = "Year") +
    coord_flip() +
    theme_1() +
    theme(legend.position = "right",
          plot.caption = element_text(hjust = 0, face = "plain"))
)

# Export plot to .png
ggsave(plot.stats.all, filename = paste0("outputs/figures/manuscript/spectral_vs_ndvi_b", buffer, "_", spectral.metric, filepath.exponentiate,
                                         filepath.brightness, filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, ".png"), width = 10, height = 6)

  
# [11] - PLOTS: SPECTRAL VS ALL OTHER METRICS - SPATIAL & TEMPORAL ----

# Generate final plots comparing spectral and all other disismilarity
if(run.spectral.plots == TRUE){
  
  # Run loop over plotting function for Functional/Taxonomic/Biomass and Spectral
  for (a in c("Functional", "Taxonomic", "Biomass")){
    
    # Run function
    create.spatial.temporal.plots.spectral(SpectralMetric = "Spectral",
                                           BetaMetric = a,
                                           SpatialData = plot.spatial,
                                           TemporalData = mantel.temporal.full,
                                           SpatialStatisticsData = mantel.spatial,
                                           TemporalStatisticsData = mantel.temporal.summ)
    
  } # End of function loop
  
} # End of if statement


# [12] - PLOTS: NDVI VS ALL OTHER METRICS - SPATIAL & TEMPORAL ----

# Generate final plots comparing spectral and all other disismilarity
if(run.ndvi.plots == TRUE){
  
  # Run loop over plotting function for Functional/Taxonomic/Biomass and Spectral
  for (a in c("Functional", "Taxonomic", "Biomass")){
    
    # Run function
    create.spatial.temporal.plots.NDVI(SpectralMetric = "NDVI",
                                       BetaMetric = a,
                                       SpatialData = plot.spatial,
                                       TemporalData = mantel.temporal.full.ndvi,
                                       SpatialStatisticsData = mantel.spatial.ndvi,
                                       TemporalStatisticsData = mantel.temporal.summ.ndvi)
    
  } # End of function loop
  
} # End of if statement


# [13] - SPECTRAL DISTANCES BY YEAR-PAIR ----

# BOXPLOT: SPECTRAL DISTANCE BY YEAR ----

# Create a dataframe for plotting spectral distances by year pair
beta.temporal.boxplot <- beta.temporal %>% 
  mutate(Years = str_replace_all(Years, pattern = "_", replacement = " - "))

# Create a boxplot of spectral distance by year pairs
(boxplot.spectral <- ggplot(data = beta.temporal.boxplot,
                            aes(x = as.factor(Years), y = Spectral_Euclidean_Dis, fill = Years)) +
   geom_boxplot() +
   scale_fill_viridis(option = "magma", begin = 0.5, end = 1, direction = 1, discrete = TRUE) +
   # scale_y_continuous(limits = c(0, 35), expand = c(0,0)) +
   labs(x = "\nYear Pair",
        y = "Spectral\n(Euclidean) Distance\n") +
   theme_1() +
   theme(legend.position = "none"))

# Export the plot
ggsave(boxplot.spectral, filename = "outputs/figures/manuscript/spectral_distances_year_pair.png",
       width = 12, height = 6)
