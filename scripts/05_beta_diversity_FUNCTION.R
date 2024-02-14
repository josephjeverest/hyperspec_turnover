# 05a - Functions for calculating Beta Diversity - Taxonomic, Functional & Spectral
# Joseph Everest
# February 2023, modified April 2023, May 2023, September 2023


# LOAD PACKAGES & THEMES ----

# Load packages
library(tidyverse)
library(viridis)
library(gridExtra)
library(reshape2)
library(vegan)
library(biotools)

# Load themes
source("hyperspectral/scripts/EX1_ggplot_themes.R")


# FUNCTION: TAXONOMIC BETA DIVERSTIY (BRAY-CURTIS) - SPATIAL ----

calc.beta.taxonomic.spatial <- function(composition.years){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  taxonomic.output.df <- data.frame()
  
  # Create empty list for assimiliating dissimilarity matrix as matrix
  taxonomic.output.list <- list()
  
  
  # Loop to sort dataframe and calculate bray-curtis dissimilarity
  for (i in composition.years){
    
    # Prepare the composition input dataframe
    composition.input.1 <- composition.traits %>% 
      filter(YEAR == i) %>% # Filter to retain single year
      dplyr::select(PLOT, SPECIES, RelativeCover) %>% # Select useful columns
      pivot_wider(names_from = "SPECIES", values_from = "RelativeCover") %>%
      mutate_all(~replace(., is.na(.), 0))
    
    # Create a key of what plot is each row number
    composition.key <- composition.input.1 %>% 
      mutate(ID = row.names(.),
             ID = as.numeric(ID)) %>% 
      dplyr::select(ID, PLOT)
    
    # Finalise the composition input dataframe
    composition.input.2 <- composition.input.1 %>% 
      dplyr::select(-PLOT)
    
    # Calculate taxonomic Bray-Curtis dissimilarity (https://rdrr.io/cran/vegan/man/vegdist.html)
    composition.bc.1 <- vegdist(composition.input.2,
                                method = "bray",
                                binary = FALSE)
    
    # Append matrix to list for outputting
    taxonomic.output.list[[paste0(i)]] <- composition.bc.1

    # Convert to dataframe
    composition.bc.2 <- melt(as.matrix(composition.bc.1), varnames = c("row", "col")) %>% 
      rename(Dissimilarity = value)
    
    # Remove duplicate rows in dataframe (from conversion to dataframe from distance matrix)
    composition.bc.3 <- composition.bc.2 %>% 
      mutate(Temp_ID = row.names(.)) %>% 
      relocate(Temp_ID, .before = ) %>% 
      pivot_longer(names_to = "PlotIDType", values_to = "PlotID", cols = c("row", "col")) %>% 
      group_by(Temp_ID) %>%
      mutate(minPlotID = min(PlotID), maxPlotID = max(PlotID)) %>% 
      ungroup() %>% 
      distinct(minPlotID, maxPlotID, .keep_all = TRUE) %>% 
      filter(minPlotID != maxPlotID) %>% 
      dplyr::select(-c(PlotIDType, PlotID, Temp_ID))
    
    # Add actual plot names back in using plot name to ID key
    composition.bc.4 <- left_join(composition.bc.3, composition.key, by = c("minPlotID" = "ID")) %>% 
      rename(PLOT_1 = PLOT) %>% 
      left_join(., composition.key, by = c("maxPlotID" = "ID")) %>%
      rename(PLOT_2 = PLOT) %>% 
      dplyr::select(-c(minPlotID, maxPlotID)) %>% 
      mutate(Year = i,
             Type = "Taxonomic",
             Method = "Bray-Curtis") %>% 
      relocate(Year, Type, Method, PLOT_1, PLOT_2, .before = ) %>% 
      arrange(PLOT_1, PLOT_2)
    
    # Join dataframe output to overall output dataframe
    taxonomic.output.df <- rbind(taxonomic.output.df, composition.bc.4)
    
    # Remove intermediate objects
    rm(composition.input.1, composition.input.2, composition.key, composition.bc.1,
       composition.bc.2, composition.bc.3, composition.bc.4)
    
    
  } # End of loop
  
  
  # Export dataframe of Bray-Curtis dissimilarity across all years
  write.csv(taxonomic.output.df, file = paste0("outputs/output_beta_taxonomic_spatial", filepath.brightness,
                                               filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Export list of Bray-Curtis dissimilarity across all years
  save(taxonomic.output.list, file = paste0("outputs/output_beta_taxonomic_spatial", filepath.brightness,
                                            filepath.smoothing, filepath.37, filepath.top.hits, ".RData"))
  
  # Remove intermediate objects
  rm(taxonomic.output.list, taxonomic.output.df)
  
  
} # End of function


# FUNCTION: TAXONOMIC BETA DIVERSTIY (BRAY-CURTIS) - TEMPORAL ----

calc.beta.taxonomic.temporal <- function(composition.year.pairs, composition.plots){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  taxonomic.output.df <- data.frame()

  
  # Run loop for each pairwise year combination
  for (j in composition.year.pairs){
    
    
    # Determine start and end year
    year.pair <- data.frame(j) %>% 
      separate(j, into = c("Start", "End"), sep = "_")
    
    # Generate string for start and end year
    year.start <- year.pair[1,1]
    year.end <- year.pair[1,2]
    
    
    # Loop to sort dataframe and calculate bray-curtis dissimilarity
    for (i in composition.plots){
      
      
      # Prepare the composition input dataframe
      composition.input <- composition.traits %>% 
        filter(YEAR %in% c(year.start, year.end),
               PLOT == i) %>% 
        dplyr::select(PlotYear, SPECIES, RelativeCover) %>% # Select useful columns
        pivot_wider(names_from = "SPECIES", values_from = "RelativeCover") %>%
        mutate_all(~replace(., is.na(.), 0)) %>% 
        dplyr::select(-PlotYear)
      
      # Calculate taxonomic Bray-Curtis dissimilarity (https://rdrr.io/cran/vegan/man/vegdist.html)
      composition.bc.1 <- vegdist(composition.input, method = "bray", binary = FALSE)
      
      # Turn into dataframe and select required cell
      composition.bc.2 <- melt(as.matrix(composition.bc.1), varnames = c("row", "col")) %>% 
        dplyr::select(value) %>% 
        filter(value != 0) %>% 
        distinct() %>% 
        rename(Dissimilarity = value)
      
      # Add in supporting information
      composition.bc.3 <- composition.bc.2 %>% 
        mutate(PLOT = i,
               Years = j,
               Type = "Taxonomic",
               Method = "Bray-Curtis") %>% 
        relocate(Dissimilarity, .after = Method)
      
      # Append to main output
      taxonomic.output.df <- rbind(taxonomic.output.df, composition.bc.3)
      
      # Remove intermediate objects
      rm(composition.input, composition.bc.1, composition.bc.2, composition.bc.3)
  
      
    } # End of plot loop
   
     
  } # End of pairwise year loop
  
  # Determine how many plots are in the original input dataframe
  plot.num.input <- length(unique(composition.traits$PLOT))
  
  # Determine how many plots are in the output dataframe
  plot.num.output <- length(unique(taxonomic.output.df$PLOT))
  
  # Determine if same number of plots
  plot.num.check <- ifelse(length(unique(composition.traits$PLOT)) == 
                             length(unique(taxonomic.output.df$PLOT)),
                           TRUE, FALSE)
  
  # Add if statement to add row for missing plot
  if (plot.num.check == FALSE){
    
    
    # Determine which plot is missing
    plot.num.missing <- as.numeric(setdiff(unique(composition.traits$PLOT),
                                           unique(taxonomic.output.df$PLOT)))
    
    # Add row for each year pair
    for (k in composition.year.pairs){
      
      
      # Create row with dissimilarity = 0 (for when plots are identical, function leaves them out)
      plot.num.input <- data.frame("PLOT" = plot.num.missing,
                                   "Years" = k,
                                   "Type" = "Taxonomic",
                                   "Method" = "Bray-Curtis",
                                   "Dissimilarity" = 0)
      
      # Append to main output
      taxonomic.output.df <- rbind(taxonomic.output.df, plot.num.input) %>% 
        arrange(PLOT)
      
      
    } # End of for loop
    
    
  } # End of if statement
  

  # Export dataframe of Bray-Curtis dissimilarity across all years
  write.csv(taxonomic.output.df, file = paste0("outputs/output_beta_taxonomic_temporal", filepath.brightness,
                                               filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Remove intermediate objects
  rm(taxonomic.output.df)
  
  
} # End of function


# FUNCTION: FUNCTIONAL BETA DIVERSITY - SPATIAL ----

  # Built on the work of Ricotta and Pavoine (2022) (https://www.sciencedirect.com/science/article/pii/S0304380022000084#sec0008)

calc.beta.functional.spatial <- function(composition.years){
  
  
  # Import functional dissimilarity functions
  source("scripts/EX1_adiv_functions.txt")
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  functional.output.df <- data.frame()
  
  # Create empty list for assimiliating dissimilarity matrix as matrix
  functional.output.list <- list()
  
  
  # Run loop to calculate for each year
  for (i in composition.years){
    
    # Filter saddle traits to correct year
    composition.traits.year <- filter(composition.traits, YEAR == i)
    
    # Modify species name so unique for each plot (allows for intraspecific trait variation)
    composition.traits.sp <- composition.traits.year %>% 
      mutate(SPECIES_unique = paste0(SPECIES, "_", PlotYear)) %>% 
      relocate(SPECIES_unique, .after = SPECIES)
    
    # Transform abundance data so plot is the row name and columns are species
    saddle.abundance <- composition.traits.sp %>% 
      dplyr::select(PLOT, SPECIES_unique, RelativeCover) %>% 
      pivot_wider(names_from = "SPECIES_unique", values_from = "RelativeCover") %>% 
      column_to_rownames(var = "PLOT") %>% # Make plot name the row names
      dplyr::select(order(colnames(.))) %>% # Order species columns alphabetically
      mutate_all(~replace(., is.na(.), 0)) # Replace all NAs with 0
    
    # Transform trait data so have a trait per species
    saddle.traits <- composition.traits.sp %>% 
      dplyr::select(SPECIES_unique, Height, SLA, Chlorophyll, D13C, D15N, LDMC, LeafN, LeafC) %>% # Don't iunclude leaf area as SLA derived from it
      distinct() %>% 
      arrange(SPECIES_unique) %>% 
      column_to_rownames(var = "SPECIES_unique")
    
    # Calculate functional distance between species using Euclidean distance between each species' traits
    saddle.fdis <- dist(scale(saddle.traits)) # Standardized to zero mean and unit standard deviation
    
    # Scale to the unit range by dividing distances by the maximum value in the distance matrix
    saddle.fdis.scaled <- saddle.fdis / max(saddle.fdis)
    
    # Run the functional dissimilarity function from Ricotta and Pavoine (2022) (https://www.sciencedirect.com/science/article/pii/S0304380022000084#sec0008)
    functional.ds.1 <- fundisparam(comm = saddle.abundance,
                                   dis = saddle.fdis.scaled,
                                   method = "D",
                                   abundance = "relative",
                                   alpha = 2,
                                   tol = 1e-8)
    
    # Append matrix to list for outputting
    functional.output.list[[paste0(i)]] <- functional.ds.1
    
    # Convert to dataframe
    functional.ds.2 <- melt(as.matrix(functional.ds.1), varnames = c("row", "col")) %>% 
      rename(Dissimilarity = value)
    
    # Remove duplicate rows in dataframe (from conversion to dataframe from distance matrix)
    functional.ds.3 <- functional.ds.2 %>% 
      mutate(Temp_ID = row.names(.)) %>% 
      relocate(Temp_ID, .before = ) %>% 
      pivot_longer(names_to = "PlotIDType", values_to = "PlotID", cols = c("row", "col")) %>% 
      group_by(Temp_ID) %>%
      mutate(minPlotID = min(PlotID), maxPlotID = max(PlotID)) %>% 
      ungroup() %>% 
      distinct(minPlotID, maxPlotID, .keep_all = TRUE) %>% 
      filter(minPlotID != maxPlotID) %>% 
      dplyr::select(minPlotID, maxPlotID, Dissimilarity)
    
    # Prepare output for combination with other dataframes
    functional.ds.4 <- functional.ds.3 %>% 
      rename(PLOT_1 = minPlotID, PLOT_2 = maxPlotID) %>% 
      mutate(Year = i,
             Type = "Functional",
             Method = "FDis") %>% 
      relocate(Year, Type, Method, PLOT_1, PLOT_2, .before = ) %>% 
      arrange(PLOT_1, PLOT_2)
    
    # Join dataframe output to overall output dataframe
    functional.output.df <- rbind(functional.output.df, functional.ds.4)
    
    # Remove intermediate objects
    rm(composition.traits.year, saddle.abundance, saddle.traits, saddle.fdis, saddle.fdis.scaled,
       functional.ds.1, functional.ds.2, functional.ds.3, functional.ds.4)
    
    
  } # End loop
  
  
  # Export dataframe of functional dissimilarity across all years
  write.csv(functional.output.df, file = paste0("outputs/output_beta_functional_spatial", filepath.brightness,
                                                filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Export list of functional dissimilarity across all years
  save(functional.output.list, file = paste0("outputs/output_beta_functional_spatial", filepath.brightness,
                                             filepath.smoothing, filepath.37, filepath.top.hits, ".RData"))
  
  # Remove intermediate objects
  rm(functional.output.list, functional.output.df)
  
  
} # End of function





# FUNCTION: FUNCTIONAL BETA DIVERSTIY - TEMPORAL ----

calc.beta.functional.temporal <- function(composition.year.pairs, composition.plots){
  
  
  # Import functional dissimilarity functions
  source("scripts/EX1_adiv_functions.txt")
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  functional.output.df <- data.frame()
  
  
  # Run loop for each pairwise year combination
  for (j in composition.year.pairs){
    
    
    # Determine start and end year
    year.pair <- data.frame(j) %>% 
      separate(j, into = c("Start", "End"), sep = "_")
    
    # Generate string for start and end year
    year.start <- year.pair[1,1]
    year.end <- year.pair[1,2]
  
  
    # Loop to sort dataframe and calculate bray-curtis dissimilarity
    for (i in composition.plots){
      
      
      # Filter saddle traits to correct plot
      composition.traits.year <- composition.traits %>% 
        filter(YEAR %in% c(year.start, year.end), PLOT == i)
      
      # Modify species name so unique for each plot (allows for intraspecific trait variation)
      composition.traits.sp <- composition.traits.year %>% 
        mutate(SPECIES_unique = paste0(SPECIES, "_", PlotYear)) %>% 
        relocate(SPECIES_unique, .after = SPECIES)
      
      
      # Run an if loop for if there is only a single identical species in the plot in each year (e.g. 37 is one shrub when top hits only)
      if (length(unique(composition.traits.sp$SPECIES)) != 1){
        
        
        # Transform abundance data so plot is the row name and columns are species
        saddle.abundance <- composition.traits.sp %>% 
          dplyr::select(YEAR, SPECIES_unique, RelativeCover) %>% 
          pivot_wider(names_from = "SPECIES_unique", values_from = "RelativeCover") %>% 
          column_to_rownames(var = "YEAR") %>% # Make plot name the row names
          dplyr::select(order(colnames(.))) %>% # Order species columns alphabetically
          mutate_all(~replace(., is.na(.), 0)) # Replace all NAs with 0
        
        # Transform trait data so have a trait per species
        saddle.traits <- composition.traits.sp %>% 
          dplyr::select(SPECIES_unique, Height, SLA, Chlorophyll, D13C, D15N, LDMC, LeafN, LeafC) %>% # Don't iunclude leaf area as SLA derived from it
          distinct() %>% 
          arrange(SPECIES_unique) %>% 
          column_to_rownames(var = "SPECIES_unique")
        
        # Calculate functional distance between species using Euclidean distance between each species' traits
        saddle.fdis <- dist(scale(saddle.traits)) # Standardized to zero mean and unit standard deviation
        
        # Scale to the unit range by dividing distances by the maximum value in the distance matrix
        saddle.fdis.scaled <- saddle.fdis / max(saddle.fdis)
        
        # Run the functional dissimilarity function from Ricotta and Pavoine (2022) (https://www.sciencedirect.com/science/article/pii/S0304380022000084#sec0008)
        functional.ds.1 <- fundisparam(comm = saddle.abundance,
                                       dis = saddle.fdis.scaled,
                                       method = "D",
                                       abundance = "relative",
                                       alpha = 2,
                                       tol = 1e-8)
        
        # Turn into dataframe and select required cell
        functional.ds.2 <- melt(as.matrix(functional.ds.1), varnames = c("row", "col")) %>% 
          dplyr::select(value) %>% 
          filter(value != 0) %>% 
          distinct() %>% 
          rename(Dissimilarity = value)
        
        # Add in supporting information
        functional.ds.3 <- functional.ds.2 %>% 
          mutate(PLOT = i,
                 Years = j,
                 Type = "Functional",
                 Method = "FDis") %>% 
          relocate(Dissimilarity, .after = Method)
        
        # Append to main output
        functional.output.df <- rbind(functional.output.df, functional.ds.3)
        
        
      } else {
        
        
        # Generate row output to append to dataframe
        functional.ds.single.species <- data.frame("PLOT" = i,
                                                   "Years" = j,
                                                   "Type" = "Functional",
                                                   "Method" = "FDis",
                                                   "Dissimilarity" = 0)
        
        # Append to main output
        functional.output.df <- rbind(functional.output.df, functional.ds.single.species)
        
        
      } # End of if else statement for single species
      
      # Remove intermediate objects
      rm(composition.traits.year, saddle.abundance, saddle.traits, saddle.fdis,
         saddle.fdis.scaled, functional.ds.1, functional.ds.2, functional.ds.3,
         functional.ds.single.species)
      
  
    } # End of plot loop
  
    
  } # End of pairwise-year loop
  
  
  # Export dataframe of Bray-Curtis dissimilarity across all years
  write.csv(functional.output.df, file = paste0("outputs/output_beta_functional_temporal", filepath.brightness,
                                                filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Remove intermediate objects
  rm(functional.output.df)
  
  
} # End of function


# FUNCTION: BIOMASS 'BETA' DIVERSITY - SPATIAL ----

calc.beta.biomass.spatial <- function(composition.years){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  biomass.output.df <- data.frame()
  
  # Create empty list for assimilating dissimilarity matrix outputs as lists
  biomass.output.list <- list()
  
  
  # Run a loop calculating each output per year
  for (i in composition.years){
    
    
    # Filter the input dataframe to the correct year
    biomass.cut <- filter(biomass, Year == i)
    
    # Create a key of the plots to row numbers
    biomass.key <- biomass.cut %>%
      dplyr::select(PLOT) %>% 
      mutate(ID = as.numeric(row.names(.)))    
    
    # Create input dataframe for Euclidean distance calculations
    biomass.input <- biomass.cut %>% 
      dplyr::select(NPP)
    
    # Calculate distances
    biomass.distance.1 <- dist(biomass.input, method = "euclidean", diag = FALSE, upper = FALSE)
    
    # Append matrix to list for outputting
    biomass.output.list[[paste0(i)]] <- biomass.distance.1 
    
    # Convert to dataframe
    biomass.distance.2 <- melt(as.matrix(biomass.distance.1), varnames = c("row", "col")) %>% 
      rename(Dissimilarity = value)
    
    # Remove duplicate rows in dataframe (from conversion to dataframe from distance matrix)
    biomass.distance.3 <- biomass.distance.2 %>% 
      mutate(Temp_ID = row.names(.)) %>% 
      relocate(Temp_ID, .before = ) %>% 
      pivot_longer(names_to = "PlotIDType", values_to = "PlotID", cols = c("row", "col")) %>% 
      group_by(Temp_ID) %>%
      mutate(minPlotID = min(PlotID), maxPlotID = max(PlotID)) %>% 
      ungroup() %>% 
      distinct(minPlotID, maxPlotID, .keep_all = TRUE) %>% 
      filter(minPlotID != maxPlotID) %>% 
      dplyr::select(minPlotID, maxPlotID, Dissimilarity)
    
    # Rename to actual plot numbers and prepare output for combination with other dataframes
    biomass.distance.4 <- left_join(biomass.distance.3, biomass.key, by = c("minPlotID" = "ID")) %>% 
      rename(PLOT_1 = PLOT) %>% 
      left_join(., biomass.key, by = c("maxPlotID" = "ID")) %>%
      rename(PLOT_2 = PLOT) %>% 
      dplyr::select(-c(minPlotID, maxPlotID)) %>% 
      mutate(Year = i,
             Type = "Biomass",
             Method = "Euclidean") %>% 
      relocate(Year, Type, Method, PLOT_1, PLOT_2, .before = ) %>% 
      arrange(PLOT_1, PLOT_2)

    # Join dataframe output to overall output dataframe
    biomass.output.df <- rbind(biomass.output.df, biomass.distance.4)
    
    
  } # End of year loop
  
  
  # Export dataframe of spectral distance across all years
  write.csv(biomass.output.df, file = paste0("outputs/output_beta_biomass_spatial", filepath.brightness,
                                             filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Export list of spectral distance across all years
  save(biomass.output.list, file = paste0("outputs/output_beta_biomass_spatial", filepath.brightness,
                                          filepath.smoothing, filepath.37, filepath.top.hits, ".RData"))
  
  
} # End of function


# FUNCTION: BIOMASS BETA DIVERSTIY - TEMPORAL ----

calc.beta.biomass.temporal <- function(composition.year.pairs, composition.plots){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  biomass.output.df <- data.frame()
  
  
  # Run loop for each pairwise year combination
  for (j in composition.year.pairs){
    
    
    # Determine start and end year
    year.pair <- data.frame(j) %>% 
      separate(j, into = c("Start", "End"), sep = "_")
    
    # Generate string for start and end year
    year.start <- year.pair[1,1]
    year.end <- year.pair[1,2]
  
  
    # Loop to sort dataframe and calculate bray-curtis dissimilarity
    for (i in composition.plots){
      
      
      # Filter the input dataframe to the correct year
      biomass.input <- biomass %>% 
        filter(Year %in% c(year.start, year.end), PLOT == i) %>% 
        dplyr::select(NPP)
      
      # Calculate distances
      biomass.distance.1 <- dist(biomass.input, method = "euclidean", diag = FALSE, upper = FALSE)
      
      # Turn into dataframe and select required cell
      biomass.distance.2 <- melt(as.matrix(biomass.distance.1), varnames = c("row", "col")) %>% 
        dplyr::select(value) %>% 
        filter(value != 0) %>% 
        distinct() %>% 
        rename(Dissimilarity = value)
  
      # Add in supporting information
      biomass.distance.3 <- biomass.distance.2 %>% 
        mutate(PLOT = i,
               Years = j,
               Type = "Biomass",
               Method = "Euclidean") %>% 
        relocate(Dissimilarity, .after = Method)
  
      # Append to main output
      biomass.output.df <- rbind(biomass.output.df, biomass.distance.3)
      
      # Remove intermediate objects
      rm(biomass.input, biomass.distance.1, biomass.distance.2, biomass.distance.3)
      
      
    } # End of plot loop
  
    
  } # End of pairwise-year loop
  
  
  # Export dataframe of Bray-Curtis dissimilarity across all years
  write.csv(biomass.output.df, file = paste0("outputs/output_beta_biomass_temporal", filepath.brightness,
                                             filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Remove intermediate objects
  rm(biomass.output.df)
  
  
} # End of function


# FUNCTION: SPECTRAL BETA DIVERSITY - SPATIAL ----

calc.beta.spectral.spatial <- function(composition.years){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  spectra.output.df <- data.frame()
  
  # Create empty lists for assimiliating dissimilarity matrices as matrix
  spectra.output.list.euc <- list()
  spectra.output.list.man <- list()
  spectra.output.list.sam <- list()
  
  
  # Run spectral euclidean distance calculations as a loop
  for (i in composition.years){
    
    
    # Filter to get correct year
    spectra.input.1 <- filter(spectra, Year == i) %>% 
      mutate(Wavelength = round(Wavelength, digits = 1)) # Round to 1.d.p. as some off by 0.0001 in 2020
    
    # Create an NIR key to know what dissimilarity values to put as NA (NIR < 0.2 likely shaded)
    spectra.key.NIR <- spectra.input.1 %>% 
      dplyr::select(PLOT, mean_NIR) %>% 
      mutate(mean_NIR = round(mean_NIR, digits = 5)) %>% # Strange error in plot 201
      arrange(PLOT) %>%
      distinct() %>% 
      mutate(mean_NIR = ifelse(mean_NIR >= 0.2, "Not_Shaded", "Shaded"))
    
    # Create an NDVI key to know what dissimilarity values to put as NA (NDVI >= 0.2 retain)
    spectra.key.NDVI <- spectra.input.1 %>% 
      dplyr::select(PLOT, NDVI_broad) %>% 
      mutate(NDVI_broad = round(NDVI_broad, digits = 5)) %>% # Strange error in plot 201
      arrange(PLOT) %>%
      distinct() %>% 
      mutate(NDVI_broad = ifelse(NDVI_broad >= 0.2, "NDVI_keep", "NDVI_remove"))
    
    # Create an input dataframe for running the euclidean distance calculations
    spectra.input.2 <- spectra.input.1 %>% 
      dplyr::select(PLOT, Band, smooth_Reflectance) %>% 
      pivot_wider(names_from = "Band", values_from = "smooth_Reflectance")
    
    # Create a plot key for rejoining plot information below
    spectra.key.plot <- spectra.input.2 %>% 
      dplyr::select(PLOT) %>% 
      mutate(ID = as.numeric(row.names(.)))
    
    # Remove plot column for calculations
    spectra.input.3 <- dplyr::select(spectra.input.2, -PLOT)
    
    # Create vector of wavelengths
    spectra.wavelengths <- sort(unique(spectra.input.1$Wavelength))
    
    # Turn into class object 'speclib' (DON'T LOAD "hsdar" PACKAGE (confuses tidyverse)
    spectra.speclib <- hsdar::speclib(as.matrix(spectra.input.3), spectra.wavelengths)
    
    # Calculate three types of distance matrix using different methods
    spectra.distance.euc.1 <- hsdar::dist.speclib(spectra.speclib, method = "euclidean") # Euclidean distance
    spectra.distance.man.1 <- hsdar::dist.speclib(spectra.speclib, method = "manhattan") # Manhattan distance
    spectra.distance.sam.1 <- hsdar::dist.speclib(spectra.speclib, method = "sam") # Spectral Angle Metic (SAM)
    
    # Append matrices to lists for outputting
    spectra.output.list.euc[[paste0(i)]] <- spectra.distance.euc.1
    spectra.output.list.man[[paste0(i)]] <- spectra.distance.man.1
    spectra.output.list.sam[[paste0(i)]] <- spectra.distance.sam.1

    # Convert matrices to dataframes
    spectra.distance.euc.2 <- melt(as.matrix(spectra.distance.euc.1), varnames = c("row", "col")) %>% rename(Euclidean = value)
    spectra.distance.man.2 <- melt(as.matrix(spectra.distance.man.1), varnames = c("row", "col")) %>% rename(Manhattan = value)
    spectra.distance.sam.2 <- melt(as.matrix(spectra.distance.sam.1), varnames = c("row", "col")) %>% rename(SAM = value)
    
    # Join dataframes together
    spectra.distance.3 <- spectra.distance.euc.2 %>% 
      left_join(., spectra.distance.man.2, by = c("row" = "row", "col" = "col")) %>% 
      left_join(., spectra.distance.sam.2, by = c("row" = "row", "col" = "col"))

    # Remove duplicate rows in dataframe (from conversion to dataframe from distance matrix)
    spectra.distance.4 <- spectra.distance.3 %>% 
      mutate(Temp_ID = row.names(.)) %>% 
      relocate(Temp_ID, .before = ) %>% 
      pivot_longer(names_to = "PlotIDType", values_to = "PlotID", cols = c("row", "col")) %>% 
      group_by(Temp_ID) %>%
      mutate(minPlotID = min(PlotID), maxPlotID = max(PlotID)) %>% 
      ungroup() %>% 
      distinct(minPlotID, maxPlotID, .keep_all = TRUE) %>% 
      filter(minPlotID != maxPlotID) %>% 
      dplyr::select(minPlotID, maxPlotID, Euclidean, Manhattan, SAM)
    
    # Rename to actual plot numbers and prepare output for combination with other dataframes
    spectra.distance.5 <- left_join(spectra.distance.4, spectra.key.plot, by = c("minPlotID" = "ID")) %>% 
      rename(PLOT_1 = PLOT) %>% 
      left_join(., spectra.key.plot, by = c("maxPlotID" = "ID")) %>%
      rename(PLOT_2 = PLOT) %>% 
      dplyr::select(-c(minPlotID, maxPlotID)) %>% 
      pivot_longer(names_to = "Method", values_to = "Dissimilarity", cols = c("Euclidean", "Manhattan", "SAM")) %>% 
      mutate(Year = i,
             Type = "Spectral") %>% 
      relocate(Year, Type, Method, PLOT_1, PLOT_2, .before = ) %>% 
      arrange(Method, PLOT_1, PLOT_2)
    
    # Replace shaded NIR values (< 0.2) with NAs
    spectra.distance.6 <- left_join(spectra.distance.5, spectra.key.NIR, by = c("PLOT_1" = "PLOT")) %>% 
      rename(Shaded_1 = mean_NIR) %>% 
      left_join(., spectra.key.NIR, by = c("PLOT_2" = "PLOT")) %>% 
      rename(Shaded_2 = mean_NIR) %>%
      mutate(Dissimilarity = ifelse(Shaded_1 == "Shaded" | Shaded_2 == "Shaded", NA, Dissimilarity)) %>% 
      dplyr::select(-c(Shaded_1, Shaded_2))
    
    # If using NDVI threshold (NDVI >= 0.2), replace plots that don't meet it with NAs
    spectra.distance.7 <- left_join(spectra.distance.6, spectra.key.NDVI, by = c("PLOT_1" = "PLOT")) %>% 
      rename(NDVI_1 = NDVI_broad) %>% 
      left_join(., spectra.key.NDVI, by = c("PLOT_2" = "PLOT")) %>% 
      rename(NDVI_2 = NDVI_broad) %>%
      mutate(Dissimilarity = ifelse(NDVI_1 == "NDVI_remove" | NDVI_2 == "NDVI_remove", NA, Dissimilarity)) %>% 
      dplyr::select(-c(NDVI_1, NDVI_2))

    # Join dataframe output to overall output dataframe
    spectra.output.df <- rbind(spectra.output.df, spectra.distance.7)
    
    # Remove intermediate objects
    rm(spectra.input.1, spectra.input.2, spectra.input.3, spectra.key.NIR, spectra.key.plot,
       spectra.speclib, spectra.distance.euc.1, spectra.distance.euc.2, spectra.distance.man.1,
       spectra.distance.man.2, spectra.distance.sam.1, spectra.distance.sam.2, spectra.distance.3,
       spectra.distance.4, spectra.distance.5, spectra.distance.6, spectra.distance.7)


  } # End of years loop

  
  # Generate wide format data for plotting
  spectra.output.df.w <- spectra.output.df %>% 
    pivot_wider(names_from = "Method", values_from = "Dissimilarity")
  
  # Generate series of plots looking at Euclidean vs Manhattan vs SAM
  (spectra.plot.e.m <- ggplot(data = spectra.output.df.w, aes(x = Euclidean, y = Manhattan, fill = Year)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "magma", begin = 0, end = 1, direction = 1, discrete = FALSE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Euclidean), max(spectra.output.df.w$Euclidean)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$Manhattan), max(spectra.output.df.w$Manhattan)), expand = c(0,0)) +
      labs(title = "Spectral Distance: Euclidean vs Manhattan",
           subtitle = "SPATIAL: All Years",
           x = "Euclidean Distance",
           y = "Manhattan Distance",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  (spectra.plot.e.s <- ggplot(data = spectra.output.df.w, aes(x = Euclidean, y = SAM, fill = Year)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "viridis", begin = 0, end = 1, direction = 1, discrete = FALSE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Euclidean), max(spectra.output.df.w$Euclidean)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$SAM), max(spectra.output.df.w$SAM)), expand = c(0,0)) +
      labs(title = "Spectral Distance: Euclidean vs SAM",
           subtitle = "SPATIAL: All Years",
           x = "Euclidean Distance",
           y = "Spectral Angle Metric",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  (spectra.plot.m.s <- ggplot(data = spectra.output.df.w, aes(x = Manhattan, y = SAM, fill = Year)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = FALSE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Manhattan), max(spectra.output.df.w$Manhattan)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$SAM), max(spectra.output.df.w$SAM)), expand = c(0,0)) +
      labs(title = "Spectral Distance: Manhattan vs SAM",
           subtitle = "SPATIAL: All Years",
           x = "Manhattan Distance",
           y = "Spectral Angle Metric",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Create a panel of the three plots
  spectra.panel <- grid.arrange(spectra.plot.e.m, spectra.plot.e.s, spectra.plot.m.s, ncol = 3)
  
  # Output plot to .png
  ggsave(spectra.panel, filename = paste0("outputs/figures/beta_spectral_comparison_spatial_b", buffer, filepath.brightness,
                                          filepath.smoothing, filepath.37, filepath.top.hits, ".png"), width = 24, height = 6.5)
  
  # Export dataframe of spectral distance across all years
  write.csv(spectra.output.df, file = paste0("outputs/output_beta_spectral_spatial_b", buffer, filepath.brightness,
                                             filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Export lists of spectral distance across all years
  save(spectra.output.list.euc, file = paste0("outputs/output_beta_spectral_spatial_b", buffer, "_euclidean", filepath.brightness,
                                              filepath.smoothing, filepath.37, filepath.top.hits, ".RData"))
  save(spectra.output.list.man, file = paste0("outputs/output_beta_spectral_spatial_b", buffer, "_manhattan", filepath.brightness,
                                              filepath.smoothing, filepath.37, filepath.top.hits, ".RData"))
  save(spectra.output.list.sam, file = paste0("outputs/output_beta_spectral_spatial_b", buffer, "_SAM", filepath.brightness,
                                              filepath.smoothing, filepath.37, filepath.top.hits, ".RData"))
  
  # Remove intermediate objects
  rm(spectra.output.df, spectra.output.list.euc, spectra.output.list.man, spectra.output.list.sam,
     spectra.output.df.w, spectra.plot.e.m, spectra.plot.e.s, spectra.plot.m.s, spectra.panel)
  
  
} # End of function


# FUNCTION: SPECTRAL BETA DIVERSTIY - TEMPORAL ----

calc.beta.spectral.temporal <- function(composition.year.pairs, composition.plots){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  spectra.output.df <- data.frame()
  
  
  # Run loop for each pairwise year combination
  for (j in composition.year.pairs){
    

    # Determine start and end year
    year.pair <- data.frame(j) %>% 
      separate(j, into = c("Start", "End"), sep = "_")
    
    # Generate string for start and end year
    year.start <- year.pair[1,1]
    year.end <- year.pair[1,2]
  
  
    # Loop to sort dataframe and calculate bray-curtis dissimilarity
    for (i in composition.plots){
      
      
      # Filter the input dataframe to the correct year
      spectra.input.1 <- spectra %>% 
        filter(Year %in% c(year.start, year.end), PLOT == i) %>% 
        mutate(Wavelength = round(Wavelength, digits = 1)) # Round to 1.d.p. as some off by 0.0001 in 2020
      
      
      # Create an NIR key to know what dissimilarity values to put as NA (NIR < 0.2 likely shaded)
      spectra.shaded <- spectra.input.1 %>% 
        dplyr::select(PLOT, mean_NIR) %>%
        mutate(mean_NIR = round(mean_NIR, digits = 5)) %>% # Strange error in plot 201
        arrange(PLOT) %>% 
        distinct() %>% 
        mutate(mean_NIR = ifelse(mean_NIR >= 0.2, "Not_Shaded", "Shaded")) %>%
        mutate(ID = paste0("Year_", row.names(.))) %>% 
        dplyr::select(-PLOT) %>% 
        pivot_wider(names_from = "ID", values_from = "mean_NIR") %>% 
        mutate(shaded = ifelse(Year_1 == "Not_Shaded" & 
                                 Year_2 == "Not_Shaded", "Not Shaded", "Shaded")) %>% 
        dplyr::select(-c(Year_1, Year_2))
      
      # Determine if either of the plots are shaded
      spectra.shaded.y.n <- unique(spectra.shaded$shaded)
      
      # Create an NDVI key to determine if have sufficient NDVI to retain
      spectra.ndvi <- spectra.input.1 %>% 
        dplyr::select(PLOT, NDVI_broad) %>% 
        mutate(NDVI_broad = round(NDVI_broad, digits = 5)) %>% # Strange error in plot 201
        arrange(PLOT) %>%
        distinct() %>% 
        mutate(NDVI_broad = ifelse(NDVI_broad >= 0.2, "NDVI_keep", "NDVI_remove")) %>% 
        mutate(ID = paste0("Year_", row.names(.))) %>% 
        dplyr::select(-PLOT) %>% 
        pivot_wider(names_from = "ID", values_from = "NDVI_broad") %>% 
        mutate(ndvi = ifelse(Year_1 == "NDVI_keep" & 
                                 Year_2 == "NDVI_keep", "NDVI_keep", "NDVI_remove")) %>% 
        dplyr::select(-c(Year_1, Year_2))
      
      # Determine if want to retain both plots based on NDVI
      spectra.ndvi.y.n <- unique(spectra.ndvi$ndvi)
  
      
      # Run an if statement to determine if both plots not shaded
      if (spectra.shaded.y.n == "Not Shaded"){
        
        
        # Create an input dataframe for running the euclidean distance calculations
        spectra.input.2 <- spectra.input.1 %>% 
          dplyr::select(Year, Band, smooth_Reflectance) %>% 
          pivot_wider(names_from = "Band", values_from = "smooth_Reflectance") %>% 
          dplyr::select(-Year)
        
        # Create vector of wavelengths
        spectra.wavelengths <- sort(unique(filter(spectra.input.1, Year == year.start)$Wavelength))
        
        # Turn into class object 'speclib' (DON'T LOAD "hsdar" PACKAGE (confuses tidyverse)
        spectra.speclib <- hsdar::speclib(as.matrix(spectra.input.2), spectra.wavelengths)
        
        # Calculate three types of distance matrix using different methods
        spectra.distance.euc.1 <- hsdar::dist.speclib(spectra.speclib, method = "euclidean") # Euclidean distance
        spectra.distance.man.1 <- hsdar::dist.speclib(spectra.speclib, method = "manhattan") # Manhattan distance
        spectra.distance.sam.1 <- hsdar::dist.speclib(spectra.speclib, method = "sam") # Spectral Angle Metic (SAM)
        
        # Turn matrices into dataframes and select required cell
        spectra.distance.euc.2 <- melt(as.matrix(spectra.distance.euc.1), varnames = c("row", "col")) %>% 
          dplyr::select(value) %>% filter(value != 0) %>% distinct() %>% rename(Euclidean = value)      
        spectra.distance.man.2 <- melt(as.matrix(spectra.distance.man.1), varnames = c("row", "col")) %>% 
          dplyr::select(value) %>% filter(value != 0) %>% distinct() %>% rename(Manhattan = value) 
        spectra.distance.sam.2 <- melt(as.matrix(spectra.distance.sam.1), varnames = c("row", "col")) %>% 
          dplyr::select(value) %>% filter(value != 0) %>% distinct() %>% rename(SAM = value) 
        
        # Add in supporting information
        spectra.distance.3 <- spectra.distance.euc.2 %>% 
          mutate(PLOT = i,
                 Years = j,
                 Type = "Spectral",
                 Manhattan = spectra.distance.man.2[1,1],
                 SAM = spectra.distance.sam.2[1,1]) %>%
          pivot_longer(names_to = "Method", values_to = "Dissimilarity", cols = c("Euclidean", "Manhattan", "SAM"))
  
        
        # If using NDVI thresholding AND NDVI is below 0.2...
        if (spectra.ndvi.y.n == "NDVI_remove"){
          
          
          # Replace all dissimilarity values for this plot combo to NA
          spectra.distance.4 <- spectra.distance.3 %>% 
            mutate(Dissimilarity = NA)
          
          
        } else { # When not using NDVI threshold and instead removing 10% bare ground plots
          
          
          # Carry over spectra.distance.3
          spectra.distance.4 <- spectra.distance.3
          
          
        } # End of NDVI threshold loop
        
        
        # Append to main output
        spectra.output.df <- rbind(spectra.output.df, spectra.distance.4)
        
        # Remove intermediate objects
        rm(spectra.input.1, spectra.shaded, spectra.shaded.y.n, spectra.input.2, spectra.wavelengths,
           spectra.speclib, spectra.distance.euc.1, spectra.distance.euc.2, spectra.distance.man.1,
           spectra.distance.man.2, spectra.distance.sam.1, spectra.distance.sam.2, spectra.distance.3,
           spectra.distance.4)
        

      } else {
        
        
        # If the plot is shaded in 2017 and/or 2020
        spectra.distance.shaded <- data.frame(PLOT = i, Years = j, Type = "Spectral",
                                              Method = "Euclidean", Dissimilarity = NA) %>% 
          add_row(PLOT = i, Years = j, Type = "Spectral",
                  Method = "Manhattan", Dissimilarity = NA) %>% 
          add_row(PLOT = i, Years = j, Type = "Spectral",
                  Method = "SAM", Dissimilarity = NA)
        
        # Append to main output
        spectra.output.df <- rbind(spectra.output.df, spectra.distance.shaded)
        
        # Remove intermediate objects
        # rm(spectra.distance.shaded)
        
        
      } # End of if statement re. shaded or not shaded
      
      
    } # End of plot loop
    
    
  } # End of pairwise-year loop
  
  
  # Export dataframe of Bray-Curtis dissimilarity across all years
  write.csv(spectra.output.df, file = paste0("outputs/output_beta_spectral_temporal_b", buffer, filepath.brightness,
                                             filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  
  # Generate wide format data for plotting
  spectra.output.df.w <- spectra.output.df %>% 
    pivot_wider(names_from = "Method", values_from = "Dissimilarity")
  
  # Generate series of plots looking at Euclidean vs Manhattan vs SAM
  (spectra.plot.e.m <- ggplot(data = spectra.output.df.w, aes(x = Euclidean, y = Manhattan)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000", fill = c("#FF3030")) +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Euclidean), max(spectra.output.df.w$Euclidean)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$Manhattan), max(spectra.output.df.w$Manhattan)), expand = c(0,0)) +
      labs(title = "Spectral Distance: Euclidean vs Manhattan",
           subtitle = "TEMPORAL: All Years",
           x = "Euclidean Distance",
           y = "Manhattan Distance") +
      theme_1() +
      theme(legend.position = "right"))
  
  (spectra.plot.e.s <- ggplot(data = spectra.output.df.w, aes(x = Euclidean, y = SAM)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000", fill = c("#00CD00")) +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Euclidean), max(spectra.output.df.w$Euclidean)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$SAM), max(spectra.output.df.w$SAM)), expand = c(0,0)) +
      labs(title = "Spectral Distance: Euclidean vs SAM",
           subtitle = "TEMPORAL: All Years",
           x = "Euclidean Distance",
           y = "Spectral Angle Metric") +
      theme_1() +
      theme(legend.position = "right"))
  
  (spectra.plot.m.s <- ggplot(data = spectra.output.df.w, aes(x = Manhattan, y = SAM)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000", fill = c("#1C86EE")) +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Manhattan), max(spectra.output.df.w$Manhattan)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$SAM), max(spectra.output.df.w$SAM)), expand = c(0,0)) +
      labs(title = "Spectral Distance: Manhattan vs SAM",
           subtitle = "TEMPORAL: All Years",
           x = "Manhattan Distance",
           y = "Spectral Angle Metric") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Create a panel of the three plots
  spectra.panel <- grid.arrange(spectra.plot.e.m, spectra.plot.e.s, spectra.plot.m.s, ncol = 3)
  
  # Output plot to .png
  ggsave(spectra.panel, filename = paste0("outputs/figures/beta_spectral_comparison_temporal_b", buffer, filepath.brightness,
                                          filepath.smoothing, filepath.37, filepath.top.hits, ".png"), width = 24, height = 6.5)
  
  # Remove intermediate objects
  rm(spectra.output.df, spectra.output.df.w, spectra.plot.e.m, spectra.plot.e.s,
     spectra.plot.m.s, spectra.panel)
  
  
} # End of function


# FUNCTION: SPECTRAL BETA DIVERSITY BY PCA - SPATIAL ----

calc.beta.spectral.pca.spatial <- function(composition.years){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  spectra.output.df <- data.frame()
  
  # Create empty lists for assimiliating dissimilarity matrices as matrix
  spectra.output.list.euc <- list()
  spectra.output.list.man <- list()
  spectra.output.list.sam <- list()
  
  
  # Run spectral euclidean distance calculations as a loop
  for (i in composition.years){
    
    
    # Filter to get correct year
    spectra.input.1 <- filter(spectra.PCA, Year == i) %>% 
      mutate(PC = str_remove(PC, pattern = "PC"),
             PC = as.numeric(PC))
    
    # Create an NIR key to know what dissimilarity values to put as NA (NIR < 0.2 likely shaded)
    spectra.key.NIR <- spectra.input.1 %>% 
      dplyr::select(PLOT, mean_NIR) %>% 
      mutate(mean_NIR = round(mean_NIR, digits = 5)) %>% # Strange error in plot 201
      arrange(PLOT) %>%
      distinct() %>% 
      mutate(mean_NIR = ifelse(mean_NIR >= 0.2, "Not_Shaded", "Shaded"))
    
    # Create an NDVI key to know what dissimilarity values to put as NA (NDVI >= 0.2 retain)
    spectra.key.NDVI <- spectra.input.1 %>% 
      dplyr::select(PLOT, NDVI_broad) %>% 
      mutate(NDVI_broad = round(NDVI_broad, digits = 5)) %>% # Strange error in plot 201
      arrange(PLOT) %>%
      distinct() %>% 
      mutate(NDVI_broad = ifelse(NDVI_broad >= 0.2, "NDVI_keep", "NDVI_remove"))
    
    # Create an input dataframe for running the euclidean distance calculations
    spectra.input.2 <- spectra.input.1 %>% 
      dplyr::select(PLOT, PC, PC_value) %>% 
      pivot_wider(names_from = "PC", values_from = "PC_value")
    
    # Create a plot key for rejoining plot information below
    spectra.key.plot <- spectra.input.2 %>% 
      dplyr::select(PLOT) %>% 
      mutate(ID = as.numeric(row.names(.)))
    
    # Remove plot column for calculations
    spectra.input.3 <- dplyr::select(spectra.input.2, -PLOT)
    
      # Create vector of wavelengths
    spectra.PCs <- sort(unique(spectra.input.1$PC))
    
    # Turn into class object 'speclib' (DON'T LOAD "hsdar" PACKAGE (confuses tidyverse)
    spectra.speclib <- hsdar::speclib(as.matrix(spectra.input.3), spectra.PCs)
    
    # Calculate three types of distance matrix using different methods
    spectra.distance.euc.1 <- hsdar::dist.speclib(spectra.speclib, method = "euclidean") # Euclidean distance
    spectra.distance.man.1 <- hsdar::dist.speclib(spectra.speclib, method = "manhattan") # Manhattan distance
    spectra.distance.sam.1 <- hsdar::dist.speclib(spectra.speclib, method = "sam") # Spectral Angle Metic (SAM)
    
    # Append matrices to lists for outputting
    spectra.output.list.euc[[paste0(i)]] <- spectra.distance.euc.1
    spectra.output.list.man[[paste0(i)]] <- spectra.distance.man.1
    spectra.output.list.sam[[paste0(i)]] <- spectra.distance.sam.1
    
    # Convert matrices to dataframes
    spectra.distance.euc.2 <- melt(as.matrix(spectra.distance.euc.1), varnames = c("row", "col")) %>% rename(Euclidean = value)
    spectra.distance.man.2 <- melt(as.matrix(spectra.distance.man.1), varnames = c("row", "col")) %>% rename(Manhattan = value)
    spectra.distance.sam.2 <- melt(as.matrix(spectra.distance.sam.1), varnames = c("row", "col")) %>% rename(SAM = value)
    
    # Join dataframes together
    spectra.distance.3 <- spectra.distance.euc.2 %>% 
      left_join(., spectra.distance.man.2, by = c("row" = "row", "col" = "col")) %>% 
      left_join(., spectra.distance.sam.2, by = c("row" = "row", "col" = "col"))
    
    # Remove duplicate rows in dataframe (from conversion to dataframe from distance matrix)
    spectra.distance.4 <- spectra.distance.3 %>% 
      mutate(Temp_ID = row.names(.)) %>% 
      relocate(Temp_ID, .before = ) %>% 
      pivot_longer(names_to = "PlotIDType", values_to = "PlotID", cols = c("row", "col")) %>% 
      group_by(Temp_ID) %>%
      mutate(minPlotID = min(PlotID), maxPlotID = max(PlotID)) %>% 
      ungroup() %>% 
      distinct(minPlotID, maxPlotID, .keep_all = TRUE) %>% 
      filter(minPlotID != maxPlotID) %>% 
      dplyr::select(minPlotID, maxPlotID, Euclidean, Manhattan, SAM)
    
    # Rename to actual plot numbers and prepare output for combination with other dataframes
    spectra.distance.5 <- left_join(spectra.distance.4, spectra.key.plot, by = c("minPlotID" = "ID")) %>% 
      rename(PLOT_1 = PLOT) %>% 
      left_join(., spectra.key.plot, by = c("maxPlotID" = "ID")) %>%
      rename(PLOT_2 = PLOT) %>% 
      dplyr::select(-c(minPlotID, maxPlotID)) %>% 
      pivot_longer(names_to = "Method", values_to = "Dissimilarity", cols = c("Euclidean", "Manhattan", "SAM")) %>% 
      mutate(Year = i,
             Type = "Spectral_PCA") %>% 
      relocate(Year, Type, Method, PLOT_1, PLOT_2, .before = ) %>% 
      arrange(Method, PLOT_1, PLOT_2)
    
    # Replace shaded NIR values (< 0.2) with NAs
    spectra.distance.6 <- left_join(spectra.distance.5, spectra.key.NIR, by = c("PLOT_1" = "PLOT")) %>% 
      rename(Shaded_1 = mean_NIR) %>% 
      left_join(., spectra.key.NIR, by = c("PLOT_2" = "PLOT")) %>% 
      rename(Shaded_2 = mean_NIR) %>%
      mutate(Dissimilarity = ifelse(Shaded_1 == "Shaded" | Shaded_2 == "Shaded", NA, Dissimilarity)) %>% 
      dplyr::select(-c(Shaded_1, Shaded_2))
    
    # If using NDVI threshold (NDVI >= 0.2), replace plots that don't meet it with NAs
    spectra.distance.7 <- left_join(spectra.distance.6, spectra.key.NDVI, by = c("PLOT_1" = "PLOT")) %>% 
      rename(NDVI_1 = NDVI_broad) %>% 
      left_join(., spectra.key.NDVI, by = c("PLOT_2" = "PLOT")) %>% 
      rename(NDVI_2 = NDVI_broad) %>%
      mutate(Dissimilarity = ifelse(NDVI_1 == "NDVI_remove" | NDVI_2 == "NDVI_remove", NA, Dissimilarity)) %>% 
      dplyr::select(-c(NDVI_1, NDVI_2))

    # Join dataframe output to overall output dataframe
    spectra.output.df <- rbind(spectra.output.df, spectra.distance.7)
    
    # Remove intermediate objects
    rm(spectra.input.1, spectra.input.2, spectra.input.3, spectra.key.NIR, spectra.key.plot,
       spectra.speclib, spectra.distance.euc.1, spectra.distance.euc.2, spectra.distance.man.1,
       spectra.distance.man.2, spectra.distance.sam.1, spectra.distance.sam.2, spectra.distance.3,
       spectra.distance.4, spectra.distance.5, spectra.distance.6, spectra.distance.7)
    
    
  } # End of years loop
  
  
  # Generate wide format data for plotting
  spectra.output.df.w <- spectra.output.df %>% 
    pivot_wider(names_from = "Method", values_from = "Dissimilarity")
  
  # Generate series of plots looking at Euclidean vs Manhattan vs SAM
  (spectra.plot.e.m <- ggplot(data = spectra.output.df.w, aes(x = Euclidean, y = Manhattan, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      
      # ggpmisc::stat_fit_glance(method = 'lm', # Don't load this package (breaks tidyverse)
      #                          method.args = list(formula = y ~ x),  geom = 'text',
      #                          aes(label = paste0("\nR2 = ", signif(..r.squared.., digits = 3)))) +
      
      scale_fill_viridis(option = "magma", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Euclidean), max(spectra.output.df.w$Euclidean)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$Manhattan), max(spectra.output.df.w$Manhattan)), expand = c(0,0)) +
      labs(# title = "Spectral PCA Distance: Euclidean vs Manhattan",
           # subtitle = "SPATIAL: All Years",
           x = "\nEuclidean Distance",
           y = "Manhattan Distance\n",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  (spectra.plot.e.s <- ggplot(data = spectra.output.df.w, aes(x = Euclidean, y = SAM, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      
      # ggpmisc::stat_fit_glance(method = 'lm', # Don't load this package (breaks tidyverse)
      #                          method.args = list(formula = y ~ x),  geom = 'text',
      #                          aes(label = paste0("\nR2 = ", signif(..r.squared.., digits = 3)))) +
      
      scale_fill_viridis(option = "viridis", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Euclidean), max(spectra.output.df.w$Euclidean)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$SAM), max(spectra.output.df.w$SAM)), expand = c(0,0)) +
      labs(# title = "Spectral PCA Distance: Euclidean vs SAM",
           # subtitle = "SPATIAL: All Years",
           x = "\nEuclidean Distance",
           y = "Spectral Angle Metric\n",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  (spectra.plot.m.s <- ggplot(data = spectra.output.df.w, aes(x = Manhattan, y = SAM, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      
      # ggpmisc::stat_fit_glance(method = 'lm', # Don't load this package (breaks tidyverse)
      #                          method.args = list(formula = y ~ x),  geom = 'text',
      #                          aes(label = paste0("\nR2 = ", signif(..r.squared.., digits = 3)))) +
      
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Manhattan), max(spectra.output.df.w$Manhattan)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$SAM), max(spectra.output.df.w$SAM)), expand = c(0,0)) +
      labs(# title = "Spectral PCA Distance: Manhattan vs SAM",
           # subtitle = "SPATIAL: All Years",
           x = "\nManhattan Distance",
           y = "Spectral Angle Metric\n",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Create a panel of the three plots
  spectra.panel <- grid.arrange(spectra.plot.e.m, spectra.plot.e.s, spectra.plot.m.s, ncol = 3)
  
  # Output plot to .png
  ggsave(spectra.panel, filename = paste0("outputs/figures/beta_spectral_comparison_spatial_b", buffer, filepath.brightness,
                                          filepath.smoothing, "_PCA", filepath.37, filepath.top.hits, ".png"), width = 18, height = 5)
  
  # Export dataframe of spectral distance across all years
  write.csv(spectra.output.df, file = paste0("outputs/output_beta_spectral_spatial_b", buffer, filepath.brightness,
                                             filepath.smoothing, "_PCA", filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Export lists of spectral distance across all years
  save(spectra.output.list.euc, file = paste0("outputs/output_beta_spectral_spatial_b", buffer, "_euclidean", filepath.brightness,
                                              filepath.smoothing, "_PCA", filepath.37, filepath.top.hits, ".RData"))
  save(spectra.output.list.man, file = paste0("outputs/output_beta_spectral_spatial_b", buffer, "_manhattan", filepath.brightness,
                                              filepath.smoothing, "_PCA", filepath.37, filepath.top.hits, ".RData"))
  save(spectra.output.list.sam, file = paste0("outputs/output_beta_spectral_spatial_b", buffer, "_SAM", filepath.brightness,
                                              filepath.smoothing, "_PCA", filepath.37, filepath.top.hits, ".RData"))
  
  # Remove intermediate objects
  rm(spectra.output.df, spectra.output.list.euc, spectra.output.list.man, spectra.output.list.sam,
     spectra.output.df.w, spectra.plot.e.m, spectra.plot.e.s, spectra.plot.m.s, spectra.panel)
  
  
} # End of function


# FUNCTION: SPECTRAL BETA DIVERSTIY BY PCA - TEMPORAL ----

calc.beta.spectral.pca.temporal <- function(composition.year.pairs, composition.plots){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  spectra.output.df <- data.frame()
  
  
  # Run loop for each pairwise year combination
  for (j in composition.year.pairs){
    
    
    # Determine start and end year
    year.pair <- data.frame(j) %>% 
      separate(j, into = c("Start", "End"), sep = "_")
    
    # Generate string for start and end year
    year.start <- year.pair[1,1]
    year.end <- year.pair[1,2]
    
    
    # Loop to sort dataframe and calculate bray-curtis dissimilarity
    for (i in composition.plots){
      
      
      # Filter the input dataframe to the correct year
      spectra.input.1 <- spectra.PCA %>% 
        filter(Year %in% c(year.start, year.end), PLOT == i) %>% 
        mutate(PC = str_remove(PC, pattern = "PC"),
               PC = as.numeric(PC))
      
      # Create an NIR key to know what dissimilarity values to put as NA (NIR < 0.2 likely shaded)
      spectra.shaded <- spectra.input.1 %>% 
        dplyr::select(PLOT, mean_NIR) %>%
        mutate(mean_NIR = round(mean_NIR, digits = 5)) %>% # Strange error in plot 201
        arrange(PLOT) %>% 
        distinct() %>% 
        mutate(mean_NIR = ifelse(mean_NIR >= 0.2, "Not_Shaded", "Shaded")) %>%
        mutate(ID = paste0("Year_", row.names(.))) %>% 
        dplyr::select(-PLOT) %>% 
        pivot_wider(names_from = "ID", values_from = "mean_NIR") %>% 
        mutate(shaded = ifelse(Year_1 == "Not_Shaded" & 
                                 Year_2 == "Not_Shaded", "Not Shaded", "Shaded")) %>% 
        dplyr::select(-c(Year_1, Year_2))
      
      # Determine if either of the plots are shaded
      spectra.shaded.y.n <- unique(spectra.shaded$shaded)
      
      # Create an NDVI key to determine if have sufficient NDVI to retain
      spectra.ndvi <- spectra.input.1 %>% 
        dplyr::select(PLOT, NDVI_broad) %>% 
        mutate(NDVI_broad = round(NDVI_broad, digits = 5)) %>% # Strange error in plot 201
        arrange(PLOT) %>%
        distinct() %>% 
        mutate(NDVI_broad = ifelse(NDVI_broad >= 0.2, "NDVI_keep", "NDVI_remove")) %>% 
        mutate(ID = paste0("Year_", row.names(.))) %>% 
        dplyr::select(-PLOT) %>% 
        pivot_wider(names_from = "ID", values_from = "NDVI_broad") %>% 
        mutate(ndvi = ifelse(Year_1 == "NDVI_keep" & 
                               Year_2 == "NDVI_keep", "NDVI_keep", "NDVI_remove")) %>% 
        dplyr::select(-c(Year_1, Year_2))
      
      # Determine if want to retain both plots based on NDVI
      spectra.ndvi.y.n <- unique(spectra.ndvi$ndvi)
      
      
      # Run an if statement to determine if both plots not shaded
      if (spectra.shaded.y.n == "Not Shaded"){
        
        
        # Create an input dataframe for running the euclidean distance calculations
        spectra.input.2 <- spectra.input.1 %>% 
          dplyr::select(Year, PC, PC_value) %>% 
          pivot_wider(names_from = "PC", values_from = "PC_value") %>% 
          dplyr::select(-Year)
        
        # Create vector of wavelengths
        spectra.PCs <- sort(unique(filter(spectra.input.1, Year == year.start)$PC))
        
        # Turn into class object 'speclib' (DON'T LOAD "hsdar" PACKAGE (confuses tidyverse)
        spectra.speclib <- hsdar::speclib(as.matrix(spectra.input.2), spectra.PCs)
        
        # Calculate three types of distance matrix using different methods
        spectra.distance.euc.1 <- hsdar::dist.speclib(spectra.speclib, method = "euclidean") # Euclidean distance
        spectra.distance.man.1 <- hsdar::dist.speclib(spectra.speclib, method = "manhattan") # Manhattan distance
        spectra.distance.sam.1 <- hsdar::dist.speclib(spectra.speclib, method = "sam") # Spectral Angle Metic (SAM)
        
        # Turn matrices into dataframes and select required cell
        spectra.distance.euc.2 <- melt(as.matrix(spectra.distance.euc.1), varnames = c("row", "col")) %>% 
          dplyr::select(value) %>% filter(value != 0) %>% distinct() %>% rename(Euclidean = value)      
        spectra.distance.man.2 <- melt(as.matrix(spectra.distance.man.1), varnames = c("row", "col")) %>% 
          dplyr::select(value) %>% filter(value != 0) %>% distinct() %>% rename(Manhattan = value) 
        spectra.distance.sam.2 <- melt(as.matrix(spectra.distance.sam.1), varnames = c("row", "col")) %>% 
          dplyr::select(value) %>% filter(value != 0) %>% distinct() %>% rename(SAM = value) 
        
        # Add in supporting information
        spectra.distance.3 <- spectra.distance.euc.2 %>% 
          mutate(PLOT = i,
                 Years = j,
                 Type = "Spectral_PCA",
                 Manhattan = spectra.distance.man.2[1,1],
                 SAM = spectra.distance.sam.2[1,1]) %>%
          pivot_longer(names_to = "Method", values_to = "Dissimilarity", cols = c("Euclidean", "Manhattan", "SAM"))
        
        
        # If using NDVI thresholding AND NDVI is below 0.2...
        if (spectra.ndvi.y.n == "NDVI_remove"){
          
          
          # Replace all dissimilarity values for this plot combo to NA
          spectra.distance.4 <- spectra.distance.3 %>% 
            mutate(Dissimilarity = NA)
          
          
        } else { # When not using NDVI threshold and instead removing 10% bare ground plots
          
          
          # Carry over spectra.distance.3
          spectra.distance.4 <- spectra.distance.3
          
          
        } # End of NDVI threshold loop
        
        
        # Append to main output
        spectra.output.df <- rbind(spectra.output.df, spectra.distance.4)
        
        # Remove intermediate objects
        rm(spectra.input.1, spectra.shaded, spectra.shaded.y.n, spectra.input.2, spectra.wavelengths,
           spectra.speclib, spectra.distance.euc.1, spectra.distance.euc.2, spectra.distance.man.1,
           spectra.distance.man.2, spectra.distance.sam.1, spectra.distance.sam.2, spectra.distance.3,
           spectra.distance.4)
        
        
      } else {
        
        
        # If the plot is shaded in 2017 and/or 2020
        spectra.distance.shaded <- data.frame(PLOT = i, Years = j, Type = "Spectral_PCA",
                                              Method = "Euclidean", Dissimilarity = NA) %>% 
          add_row(PLOT = i, Years = j, Type = "Spectral_PCA",
                  Method = "Manhattan", Dissimilarity = NA) %>% 
          add_row(PLOT = i, Years = j, Type = "Spectral_PCA",
                  Method = "SAM", Dissimilarity = NA)
        
        # Append to main output
        spectra.output.df <- rbind(spectra.output.df, spectra.distance.shaded)
        
        # Remove intermediate objects
        rm(spectra.distance.shaded)
        
        
      } # End of if statement re. shaded or not shaded
      
      
    } # End of plot loop
    
    
  } # End of pairwise-year loop
  
  
  # Export dataframe of Bray-Curtis dissimilarity across all years
  write.csv(spectra.output.df, file = paste0("outputs/output_beta_spectral_temporal_b", buffer, filepath.brightness,
                                             filepath.smoothing, "_PCA", filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  
  # Generate wide format data for plotting
  spectra.output.df.w <- spectra.output.df %>% 
    pivot_wider(names_from = "Method", values_from = "Dissimilarity")
  
  # Generate series of plots looking at Euclidean vs Manhattan vs SAM
  (spectra.plot.e.m <- ggplot(data = spectra.output.df.w, aes(x = Euclidean, y = Manhattan)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000", fill = c("#FF3030")) +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Euclidean), max(spectra.output.df.w$Euclidean)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$Manhattan), max(spectra.output.df.w$Manhattan)), expand = c(0,0)) +
      labs(title = "Spectral PCA Distance: Euclidean vs Manhattan",
           subtitle = "TEMPORAL: All Years",
           x = "Euclidean Distance",
           y = "Manhattan Distance") +
      theme_1() +
      theme(legend.position = "right"))
  
  (spectra.plot.e.s <- ggplot(data = spectra.output.df.w, aes(x = Euclidean, y = SAM)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000", fill = c("#00CD00")) +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Euclidean), max(spectra.output.df.w$Euclidean)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$SAM), max(spectra.output.df.w$SAM)), expand = c(0,0)) +
      labs(title = "Spectral PCA Distance: Euclidean vs SAM",
           subtitle = "TEMPORAL: All Years",
           x = "Euclidean Distance",
           y = "Spectral Angle Metric") +
      theme_1() +
      theme(legend.position = "right"))
  
  (spectra.plot.m.s <- ggplot(data = spectra.output.df.w, aes(x = Manhattan, y = SAM)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000", fill = c("#1C86EE")) +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_x_continuous(limits = c(min(spectra.output.df.w$Manhattan), max(spectra.output.df.w$Manhattan)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(spectra.output.df.w$SAM), max(spectra.output.df.w$SAM)), expand = c(0,0)) +
      labs(title = "Spectral PCA Distance: Manhattan vs SAM",
           subtitle = "TEMPORAL: All Years",
           x = "Manhattan Distance",
           y = "Spectral Angle Metric") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Create a panel of the three plots
  spectra.panel <- grid.arrange(spectra.plot.e.m, spectra.plot.e.s, spectra.plot.m.s, ncol = 3)
  
  # Output plot to .png
  ggsave(spectra.panel, filename = paste0("outputs/figures/beta_spectral_comparison_temporal_b", buffer, filepath.brightness,
                                          filepath.smoothing, "_PCA", filepath.37, filepath.top.hits, ".png"), width = 24, height = 6.5)
  
  # Remove intermediate objects
  rm(spectra.output.df, spectra.output.df.w, spectra.plot.e.m, spectra.plot.e.s,
     spectra.plot.m.s, spectra.panel)
  
  
} # End of function