# 08 - Functions to extract NDVI and compare to biomass
# Joseph Everest
# March 2023


# FUNCTION: EXTRACT NDVI FROM ALREADY EXTRACTED HYPERSPECTRA ----

extract.ndvi.buffer <- function(spectra){
  
  # Generate broad-band NDVI using Sentinel-2 bands
  spectra.ndvi <- spectra %>% 
    dplyr::select(Plot, Year, Wavelength, scale_Reflectance) %>% 
    mutate(Sentinel_Band = NA, # Label wavelengths in Sentinel-2 bands
           Sentinel_Band = ifelse(Wavelength > 633 & Wavelength < 695, "Band_4", Sentinel_Band),
           Sentinel_Band = ifelse(Wavelength > 726 & Wavelength < 938, "Band_8", Sentinel_Band)) %>% 
    filter(Sentinel_Band %in% c("Band_4", "Band_8")) %>% # Retain only Sentintel-2 bands
    dplyr::select(-Wavelength) %>% 
    group_by(Plot, Year, Sentinel_Band) %>% 
    summarise(mean_Reflectance = mean(scale_Reflectance)) %>% # Calculte mean reflectance for each Sentinel-2 band
    ungroup() %>% 
    pivot_wider(names_from = "Sentinel_Band", values_from = "mean_Reflectance") %>% 
    mutate(NDVI = (Band_8 - Band_4) / (Band_8 + Band_4),
           exp_NDVI = exp(NDVI)) # Add an exponentiated NDVI column
  
  # Write output to csv
  write.csv(spectra.ndvi, file = "outputs/output_saddle_ndvi.csv", row.names = FALSE)
  
  # Remove intermediate objects
  rm(spectra.ndvi)
  
  
}


# FUNCTION: NDVI 'BETA' DIVERSITY - SPATIAL ----

calc.beta.ndvi.spatial <- function(composition.years){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  ndvi.output.df <- data.frame()
  
  # Create empty list for assimilating dissimilarity matrix outputs as lists
  ndvi.output.list <- list()
  
  
  # Run a loop calculating each output per year
  for (i in composition.years){
    
    
    # Filter the input dataframe to the correct year
    ndvi.cut <- filter(saddle.ndvi, Year == i)
    
    # Create a key of the plots to row numbers
    ndvi.key <- ndvi.cut %>%
      dplyr::select(Plot) %>% 
      mutate(ID = as.numeric(row.names(.)))    
    
    # Create input dataframe for Euclidean distance calculations
    ndvi.input <- ndvi.cut %>% 
      dplyr::select(exp_NDVI)
    
    # Calculate distances
    ndvi.distance.1 <- dist(ndvi.input, method = "euclidean", diag = FALSE, upper = FALSE)
    
    # Append matrix to list for outputting
    ndvi.output.list[[paste0(i)]] <- ndvi.distance.1 
    
    # Convert to dataframe
    ndvi.distance.2 <- melt(as.matrix(ndvi.distance.1), varnames = c("row", "col")) %>% 
      rename(Dissimilarity = value)
    
    # Remove duplicate rows in dataframe (from conversion to dataframe from distance matrix)
    ndvi.distance.3 <- ndvi.distance.2 %>% 
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
    ndvi.distance.4 <- left_join(ndvi.distance.3, ndvi.key, by = c("minPlotID" = "ID")) %>% 
      rename(PLOT_1 = Plot) %>% 
      left_join(., ndvi.key, by = c("maxPlotID" = "ID")) %>%
      rename(PLOT_2 = Plot) %>% 
      dplyr::select(-c(minPlotID, maxPlotID)) %>% 
      mutate(Year = i,
             Type = "NDVI",
             Method = "Euclidean") %>% 
      relocate(Year, Type, Method, PLOT_1, PLOT_2, .before = ) %>% 
      arrange(PLOT_1, PLOT_2)
    
    # Join dataframe output to overall output dataframe
    ndvi.output.df <- rbind(ndvi.output.df, ndvi.distance.4)

    
  } # End of year if statement
  
  
  # Export dataframe of spectral distance across all years
  write.csv(ndvi.output.df, file = "outputs/output_beta_ndvi_spatial.csv", row.names = FALSE)
  
  # Export list of spectral distance across all years
  save(ndvi.output.list, file = "outputs/output_beta_ndvi_spatial.RData")
  
  
} # End of function


# FUNCTION: NDVI BETA DIVERSTIY - TEMPORAL ----

calc.beta.ndvi.temporal <- function(composition.year.pairs, composition.plots){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  ndvi.output.df <- data.frame()
  
  
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
      ndvi.input <- saddle.ndvi %>% 
        filter(Year %in% c(year.start, year.end), Plot == i) %>% 
        dplyr::select(exp_NDVI)
      
      # Calculate distances
      ndvi.distance.1 <- dist(ndvi.input, method = "euclidean", diag = FALSE, upper = FALSE)
      
      # Turn into dataframe and select required cell
      ndvi.distance.2 <- melt(as.matrix(ndvi.distance.1), varnames = c("row", "col")) %>% 
        dplyr::select(value) %>% 
        filter(value != 0) %>% 
        distinct() %>% 
        rename(Dissimilarity = value)
      
      # Add in supporting information
      ndvi.distance.3 <- ndvi.distance.2 %>% 
        mutate(Plot = i,
               Years = j,
               Type = "NDVI",
               Method = "Euclidean") %>% 
        relocate(Dissimilarity, .after = Method)
      
      # Append to main output
      ndvi.output.df <- rbind(ndvi.output.df, ndvi.distance.3)
      
      # Remove intermediate objects
      rm(ndvi.input, ndvi.distance.1, ndvi.distance.2, ndvi.distance.3)
      
      
    } # End of plot loop
    
    
  } # End of pairwise-year loop
  
  
  # Export dataframe of Bray-Curtis dissimilarity across all years
  write.csv(ndvi.output.df, file = "outputs/output_beta_ndvi_temporal.csv", row.names = FALSE)
  
  # Remove intermediate objects
  rm(ndvi.output.df)
  
  
} # End of function


# FUNCTION: SPATIAL VISUALISATIONS ----

beta.visualisations.spatial.ndvi <- function(beta.spatial.full){

  
  # Create string for plot subtitle
  plots.included <- "PLOTS: NDVI Threshold (O.2)"
  
  # Visualize relationship between ndvi and taxonomic beta diversity
  (plot.beta.spatial.n.t <- ggplot(data = beta.spatial.full, aes(x = NDVI_Dis, y = Taxonomic_Dis, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "magma", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.spatial.full$NDVI_Dis), max(beta.spatial.full$NDVI_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.spatial.full$Taxonomic_Dis), max(beta.spatial.full$Taxonomic_Dis)), expand = c(0,0)) +
      labs(title = "NDVI vs Taxonomic Dissimilarity [SPATIAL]",
           subtitle = plots.included,
           x = "NDVI (Euclidean) Dissimilarity",
           y = "Taxonomic (Bray-Curtis) Dissimilarity",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between ndvi and taxonomic beta diversity
  (plot.beta.spatial.n.f <- ggplot(data = beta.spatial.full, aes(x = NDVI_Dis, y = Functional_Dis, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "viridis", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.spatial.full$NDVI_Dis), max(beta.spatial.full$NDVI_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.spatial.full$Functional_Dis), max(beta.spatial.full$Functional_Dis)), expand = c(0,0)) +
      labs(title = "NDVI vs Functional Dissimilarity [SPATIAL]",
           subtitle = plots.included,
           x = "NDVI (Euclidean) Dissimilarity",
           y = "Functional (FDis) Dissimilarity",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between ndvi and taxonomic beta diversity
  (plot.beta.spatial.n.b <- ggplot(data = beta.spatial.full, aes(x = NDVI_Dis, y = Biomass_Dis, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.spatial.full$NDVI_Dis), max(beta.spatial.full$NDVI_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.spatial.full$Biomass_Dis), max(beta.spatial.full$Biomass_Dis)), expand = c(0,0)) +
      labs(title = "NDVI vs Biomass Dissimilarity [SPATIAL]",
           subtitle = plots.included,
           x = "NDVI (Euclidean) Dissimilarity",
           y = "Biomass (Euclidean) Dissimilarity",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Generate a panel of the outputs
  panel.spatial <- grid.arrange(plot.beta.spatial.n.t, plot.beta.spatial.n.f, plot.beta.spatial.n.b, ncol = 3)
  
  # Export plot
  ggsave(panel.spatial, file = "outputs/figures/beta_spatial_ndvi.png", width = 21, height = 5)
  
  
} # End of function


# FUNCTION: TEMPORAL VISUALISATIONS ----

beta.visualisations.temporal.ndvi <- function(beta.temporal.full){

  
  # Create string for plot subtitle
  plots.included <- "PLOTS: NDVI Threshold (O.2)"
  
  # Visualize relationship between ndvi and taxonomic beta diversity
  (plot.beta.temporal.n.t <- ggplot(data = beta.temporal.full, aes(x = NDVI_Dis, y = Taxonomic_Dis, fill = Years)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "magma", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.temporal.full$NDVI_Dis), max(beta.temporal.full$NDVI_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.temporal.full$Taxonomic_Dis), max(beta.temporal.full$Taxonomic_Dis)), expand = c(0,0)) +
      labs(title = "NDVI vs Taxonomic Dissimilarity [TEMPORAL]",
           subtitle = plots.included,
           x = "NDVI (Euclidean) Dissimilarity",
           y = "Taxonomic (Bray-Curtis) Dissimilarity",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between ndvi and taxonomic beta diversity
  (plot.beta.temporal.n.f <- ggplot(data = beta.temporal.full, aes(x = NDVI_Dis, y = Functional_Dis, fill = Years)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "viridis", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.temporal.full$NDVI_Dis), max(beta.temporal.full$NDVI_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.temporal.full$Functional_Dis), max(beta.temporal.full$Functional_Dis)), expand = c(0,0)) +
      labs(title = "NDVI vs Functional Dissimilarity [TEMPORAL]",
           subtitle = plots.included,
           x = "NDVI (Euclidean) Dissimilarity",
           y = "Functional (FDis) Dissimilarity",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between ndvi and taxonomic beta diversity
  (plot.beta.temporal.n.b <- ggplot(data = beta.temporal.full, aes(x = NDVI_Dis, y = Biomass_Dis, fill = Years)) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.temporal.full$NDVI_Dis), max(beta.temporal.full$NDVI_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.temporal.full$Biomass_Dis), max(beta.temporal.full$Biomass_Dis)), expand = c(0,0)) +
      labs(title = "NDVI vs Biomass Dissimilarity [TEMPORAL]",
           subtitle = plots.included,
           x = "NDVI (Euclidean) Dissimilarity",
           y = "Biomass (Euclidean) Dissimilarity",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Generate a panel of the outputs
  panel.temporal <- grid.arrange(plot.beta.temporal.n.t, plot.beta.temporal.n.f, plot.beta.temporal.n.b, ncol = 3)
  
  # Export plot
  ggsave(panel.temporal, file = "outputs/figures/beta_temporal_ndvi.png", width = 21, height = 5)
  
  
} # End of function


# FUNCTION: MANTEL TESTS (SPATIAL) ----

run.mantel.spatial.ndvi <- function(saddle.years){
  
  
  # Load in taxonomic, functional and biomass distance matries
  matrices.taxonomic <- get(load("outputs/output_beta_taxonomic_spatial.RData"))
  matrices.functional <- get(load("outputs/output_beta_functional_spatial.RData"))
  matrices.biomass <- get(load("outputs/output_beta_biomass_spatial.RData"))
  matrices.ndvi <- get(load("outputs/output_beta_ndvi_spatial.RData"))
  
  # Remove duplicate objects
  rm(taxonomic.output.list, functional.output.list, biomass.output.list, ndvi.output.list)
  
  # Create list to append mantel test outputs too
  mantel.tests <- list()
  
  # Create dataframe to append mantel test statistics too
  mantel.output <- data.frame()
  
  
  # Run loop to generate mantel tests for all four years across all three matrix types
  for (i in saddle.years){
    
    # Convert the year index 'i' to character
    i <- as.character(i)
    
    # Run mantel tests for all three combinations
    mantel.nt <- mantel(matrices.ndvi[[i]], matrices.taxonomic[[i]], method = "spearman", permutations = 9999, na.rm = TRUE)
    mantel.nf <- mantel(matrices.ndvi[[i]], matrices.functional[[i]], method = "spearman", permutations = 9999, na.rm = TRUE)
    mantel.nb <- mantel(matrices.ndvi[[i]], matrices.biomass[[i]], method = "spearman", permutations = 9999, na.rm = TRUE)
    
    # Save outputs to list
    mantel.tests[[paste0("NDVI_Taxonomic_", i)]] <- mantel.nt
    mantel.tests[[paste0("NDVI_Functional_", i)]] <- mantel.nf
    mantel.tests[[paste0("NDVI_Biomass_", i)]] <- mantel.nb
    
    # Extract useful values from the output
    mantel.nt.output <- data.frame("Year" = i, "Matrix_1" = "NDVI", "Matrix_2" = "Taxonomic", "R_value" = mantel.nt$statistic, "p_value" = mantel.nt$signif)
    mantel.nf.output <- data.frame("Year" = i, "Matrix_1" = "NDVI", "Matrix_2" = "Functional", "R_value" = mantel.nf$statistic, "p_value" = mantel.nf$signif)
    mantel.nb.output <- data.frame("Year" = i, "Matrix_1" = "NDVI", "Matrix_2" = "Biomass", "R_value" = mantel.nb$statistic, "p_value" = mantel.nb$signif)
    
    # Bind dataframes to output
    mantel.output <- rbind(mantel.output, mantel.nt.output, mantel.nf.output, mantel.nb.output)
    
    # Remove intermediate objects
    rm(mantel.nt, mantel.nf, mantel.nb, mantel.nt.output, mantel.nf.output, mantel.nb.output)
       
    
  } # End of year loop
  
  
  # Modify dataframe for creating correlation plots
  mantel.output.tidy <- mantel.output %>% 
    mutate(Matrix_1 = as.factor(Matrix_1),
           Matrix_2 = as.factor(Matrix_2))
  
  # Change factor order
  mantel.output.tidy$Matrix_1 <- factor(mantel.output.tidy$Matrix_1,
                                        levels = c("NDVI", "Taxonomic", "Functional", "Biomass"))
  mantel.output.tidy$Matrix_2 <- factor(mantel.output.tidy$Matrix_2,
                                        levels = c("NDVI", "Taxonomic", "Functional", "Biomass"))
  
  # Plot correlation matrices for R-values
  (mantel.plot.r <- ggplot(data = mantel.output.tidy, aes(x = Matrix_1, y = Matrix_2, fill = R_value)) +
      geom_tile() +
      scale_fill_gradient(low = "#FFFFFF", high = "#00CD00", limits = c(-0.1, 0.8)) +
      geom_text(aes(x = Matrix_1, y = Matrix_2, label = round(R_value, 4)),
                color = "black", fontface = "bold", size = 5) +
      facet_wrap(~ Year) +
      labs(title = "Mantel Test Outputs: R-values",
           subtitle = "SPATIAL",
           x = "", y = "", fill = "R-Value") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Plot correlation matrices for p-values
  (mantel.plot.p <- ggplot(data = mantel.output.tidy, aes(x = Matrix_1, y = Matrix_2, fill = p_value)) +
      geom_tile() +
      scale_fill_gradient2(low = "#F20000", high = "#FFFFFF", midpoint = 0.2, limits = c(0, 0.8)) +
      geom_text(aes(x = Matrix_1, y = Matrix_2, label = round(p_value, 4)),
                color = "black", fontface = "bold", size = 5) +
      facet_wrap(~ Year) +
      labs(title = "Mantel Test Outputs: p-values",
           subtitle = "SPATIAL",
           x = "", y = "", fill = "p-Value") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Create panel of both outputs
  mantel.panel <- grid.arrange(mantel.plot.r, mantel.plot.p, ncol = 2)
  
  # Export outputs to image file
  ggsave(mantel.panel, filename = "outputs/figures/correlations_spatial_ndvi.png", width = 20, height = 10)
  
  # Save the list of mantel test outputs
  save(mantel.tests, file = "outputs/output_statistics_spatial_ndvi.RData")
  
  # Save mantel test output statistics
  write.csv(mantel.output.tidy, file = "outputs/output_statistics_spatial_ndvi.csv", row.names = FALSE)
  
  # Remove intermediate objects
  rm(mantel.output, mantel.output.tidy, mantel.tests, mantel.plot.r, mantel.plot.p, mantel.panel)
  
  
} # End of function


# FUNCTION: MANTEL TESTS (TEMPORAL) ----

run.mantel.temporal.ndvi <- function(composition.plots, composition.years){
  
  
  # Create list to append mantel test outputs too
  mantel.tests <- list()
  
  # Create dataframe to append mantel test statistics too
  mantel.output <- data.frame()

  
  # Run loop for each plot
  for (i in composition.plots){
    
    
    # Cut the year column in two
    beta.temporal.years <- beta.temporal.full %>% 
      separate(Years, into = c("Year_1", "Year_2"), sep = "_")
    
    # Filter to single plot (would be a loop in function)
    beta.temporal.cut <- beta.temporal.years %>% 
      filter(Plot == i)
    
    # Cut to just the required distance metric
    beta.temporal.df.n <- dplyr::select(beta.temporal.cut, Year_1, Year_2, NDVI_Dis)
    beta.temporal.df.t <- dplyr::select(beta.temporal.cut, Year_1, Year_2, Taxonomic_Dis)
    beta.temporal.df.f <- dplyr::select(beta.temporal.cut, Year_1, Year_2, Functional_Dis)
    beta.temporal.df.b <- dplyr::select(beta.temporal.cut, Year_1, Year_2, Biomass_Dis)
    
    # Determine if the spectral values are all NA (NIR and NDVI mask)
    ndvi.na.sum <- sum(beta.temporal.df.n$NDVI_Dis)
    
    
    # Run an if loop for whether or not spectral is all NA
    if (!is.na(ndvi.na.sum) == TRUE){ # Not NA (i.e. has values to run)
      
      
      # Convert to a dist() class matrix (without recalculating distance)
      beta.temporal.matrix.n <- with(beta.temporal.df.n,
                                     structure(NDVI_Dis, Size = length(composition.years), Labels = composition.years,
                                               Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
      beta.temporal.matrix.t <- with(beta.temporal.df.t,
                                     structure(Taxonomic_Dis, Size = length(composition.years), Labels = composition.years,
                                               Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
      beta.temporal.matrix.f <- with(beta.temporal.df.f,
                                     structure(Functional_Dis, Size = length(composition.years), Labels = composition.years,
                                               Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
      beta.temporal.matrix.b <- with(beta.temporal.df.b,
                                     structure(Biomass_Dis, Size = length(composition.years), Labels = composition.years,
                                               Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))

      # Run mantel tests for all three combinations
      mantel.nt <- mantel(beta.temporal.matrix.n, beta.temporal.matrix.t, method = "spearman", permutations = 9999, na.rm = TRUE)
      mantel.nf <- mantel(beta.temporal.matrix.n, beta.temporal.matrix.f, method = "spearman", permutations = 9999, na.rm = TRUE)
      mantel.nb <- mantel(beta.temporal.matrix.n, beta.temporal.matrix.b, method = "spearman", permutations = 9999, na.rm = TRUE)
      
      # Save outputs to list
      mantel.tests[[paste0("NDVI_Taxonomic_", i)]] <- mantel.nt
      mantel.tests[[paste0("NDVI_Functional_", i)]] <- mantel.nf
      mantel.tests[[paste0("NDVI_Biomass_", i)]] <- mantel.nb
      
      # Extract useful values from the output
      mantel.nt.output <- data.frame("Plot" = i, "Matrix_1" = "NDVI", "Matrix_2" = "Taxonomic", "R_value" = mantel.nt$statistic, "p_value" = mantel.nt$signif)
      mantel.nf.output <- data.frame("Plot" = i, "Matrix_1" = "NDVI", "Matrix_2" = "Functional", "R_value" = mantel.nf$statistic, "p_value" = mantel.nf$signif)
      mantel.nb.output <- data.frame("Plot" = i, "Matrix_1" = "NDVI", "Matrix_2" = "Biomass", "R_value" = mantel.nb$statistic, "p_value" = mantel.nb$signif)
      
      # Bind dataframes to output
      mantel.output <- rbind(mantel.output, mantel.nt.output, mantel.nf.output, mantel.nb.output)
      
      # Remove intermediate objects
      rm(mantel.nt, mantel.nf, mantel.nb, mantel.nt.output, mantel.nf.output, mantel.nb.output)
      
      
    } else { # is NA (i.e. can't run NDVI)
      
      
      # Blank outputs
      mantel.nt.output <- data.frame("Plot" = i, "Matrix_1" = "NDVI", "Matrix_2" = "Taxonomic", "R_value" = NA, "p_value" = NA)
      mantel.nf.output <- data.frame("Plot" = i, "Matrix_1" = "NDVI", "Matrix_2" = "Functional", "R_value" = NA, "p_value" = NA)
      mantel.nb.output <- data.frame("Plot" = i, "Matrix_1" = "NDVI", "Matrix_2" = "Biomass", "R_value" = NA, "p_value" = NA)
      
      # Bind dataframes to output
      mantel.output <- rbind(mantel.output, mantel.nt.output, mantel.nf.output, mantel.nb.output)
      
      # Remove intermediate objects
      rm(mantel.nt.output, mantel.nf.output, mantel.nb.output)
      
      
    } # End of if statement for spectral NAs
    
    
    # Remove intermediate objects
    rm(beta.temporal.years, beta.temporal.cut, beta.temporal.df.n,
       beta.temporal.df.t, beta.temporal.df.f, beta.temporal.df.b)
    
    
  } # End of plot loop
  
  
  # Create tidy dataframe output with mean values
  mantel.output.tidy <- mantel.output %>% 
    mutate(Paired_Metrics = paste0(Matrix_1, ":", Matrix_2)) %>% 
    mutate(R2 = R_value^2) %>% 
    group_by(Paired_Metrics) %>% 
    mutate(mean_R = round(mean(abs(R_value), na.rm = TRUE), digits = 3),
           mean_R2 = round(mean(R2, na.rm = TRUE), digits = 3)) %>% 
    ungroup() %>% 
    dplyr::select(Paired_Metrics, mean_R, mean_R2) %>% 
    distinct()
  
  # Output temporal statistics
  write.csv(mantel.output.tidy, file = "outputs/output_statistics_temporal_ndvi.csv", row.names = FALSE)
  
  # Modify dataframes for plotting
  mantel.plotting <- mantel.output %>% 
    mutate(Paired_Metrics = paste0(Matrix_1, ":", Matrix_2)) %>% 
    group_by(Paired_Metrics) %>% 
    mutate(mean_R = round(mean(R_value, na.rm = TRUE), digits = 3)) %>% 
    ungroup() %>% 
    mutate(y_pos_box = -0.9,
           y_pos_rain = 1.2)
  
  # Generate boxplot for R values
  (mantel.boxplot.R <- ggplot(data = mantel.plotting, aes(x = Paired_Metrics, y = R_value, fill = Paired_Metrics)) +
      geom_boxplot() +
      geom_hline(aes(yintercept = 0), colour = "#000000", linetype = "dashed", size = 0.7) +
      geom_text(aes(x = Paired_Metrics, y = y_pos_box, label = paste0("Mean: ", mean_R), group = Paired_Metrics),
                size = 4.5) +
      scale_y_continuous(limits = c(-1, 1)) +
      scale_fill_viridis(option = "magma", begin = 0.3, end = 1, direction = -1, discrete = TRUE) +
      labs(title = "Mantel Test: R Values",
           subtitle = "By Plot (All Years)",
           x = "Paired Metrics",
           y = "R Value (Mantel Test)") +
      theme_1() +
      theme(legend.position = "none"))
  
  # Export boxplot of statistics
  ggsave(mantel.boxplot.R, filename = "outputs/figures/correlations_temporal_ndvi_box.png", width = 16, height = 8)
  
  # Generate raincloud plot for R values
  (mantel.raincloud.R <- ggplot(data = mantel.plotting, aes(x = Paired_Metrics, y = R_value, fill = Paired_Metrics)) +
      PupillometryR::geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
      geom_point(aes(colour = Paired_Metrics), position = position_jitter(width = 0.15), size = 1, alpha = 0.5) +
      geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
      geom_text(aes(x = Paired_Metrics, y = y_pos_rain, label = paste0("R: ", mean_R), group = Paired_Metrics),
                size = 4.5) +
      geom_hline(aes(yintercept = 0), colour = "#000000", linetype = "dashed", size = 0.7) +
      scale_y_continuous(limits = c(-1, 1.35)) +
      scale_colour_viridis(option = "magma", begin = 0.3, end = 1, direction = 1, discrete = TRUE) +
      scale_fill_viridis(option = "magma", begin = 0.3, end = 1, direction = 1, discrete = TRUE) +
      labs(title = "Mantel Test: R Values",
           subtitle = "By Plot (All Years)",
           x = "Paired Metrics",
           y = "R Value (Mantel Test)") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none"))
  
  # Export raincloud plot of statistics
  ggsave(mantel.raincloud.R, filename = "outputs/figures/correlations_temporal_ndvi_rain.png", width = 16, height = 8)
  
  
} # End of function