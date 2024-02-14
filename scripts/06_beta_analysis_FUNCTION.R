# 07a - Functions for statistically analysing output beta statistics
# Joseph Everest
# February 2023


# LOAD PACKAGES & THEMES----

# Load packages
library(tidyverse)
library(gridExtra)
library(vegan)

# Load themes
source("scripts/EX1_ggplot_themes.R")


# FUNCTION: MANTEL TESTS (SPATIAL) ----

run.mantel.spatial <- function(saddle.years){
  
  
  # Load in taxonomic, functional and biomass distance matries
  matrices.taxonomic <- get(load(paste0("outputs/output_beta_taxonomic_spatial", filepath.brightness, filepath.37, filepath.top.hits, ".RData")))
  matrices.functional <- get(load(paste0("outputs/output_beta_functional_spatial", filepath.brightness, filepath.37, filepath.top.hits, ".RData")))
  matrices.biomass <- get(load(paste0("outputs/output_beta_biomass_spatial", filepath.brightness, filepath.37, filepath.top.hits, ".RData")))
  
  # Remove duplicate objects
  rm(taxonomic.output.list, functional.output.list, biomass.output.list)
  
  # Run if statement to determine which spectral matrices to import
  if (PCA == "No"){
    
    # Load in matrices not PCAed
    matrices.spectral <- get(load(paste0("outputs/output_beta_spectral_spatial_b", buffer, "_", spectral.metric, filepath.brightness,
                                         filepath.smoothing, filepath.37, filepath.top.hits, ".RData")))
    
  } else {
    
    # Load in Euclidean matrices that are PCAed
    matrices.spectral <- get(load(paste0("outputs/output_beta_spectral_spatial_b", buffer, "_", spectral.metric, filepath.brightness,
                                         filepath.smoothing, "_PCA", filepath.37, filepath.top.hits, ".RData")))
    
  }

  
  # Create list to append mantel test outputs too
  mantel.tests <- list()
  
  # Create dataframe to append mantel test statistics too
  mantel.output <- data.frame()
  
  
  # Run loop to generate mantel tests for all four years across all three matrix types
  for (i in saddle.years){
    
    # Convert the year index 'i' to character
    i <- as.character(i)
    
    # Run mantel tests for all six combinations
    mantel.tf <- mantel(matrices.taxonomic[[i]], matrices.functional[[i]], method = "pearson", permutations = 9999, na.rm = TRUE)
    mantel.tb <- mantel(matrices.taxonomic[[i]], matrices.biomass[[i]], method = "pearson", permutations = 9999, na.rm = TRUE)
    mantel.ts <- mantel(matrices.taxonomic[[i]], matrices.spectral[[i]], method = "pearson", permutations = 9999, na.rm = TRUE)
    mantel.fb <- mantel(matrices.functional[[i]], matrices.biomass[[i]], method = "pearson", permutations = 9999, na.rm = TRUE)
    mantel.fs <- mantel(matrices.functional[[i]], matrices.spectral[[i]], method = "pearson", permutations = 9999, na.rm = TRUE)
    mantel.bs <- mantel(matrices.biomass[[i]], matrices.spectral[[i]], method = "pearson", permutations = 9999, na.rm = TRUE)
    
    # Save outputs to list
    mantel.tests[[paste0("Taxonomic_Functional_", i)]] <- mantel.tf
    mantel.tests[[paste0("Taxonomic_Biomass_", i)]] <- mantel.tb
    mantel.tests[[paste0("Taxonomic_Spectral", filepath.PCA, "_", i)]] <- mantel.ts
    mantel.tests[[paste0("Functional_Biomass_", i)]] <- mantel.fb
    mantel.tests[[paste0("Functional_Spectral", filepath.PCA, "_", i)]] <- mantel.fs
    mantel.tests[[paste0("Biomass_Spectral", filepath.PCA, "_", i)]] <- mantel.bs
    
    # Extract useful values from the output
    mantel.tf.output <- data.frame("Year" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = "Functional", "R_value" = mantel.tf$statistic, "p_value" = mantel.tf$signif)
    mantel.tb.output <- data.frame("Year" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = "Biomass", "R_value" = mantel.tb$statistic, "p_value" = mantel.tb$signif)
    mantel.ts.output <- data.frame("Year" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = mantel.ts$statistic, "p_value" = mantel.ts$signif)
    mantel.fb.output <- data.frame("Year" = i, "Matrix_1" = "Functional", "Matrix_2" = "Biomass", "R_value" = mantel.fb$statistic, "p_value" = mantel.fb$signif)
    mantel.fs.output <- data.frame("Year" = i, "Matrix_1" = "Functional", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = mantel.fs$statistic, "p_value" = mantel.fs$signif)
    mantel.bs.output <- data.frame("Year" = i, "Matrix_1" = "Biomass", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = mantel.bs$statistic, "p_value" = mantel.bs$signif)
    
    # Bind dataframes to output
    mantel.output <- rbind(mantel.output, mantel.tf.output, mantel.tb.output, mantel.ts.output,
                           mantel.fb.output, mantel.fs.output, mantel.bs.output)

    # Remove intermediate objects
    rm(mantel.tf, mantel.tb, mantel.ts, mantel.fb, mantel.fs, mantel.bs, 
       mantel.tf.output, mantel.tb.output, mantel.ts.output,
       mantel.fb.output, mantel.fs.output, mantel.bs.output)
    
    
  } # End of year loop
  
  
  # Modify dataframe for creating correlation plots
  mantel.output.tidy <- mantel.output %>% 
    mutate(Matrix_1 = as.factor(Matrix_1),
           Matrix_2 = as.factor(Matrix_2))
  
  # Change factor order
  mantel.output.tidy$Matrix_1 <- factor(mantel.output.tidy$Matrix_1,
                                        levels = c("Taxonomic", "Functional", "Biomass", paste0("Spectral", filepath.PCA)))
  mantel.output.tidy$Matrix_2 <- factor(mantel.output.tidy$Matrix_2,
                                        levels = c("Taxonomic", "Functional", "Biomass", paste0("Spectral", filepath.PCA)))
  
  # Plot correlation matrices for R-values
  (mantel.plot.r <- ggplot(data = mantel.output.tidy, aes(x = Matrix_1, y = Matrix_2, fill = R_value)) +
      geom_tile() +
      scale_fill_gradient(low = "#FFFFFF", high = "#00CD00", limits = c(-0.1, 0.8)) +
      geom_text(aes(x = Matrix_1, y = Matrix_2, label = round(R_value, 4)),
                color = "black", fontface = "bold", size = 5) +
      facet_wrap(~ Year) +
      labs(title = "Mantel Test Outputs: R-values",
           subtitle = paste0("SPATIAL: (Spectral", filepath.PCA, "as ", spectral.metric, ")"),
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
           subtitle = paste0("SPATIAL: (Spectral", filepath.PCA, "as ", spectral.metric, ")"),
           x = "", y = "", fill = "p-Value") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Create panel of both outputs
  mantel.panel <- grid.arrange(mantel.plot.r, mantel.plot.p, ncol = 2)
  
  # Export outputs to image file
  ggsave(mantel.panel, filename = paste0("outputs/figures/correlations_spatial_b", buffer, "_", spectral.metric, filepath.brightness,
                                         filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, ".png"), width = 20, height = 10)
  
  # Save the list of mantel test outputs
  save(mantel.tests, file = paste0("outputs/output_statistics_spatial_b", buffer, "_", spectral.metric, filepath.brightness,
                                   filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, ".RData"))
  
  # Save mantel test output statistics
  write.csv(mantel.output.tidy, file = paste0("outputs/output_statistics_spatial_b", buffer, "_", spectral.metric, filepath.brightness,
                                              filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Remove intermediate objects
  rm(mantel.output, mantel.output.tidy, mantel.tests, mantel.plot.r, mantel.plot.p, mantel.panel)
  
  
} # End of function


# FUNCTION: MANTEL TESTS (TEMPORAL) ----

run.mantel.temporal <- function(saddle.plots){
  
   # Create list to append mantel test outputs too
  mantel.tests <- list()
  
  # Create dataframe to append mantel test statistics too
  mantel.output <- data.frame()
  
  # Generate column name for spectral variable we wish to keep
  spectral.variable <- paste0("Spectral", filepath.PCA, "_", spectral.metric, "_Dis")
  
  
  # Run loop for each plot
  for (i in saddle.plots){
    
    
    # Cut the year column in two
    beta.temporal.years <- beta.temporal %>% 
      separate(Years, into = c("Year_1", "Year_2"), sep = "_")
    
    # Filter to single plot (would be a loop in function)
    beta.temporal.cut <- beta.temporal.years %>% 
      filter(PLOT == i)
    
    # Cut to just the required distance metric
    beta.temporal.df.t <- dplyr::select(beta.temporal.cut, Year_1, Year_2, Taxonomic_Dis)
    beta.temporal.df.f <- dplyr::select(beta.temporal.cut, Year_1, Year_2, Functional_Dis)
    beta.temporal.df.b <- dplyr::select(beta.temporal.cut, Year_1, Year_2, Biomass_Dis)
    beta.temporal.df.s <- dplyr::select(beta.temporal.cut, Year_1, Year_2, spectral.variable) %>% 
      rename(Spectral_Dis = paste0(spectral.variable))
    
    # Determine if the spectral values are all NA (NIR and NDVI mask)
    spectral.na.sum <- sum(beta.temporal.df.s$Spectral_Dis)
    
    
    # Run an if loop for whether or not any of the spectral values are NA
    if (!is.na(spectral.na.sum) == TRUE){ # Not NA (i.e. has all values to run)
      
      
      # Determine whether the temporal matrices can run (e.g. no temporal matrices for plot 37 as only one species in each year)
      distance.values.sum <- sum(beta.temporal.df.t$Taxonomic_Dis, na.rm = TRUE) 
      
      
      # Run if loop to exclude plots where no distance values
      if (distance.values.sum > 0){
        
        
        # Convert to a dist() class matrix (without recalculating distance)
        beta.temporal.matrix.t <- with(beta.temporal.df.t,
                                       structure(Taxonomic_Dis, Size = length(saddle.years), Labels = saddle.years,
                                                 Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
        beta.temporal.matrix.f <- with(beta.temporal.df.f,
                                       structure(Functional_Dis, Size = length(saddle.years), Labels = saddle.years,
                                                 Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
        beta.temporal.matrix.b <- with(beta.temporal.df.b,
                                       structure(Biomass_Dis, Size = length(saddle.years), Labels = saddle.years,
                                                 Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
        beta.temporal.matrix.s <- with(beta.temporal.df.s,
                                       structure(Spectral_Dis, Size = length(saddle.years), Labels = saddle.years,
                                                 Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
        
        # Run mantel tests on each two-metric permutation
        mantel.tf <- mantel(beta.temporal.matrix.t, beta.temporal.matrix.f, method = "pearson", permutations = 9999, na.rm = TRUE)
        mantel.tb <- mantel(beta.temporal.matrix.t, beta.temporal.matrix.b, method = "pearson", permutations = 9999, na.rm = TRUE)
        mantel.ts <- mantel(beta.temporal.matrix.t, beta.temporal.matrix.s, method = "pearson", permutations = 9999, na.rm = TRUE)
        mantel.fb <- mantel(beta.temporal.matrix.f, beta.temporal.matrix.b, method = "pearson", permutations = 9999, na.rm = TRUE)
        mantel.fs <- mantel(beta.temporal.matrix.f, beta.temporal.matrix.s, method = "pearson", permutations = 9999, na.rm = TRUE)
        mantel.bs <- mantel(beta.temporal.matrix.b, beta.temporal.matrix.s, method = "pearson", permutations = 9999, na.rm = TRUE)
        
        # Save outputs to list
        mantel.tests[[paste0("Taxonomic_Functional_", i)]] <- mantel.tf
        mantel.tests[[paste0("Taxonomic_Biomass_", i)]] <- mantel.tb
        mantel.tests[[paste0("Taxonomic_Spectral", filepath.PCA, "_", i)]] <- mantel.ts
        mantel.tests[[paste0("Functional_Biomass_", i)]] <- mantel.fb
        mantel.tests[[paste0("Functional_Spectral", filepath.PCA, "_", i)]] <- mantel.fs
        mantel.tests[[paste0("Biomass_Spectral", filepath.PCA, "_", i)]] <- mantel.bs
        
        # Extract useful values from the output
        mantel.tf.output <- data.frame("PLOT" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = "Functional", "R_value" = mantel.tf$statistic, "p_value" = mantel.tf$signif)
        mantel.tb.output <- data.frame("PLOT" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = "Biomass", "R_value" = mantel.tb$statistic, "p_value" = mantel.tb$signif)
        mantel.ts.output <- data.frame("PLOT" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = mantel.ts$statistic, "p_value" = mantel.ts$signif)
        mantel.fb.output <- data.frame("PLOT" = i, "Matrix_1" = "Functional", "Matrix_2" = "Biomass", "R_value" = mantel.fb$statistic, "p_value" = mantel.fb$signif)
        mantel.fs.output <- data.frame("PLOT" = i, "Matrix_1" = "Functional", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = mantel.fs$statistic, "p_value" = mantel.fs$signif)
        mantel.bs.output <- data.frame("PLOT" = i, "Matrix_1" = "Biomass", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = mantel.bs$statistic, "p_value" = mantel.bs$signif)
        
        # Bind dataframes to output
        mantel.output <- rbind(mantel.output, mantel.tf.output, mantel.tb.output, mantel.ts.output,
                               mantel.fb.output, mantel.fs.output, mantel.bs.output)
        
        
      } else {
        
        
        # Create template outputs to bind to full output
        mantel.tf.output <- data.frame("PLOT" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = "Functional", "R_value" = NA, "p_value" = NA)
        mantel.tb.output <- data.frame("PLOT" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = "Biomass", "R_value" = NA, "p_value" = NA)
        mantel.ts.output <- data.frame("PLOT" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = NA, "p_value" = NA)
        mantel.fb.output <- data.frame("PLOT" = i, "Matrix_1" = "Functional", "Matrix_2" = "Biomass", "R_value" = NA, "p_value" = NA)
        mantel.fs.output <- data.frame("PLOT" = i, "Matrix_1" = "Functional", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = NA, "p_value" = NA)
        mantel.bs.output <- data.frame("PLOT" = i, "Matrix_1" = "Biomass", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = NA, "p_value" = NA)
        
        # Bind dataframes to output
        mantel.output <- rbind(mantel.output, mantel.tf.output, mantel.tb.output, mantel.ts.output,
                               mantel.fb.output, mantel.fs.output, mantel.bs.output)
        
        
      } # End of if statement re matrixes populated with NAs (e.g. plot 37)
      
      
            # Remove intermediate objects
      rm(beta.temporal.matrix.t, beta.temporal.matrix.f, beta.temporal.matrix.b,
         beta.temporal.matrix.s, mantel.tf, mantel.tb, mantel.ts, mantel.fb,
         mantel.fs, mantel.bs, mantel.tf.output, mantel.tb.output, mantel.ts.output,
         mantel.fb.output, mantel.fs.output, mantel.bs.output)
      
      
    } else { # is NA (i.e. can't run spectral)
      
      
      # Convert to a dist() class matrix (without recalculating distance)
      beta.temporal.matrix.t <- with(beta.temporal.df.t,
                                     structure(Taxonomic_Dis, Size = length(saddle.years), Labels = saddle.years,
                                               Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
      beta.temporal.matrix.f <- with(beta.temporal.df.f,
                                     structure(Functional_Dis, Size = length(saddle.years), Labels = saddle.years,
                                               Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
      beta.temporal.matrix.b <- with(beta.temporal.df.b,
                                     structure(Biomass_Dis, Size = length(saddle.years), Labels = saddle.years,
                                               Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
      
      # Run mantel tests on each two-metric permutation
      mantel.tf <- mantel(beta.temporal.matrix.t, beta.temporal.matrix.f, method = "pearson", permutations = 9999, na.rm = TRUE)
      mantel.tb <- mantel(beta.temporal.matrix.t, beta.temporal.matrix.b, method = "pearson", permutations = 9999, na.rm = TRUE)
      mantel.fb <- mantel(beta.temporal.matrix.f, beta.temporal.matrix.b, method = "pearson", permutations = 9999, na.rm = TRUE)
      
      # Save outputs to list
      mantel.tests[[paste0("Taxonomic_Functional_", i)]] <- mantel.tf
      mantel.tests[[paste0("Taxonomic_Biomass_", i)]] <- mantel.tb
      mantel.tests[[paste0("Functional_Biomass_", i)]] <- mantel.fb
      
      # Extract useful values from the output
      mantel.tf.output <- data.frame("PLOT" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = "Functional", "R_value" = mantel.tf$statistic, "p_value" = mantel.tf$signif)
      mantel.tb.output <- data.frame("PLOT" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = "Biomass", "R_value" = mantel.tb$statistic, "p_value" = mantel.tb$signif)
      mantel.ts.output <- data.frame("PLOT" = i, "Matrix_1" = "Taxonomic", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = NA, "p_value" = NA)
      mantel.fb.output <- data.frame("PLOT" = i, "Matrix_1" = "Functional", "Matrix_2" = "Biomass", "R_value" = mantel.fb$statistic, "p_value" = mantel.fb$signif)
      mantel.fs.output <- data.frame("PLOT" = i, "Matrix_1" = "Functional", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = NA, "p_value" = NA)
      mantel.bs.output <- data.frame("PLOT" = i, "Matrix_1" = "Biomass", "Matrix_2" = paste0("Spectral", filepath.PCA), "R_value" = NA, "p_value" = NA)
      
      # Bind dataframes to output
      mantel.output <- rbind(mantel.output, mantel.tf.output, mantel.tb.output, mantel.ts.output,
                             mantel.fb.output, mantel.fs.output, mantel.bs.output)
      
      # Remove intermediate objects
      rm(beta.temporal.matrix.t, beta.temporal.matrix.f, beta.temporal.matrix.b, mantel.tf,
         mantel.tb, mantel.fb, mantel.tf.output, mantel.tb.output, mantel.ts.output,
         mantel.fb.output, mantel.fs.output, mantel.bs.output)
      
      
    } # End of if statement for spectral NAs
    
    
    # Remove intermediate objects
    rm(beta.temporal.years, beta.temporal.cut, beta.temporal.df.t, beta.temporal.df.f,
       beta.temporal.df.b, beta.temporal.df.s)
    
    
  } # End of plot loop
  
  # Export full output for later plotting
  write.csv(mantel.output, file = paste0("outputs/output_statistics_temporal_FULL_b", buffer, "_", spectral.metric, filepath.brightness,
                                         filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Create tidy dataframe output with mean values
  mantel.output.tidy <- mantel.output %>% 
    mutate(Paired_Metrics = paste0(Matrix_1, ":", Matrix_2)) %>% 
    mutate(R2 = R_value^2) %>% 
    group_by(Paired_Metrics) %>% 
    mutate(mean_R = round(mean(R_value, na.rm = TRUE), digits = 3),
           mean_R2 = round(mean(R2, na.rm = TRUE), digits = 3)) %>% 
    ungroup() %>% 
    dplyr::select(Paired_Metrics, mean_R, mean_R2) %>% 
    distinct()
  
  # Output temporal statistics
  write.csv(mantel.output.tidy, file = paste0("outputs/output_statistics_temporal_SUMMARIZED_b", buffer, "_", spectral.metric, filepath.brightness,
                                              filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
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
           subtitle = "By Plot (All Years) - Spectral PCA",
           x = "Paired Metrics",
           y = "R Value (Mantel Test)") +
      theme_1() +
      theme(legend.position = "none"))
  
  # Export boxplot of statistics
  ggsave(mantel.boxplot.R, filename = paste0("outputs/figures/correlations_temporal_b", buffer, "_", spectral.metric, filepath.brightness,
                                             filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, "_box.png"), width = 16, height = 8)
  
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
           subtitle = "By Plot (All Years)- Spectral PCA",
           x = "Paired Metrics",
           y = "R Value (Mantel Test)") +
      coord_flip() +
      theme_1() +
      theme(legend.position = "none"))
  
  # Export raincloud plot of statistics
  ggsave(mantel.raincloud.R, filename = paste0("outputs/figures/correlations_temporal_b", buffer, "_", spectral.metric, filepath.brightness,
                                               filepath.smoothing, filepath.PCA, filepath.37, filepath.top.hits, "_rain.png"), width = 16, height = 8)
  
  
} # End of function
