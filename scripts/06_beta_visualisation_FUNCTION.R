# 06 - Functions for comparing Beta Diversity - Taxonomic, Functional & Spectral
# Joseph Everest
# February 2023, adapted March 2023


# FUNCTION: SPATIAL VISUALISATIONS ----

beta.visualisations.spatial <- function(beta.spatial, spectral.metric){
  
  
  # Generate input dataframe depending on the spectral metric selected
  if (spectral.metric == "Euclidean"){

    
    # Generate input dataframe
    beta.spatial.input <- beta.spatial %>% 
      dplyr::select(Year, Taxonomic_Dis, Functional_Dis, Biomass_Dis, Spectral_PCA_Euclidean_Dis) %>% 
      rename(Spectral_Dis = Spectral_PCA_Euclidean_Dis)
    
    
  } # End of Euclidean if statement
  
  
  if (spectral.metric == "Manhattan"){
    
    
    # Generate input dataframe
    beta.spatial.input <- beta.spatial %>% 
      dplyr::select(Year, Taxonomic_Dis, Functional_Dis, Biomass_Dis, Spectral_PCA_Manhattan_Dis) %>% 
      rename(Spectral_Dis = Spectral_PCA_Manhattan_Dis)
    
    
  } # End of Manhattan if statement
  
  
  if (spectral.metric == "SAM"){
    
    
    # Generate input dataframe
    beta.spatial.input <- beta.spatial %>% 
      dplyr::select(Year, Taxonomic_Dis, Functional_Dis, Biomass_Dis, Spectral_PCA_SAM_Dis) %>% 
      rename(Spectral_Dis = Spectral_PCA_SAM_Dis)
    
    
  } # End of SAM if statement
  
  # Create string for plot subtitle
  plots.included <- "PLOTS: NDVI Threshold (O.2)"
  
  # Visualize relationship between taxonomic and functional beta diversity
  (plot.beta.spatial.t.f <- ggplot(data = beta.spatial.input, aes(x = Taxonomic_Dis, y = Functional_Dis, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "magma", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.spatial.input$Taxonomic_Dis), max(beta.spatial.input$Taxonomic_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.spatial.input$Functional_Dis), max(beta.spatial.input$Functional_Dis)), expand = c(0,0)) +
      labs(title = "Taxonomic vs Functional Dissimilarity [SPATIAL]",
           subtitle = plots.included,
           x = "Taxonomic (Bray-Curtis) Dissimilarity",
           y = "Functional (FDis) Dissimilarity",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between taxonomic and biomass beta diversity
  (plot.beta.spatial.t.b <- ggplot(data = beta.spatial.input, aes(x = Taxonomic_Dis, y = Biomass_Dis, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "viridis", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.spatial.input$Taxonomic_Dis), max(beta.spatial.input$Taxonomic_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.spatial.input$Biomass_Dis), max(beta.spatial.input$Biomass_Dis)), expand = c(0,0)) +
      labs(title = "Taxonomic vs Biomass Dissimilarity [SPATIAL]",
           subtitle = plots.included,
           x = "Taxonomic (Bray-Curtis) Dissimilarity",
           y = "Biomass (Euclidean) Dissimilarity",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between taxonomic and spectral beta diversity
  (plot.beta.spatial.t.s <- ggplot(data = beta.spatial.input, aes(x = Taxonomic_Dis, y = Spectral_Dis, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.spatial.input$Taxonomic_Dis), max(beta.spatial.input$Taxonomic_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.spatial.input$Spectral_Dis), max(beta.spatial.input$Spectral_Dis)), expand = c(0,0)) +
      labs(title = "Taxonomic vs Spectral PCA Dissimilarity [SPATIAL]",
           subtitle = plots.included,
           x = "Taxonomic (Bray-Curtis) Dissimilarity",
           y = paste0("Spectral PCA (", spectral.metric, ") Dissimilarity"),
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between functional and biomass beta diversity
  (plot.beta.spatial.f.b <- ggplot(data = beta.spatial.input, aes(x = Functional_Dis, y = Biomass_Dis, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "viridis", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.spatial.input$Functional_Dis), max(beta.spatial.input$Functional_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.spatial.input$Biomass_Dis), max(beta.spatial.input$Biomass_Dis)), expand = c(0,0)) +
      labs(title = "Functional vs Biomass Dissimilarity [SPATIAL]",
           subtitle = plots.included,
           x = "Functional (FDis) Dissimilarity",
           y = "Biomass (Euclidean) Dissimilarity",
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between functional and spectral beta diversity
  (plot.beta.spatial.f.s <- ggplot(data = beta.spatial.input, aes(x = Functional_Dis, y = Spectral_Dis, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.spatial.input$Functional_Dis), max(beta.spatial.input$Functional_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.spatial.input$Spectral_Dis), max(beta.spatial.input$Spectral_Dis)), expand = c(0,0)) +
      labs(title = "Functional vs Spectral PCA Dissimilarity [SPATIAL]",
           subtitle = plots.included,
           x = "Functional (FDis) Dissimilarity",
           y = paste0("Spectral PCA (", spectral.metric, ") Dissimilarity"),
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between biomass and spectral beta diversity
  (plot.beta.spatial.b.s <- ggplot(data = beta.spatial.input, aes(x = Biomass_Dis, y = Spectral_Dis, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.spatial.input$Biomass_Dis), max(beta.spatial.input$Biomass_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.spatial.input$Spectral_Dis), max(beta.spatial.input$Spectral_Dis)), expand = c(0,0)) +
      labs(title = "Biomass vs Spectral PCA Dissimilarity [SPATIAL]",
           subtitle = plots.included,
           x = "Biomass (Euclidean) Dissimilarity",
           y = paste0("Spectral PCA (", spectral.metric, ") Dissimilarity"),
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Create blank plot for panel
  plot.spatial.blank <- grid::grid.rect(gp = grid::gpar(col = "white"))
  
  # Generate a panel of the outputs
  panel.spatial <- grid.arrange(plot.beta.spatial.t.f, plot.beta.spatial.t.b, plot.beta.spatial.t.s,
                                plot.spatial.blank, plot.beta.spatial.f.b, plot.beta.spatial.f.s,
                                plot.spatial.blank, plot.spatial.blank, plot.beta.spatial.b.s,
                                ncol = 3)
  
  # Export plot
  ggsave(panel.spatial, file = paste0("outputs/figures/beta_spatial_", spectral.metric, filepath.37, filepath.top.hits, ".png"), width = 21, height = 15)
  
  
} # End of function


# FUNCTION: TEMPORAL VISUALISATIONS ----

beta.visualisations.temporal <- function(beta.temporal, spectral.metric){
  
  
  # Generate input dataframe depending on the spectral metric selected
  if (spectral.metric == "Euclidean"){
    
    
    # Generate input dataframe
    beta.temporal.input <- beta.temporal %>% 
      dplyr::select(Plot, Years, Taxonomic_Dis, Functional_Dis, Biomass_Dis, Spectral_PCA_Euclidean_Dis) %>% 
      rename(Spectral_Dis = Spectral_PCA_Euclidean_Dis)
    
    
  } # End of Euclidean if statement
  
  
  if (spectral.metric == "Manhattan"){
    
    
    # Generate input dataframe
    beta.temporal.input <- beta.temporal %>% 
      dplyr::select(Plot, Years, Taxonomic_Dis, Functional_Dis, Biomass_Dis, Spectral_PCA_Manhattan_Dis) %>% 
      rename(Spectral_Dis = Spectral_PCA_Manhattan_Dis)
    
    
  } # End of Manhattan if statement
  
  
  if (spectral.metric == "SAM"){
    
    
    # Generate input dataframe
    beta.temporal.input <- beta.temporal %>% 
      dplyr::select(Plot, Years, Taxonomic_Dis, Functional_Dis, Biomass_Dis, Spectral_PCA_SAM_Dis) %>% 
      rename(Spectral_Dis = Spectral_PCA_SAM_Dis)
    
    
  } # End of SAM if statement
 
  # Create string for plot subtitle
  plots.included <- "PLOTS: NDVI Threshold (O.2)"
  
  # Visualize relationship between taxonomic and functional beta diversity
  (plot.beta.temporal.input.t.f <- ggplot(data = beta.temporal.input, aes(x = Taxonomic_Dis, y = Functional_Dis, fill = Years)) +
      geom_point(shape = 21, size = 4, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "magma", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.temporal.input$Taxonomic_Dis), max(beta.temporal.input$Taxonomic_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.temporal.input$Functional_Dis), max(beta.temporal.input$Functional_Dis)), expand = c(0,0)) +
      labs(title = "Taxonomic vs Functional Dissimilarity [TEMPORAL]",
           subtitle = plots.included,
           x = "Taxonomic (Bray-Curtis) Dissimilarity",
           y = "Functional (FDis) Dissimilarity") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between taxonomic and biomass beta diversity
  (plot.beta.temporal.input.t.b <- ggplot(data = beta.temporal.input, aes(x = Taxonomic_Dis, y = Biomass_Dis, fill = Years)) +
      geom_point(shape = 21, size = 4, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "viridis", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.temporal.input$Taxonomic_Dis), max(beta.temporal.input$Taxonomic_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.temporal.input$Biomass_Dis), max(beta.temporal.input$Biomass_Dis)), expand = c(0,0)) +
      labs(title = "Taxonomic vs Biomass Dissimilarity [TEMPORAL]",
           subtitle = plots.included,
           x = "Taxonomic (Bray-Curtis) Dissimilarity",
           y = "Biomass (Euclidean) Dissimilarity") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between taxonomic and spectral beta diversity
  (plot.beta.temporal.input.t.s <- ggplot(data = beta.temporal.input, aes(x = Taxonomic_Dis, y = Spectral_Dis, fill = Years)) +
      geom_point(shape = 21, size = 4, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.temporal.input$Taxonomic_Dis), max(beta.temporal.input$Taxonomic_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.temporal.input$Spectral_Dis), max(beta.temporal.input$Spectral_Dis)), expand = c(0,0)) +
      labs(title = "Taxonomic vs PCA Spectral Dissimilarity [TEMPORAL]",
           subtitle = plots.included,
           x = "Taxonomic (Bray-Curtis) Dissimilarity",
           y = paste0("Spectral PCA (", spectral.metric, ") Dissimilarity")) +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between functional and biomass beta diversity
  (plot.beta.temporal.input.f.b <- ggplot(data = beta.temporal.input, aes(x = Functional_Dis, y = Biomass_Dis, fill = Years)) +
      geom_point(shape = 21, size = 4, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "viridis", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.temporal.input$Functional_Dis), max(beta.temporal.input$Functional_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.temporal.input$Biomass_Dis), max(beta.temporal.input$Biomass_Dis)), expand = c(0,0)) +
      labs(title = "Functional vs Biomass Dissimilarity [TEMPORAL]",
           subtitle = plots.included,
           x = "Functional (FDis) Dissimilarity",
           y = "Biomass (Euclidean) Dissimilarity") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between functional and spectral beta diversity
  (plot.beta.temporal.input.f.s <- ggplot(data = beta.temporal.input, aes(x = Functional_Dis, y = Spectral_Dis, fill = Years)) +
      geom_point(shape = 21, size = 4, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.temporal.input$Functional_Dis), max(beta.temporal.input$Functional_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.temporal.input$Spectral_Dis), max(beta.temporal.input$Spectral_Dis)), expand = c(0,0)) +
      labs(title = "Functional vs PCA Spectral Dissimilarity [TEMPORAL]",
           subtitle = plots.included,
           x = "Functional (FDis) Dissimilarity",
           y = paste0("Spectral PCA (", spectral.metric, ") Dissimilarity")) +
      theme_1() +
      theme(legend.position = "right"))
  
  # Visualize relationship between biomass and spectral beta diversity
  (plot.beta.temporal.input.b.s <- ggplot(data = beta.temporal.input, aes(x = Biomass_Dis, y = Spectral_Dis, fill = Years)) +
      geom_point(shape = 21, size = 4, alpha = 1, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(beta.temporal.input$Biomass_Dis), max(beta.temporal.input$Biomass_Dis)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(beta.temporal.input$Spectral_Dis), max(beta.temporal.input$Spectral_Dis)), expand = c(0,0)) +
      labs(title = "Biomass vs PCA Spectral Dissimilarity [TEMPORAL]",
           subtitle = plots.included,
           x = "Biomass (Euclidean) Dissimilarity",
           y = paste0("Spectral PCA (", spectral.metric, ") Dissimilarity")) +
      theme_1() +
      theme(legend.position = "right"))
  
  # Create blank plot for panel
  plot.temporal.blank <- grid::grid.rect(gp = grid::gpar(col = "white"))
  
  # Generate a panel of the outputs
  panel.temporal <- grid.arrange(plot.beta.temporal.input.t.f, plot.beta.temporal.input.t.b, plot.beta.temporal.input.t.s,
                                 plot.temporal.blank, plot.beta.temporal.input.f.b, plot.beta.temporal.input.f.s,
                                 plot.temporal.blank, plot.temporal.blank, plot.beta.temporal.input.b.s,
                                 ncol = 3)
  
  # Export plot
  ggsave(panel.temporal, file = paste0("outputs/figures/beta_temporal_", spectral.metric, filepath.37, filepath.top.hits, ".png"), width = 21, height = 15)
  
   
} # End of function
