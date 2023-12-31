# FUNCTION: Spectral Plots ----

# Create function to produce plots
create.spatial.temporal.plots.spectral <- function(SpectralMetric, BetaMetric, SpatialData, TemporalData, SpatialStatisticsData, TemporalStatisticsData){
  
  
  # Generate colour palettes for plots
  colours.Functional <- c("#FAC5ED", "#F56CB0", "#E64A90", "#CD1076")
  colours.Taxonomic <- c("#F29191", "#FF6961", "#ED2D2D", "#C20404")
  colours.Biomass <- c("#FADF96", "#F5D062", "#FFB30F", "#FFAE00")
  
  # colours.Functional <- c("#FFB8B8", "#FA7878", "#D62B2B", "#911313")
  # colours.Taxonomic <- c("#DAFFD4", "#8FC980", "#46A32A", "#167000")
  # colours.Biomass <- c("#D6F5FF", "#86C6DB", "#229EC7", "#1A4E91")
  
  colour.palette <- data.frame(cbind(colours.Functional, colours.Taxonomic, colours.Biomass)) %>% 
    rename(Functional = colours.Functional, Taxonomic = colours.Taxonomic, Biomass = colours.Biomass)
  
  # Create a list to put the spatial plot outputs into
  spatial.plot.list <- list()
  
  # Create dataframe for determining maximum axes values
  axes.max.values <- SpatialData %>% 
    dplyr::select(BetaMetric, SpectralMetric) %>% 
    rename(BetaMetric = paste0(BetaMetric),
           SpectralMetric = paste0(SpectralMetric))
  
  # Determine spatial axes maximums
  x.max <- max(axes.max.values$BetaMetric, na.rm = TRUE) * 1.1
  y.max <- max(axes.max.values$SpectralMetric, na.rm = TRUE) * 1.1
  
  # Create loop counter
  spatial.plot.counter <- 0
  
  
  # Run loop to produce plot for each year
  for (i in c(2017, 2018, 2019, 2020)){
    
    
    # Increase loop counter by 1
    spatial.plot.counter <- spatial.plot.counter + 1
    
    # Assign year variable
    Year <- i
    
    # Create input dataframe
    input.data.spatial <- filter(SpatialData, Year == i) %>% 
      dplyr::select(Year, BetaMetric, SpectralMetric) %>% 
      rename(BetaMetric = paste0(BetaMetric),
             SpectralMetric = paste0(SpectralMetric))
    
    # Create input dataframe with all years for axis
    input.data.spatial.all.years <- SpatialData %>% 
      dplyr::select(Year, BetaMetric, SpectralMetric) %>% 
      rename(BetaMetric = paste0(BetaMetric),
             SpectralMetric = paste0(SpectralMetric))
    
    # Produce yearly spatial plot
    spatial.plot <- ggplot() +
      geom_point(data = input.data.spatial, aes(x = BetaMetric, y = SpectralMetric),
                 shape = 21, size = 2, alpha = 0.75, colour = "#000000", fill = dplyr::select(colour.palette, BetaMetric)[spatial.plot.counter, ]) +
      stat_smooth(data = input.data.spatial, aes(x = BetaMetric, y = SpectralMetric),
                  method = lm, colour = "#000000", alpha = 1, se = FALSE, linetype = "solid") +
      scale_x_continuous(limits = c(0, x.max), expand = c(0,0)) +
      scale_y_continuous(limits = c(0.0000000001, y.max), expand = c(0,0)) +
      labs(title = ifelse(Year == 2017, "         a) Spatial \n", "\n"),
           subtitle = paste0(Year),
           caption = paste0("\n R: ",
                            round(filter(SpatialStatisticsData,
                                         Matrix_1 == BetaMetric,
                                         Matrix_2 == paste0(SpectralMetric, "_PCA"),
                                         Year == i)$R_value, digits = 3),
                            "; p: ",
                            ifelse(filter(SpatialStatisticsData,
                                          Matrix_1 == BetaMetric,
                                          Matrix_2 == paste0(SpectralMetric, "_PCA"),
                                          Year == i)$p_value < 0.001,
                                   "< 0.001",
                                   round(filter(SpatialStatisticsData,
                                                Matrix_1 == BetaMetric,
                                                Matrix_2 == paste0(SpectralMetric, "_PCA"),
                                                Year == i)$p_value, digits = 3)),
                            ifelse(round(filter(SpatialStatisticsData,
                                                Matrix_1 == BetaMetric,
                                                Matrix_2 == paste0(SpectralMetric, "_PCA"),
                                                Year == i)$p_value, digits = 3) < 0.001, "***",
                                   ifelse(round(filter(SpatialStatisticsData,
                                                       Matrix_1 == BetaMetric,
                                                       Matrix_2 == paste0(SpectralMetric, "_PCA"),
                                                       Year == i)$p_value, digits = 3) < 0.01, "**",
                                          ifelse(round(filter(SpatialStatisticsData,
                                                              Matrix_1 == BetaMetric,
                                                              Matrix_2 == paste0(SpectralMetric, "_PCA"),
                                                              Year == i)$p_value, digits = 3) < 0.05, "*", "")))),
           x = paste0("\n", BetaMetric, ifelse(BetaMetric %in% c("Taxonomic", "Functional"), " Dissimilarity", " Difference")),
           y = paste0(SpectralMetric, ifelse(SpectralMetric == "Spectral", " Dissimilarity", " Difference"), "\n")) +
      theme_ms() +
      theme(legend.position = "none",
            plot.caption = element_text(size = 19, vjust = 0, hjust = 0.5, face = "bold"))
    
    # Save plot to list
    spatial.plot.list[[paste0(BetaMetric, "_", SpectralMetric, "_", i)]] <- spatial.plot
    
    
  } # End of loop
 

  # Produce grid object of the four year plots
  spatial.panel.trio <- arrangeGrob(spatial.plot.list[[1]], spatial.plot.list[[2]],
                                    spatial.plot.list[[3]], spatial.plot.list[[4]],
                                    ncol = 2)
  
  # Create dataframe of temporal data - NEW STYLE (histogram)
  input.data.temporal <- TemporalData %>% 
    mutate(Matrix_2 = ifelse(Matrix_2 == "Spectral_PCA", "Spectral", Matrix_2)) %>% 
    filter(Matrix_2 == paste0(SpectralMetric),
           Matrix_1 == paste0(BetaMetric))
  
  # Determine the mean R value
  mean.temporal.R <- mean(input.data.temporal$R_value, na.rm = TRUE)
  
  # Plot temporal panel - NEW STYLE (histogram)
  (temporal.plot <- ggplot() +
      geom_histogram(data = input.data.temporal, aes(x = R_value), colour = "#000000",
                     alpha = 0.75, fill = dplyr::select(colour.palette, BetaMetric)[3,],
                     binwidth = 0.2) +
      geom_vline(aes(xintercept = mean.temporal.R), colour = "#000000", size = 0.5, linetype = "dashed") +
      labs(title = "         b) Temporal \n",
           subtitle = "2017 ~ 2020",
           caption = paste0("\nR: ", round(mean.temporal.R, digits = 3)),
           x = paste0("\n", BetaMetric, ifelse(BetaMetric %in% c("Taxonomic", "Functional"), " Dissimilarity", " Difference")),
           y = paste0(SpectralMetric, ifelse(SpectralMetric == "Spectral", " Dissimilarity", " Difference"), "\n")) +
      theme_ms()
  )
  
  # # Plot the temporal panel
  # (temporal.plot <- ggplot() +
  #     geom_point(data = input.data.temporal, aes(x = BetaMetric, y = SpectralMetric, fill = timeframe),
  #                shape = 21, size = 2, alpha = 0.75, colour = "#000000") +
  #     stat_smooth(data = input.data.temporal, aes(x = BetaMetric, y = SpectralMetric),
  #                 method = lm, colour = "#000000", alpha = 1, se = FALSE, linetype = "solid") +
  #     scale_fill_manual(values = dplyr::select(colour.palette, BetaMetric)[1:3,]) +
  #     scale_x_continuous(limits = c(0, max(input.data.temporal$BetaMetric)*1.1), expand = c(0,0)) +
  #     scale_y_continuous(limits = c(0.000000000001, max(input.data.temporal$SpectralMetric)*1.1), expand = c(0,0)) +
  #     labs(title = "         b) Temporal \n",
  #          subtitle = "2017 ~ 2020",
  #          caption = paste0("\n R: ", round(filter(TemporalStatisticsData,
  #                                                  Paired_Metrics == paste0(BetaMetric, ":", SpectralMetric, "_PCA"))$mean_R, digits = 3)),
  #          x = paste0("\n", BetaMetric, ifelse(BetaMetric %in% c("Taxonomic", "Functional"), " Dissimilarity", " Difference")),
  #          y = paste0(SpectralMetric, ifelse(SpectralMetric == "Spectral", " Dissimilarity", " Difference"), "\n"),
  #          fill = "Pariwise\nInterval\n(Years)") +
  #     theme_ms() +
  #     theme(legend.position = "right"))

  # Create blank plot for panel
  blank.plot <- grid::grid.rect(gp = grid::gpar(col = "white"))
  
  # Produce grid object of the single temporal plot
  temporal.panel.trio <- arrangeGrob(blank.plot, blank.plot,
                                     blank.plot, temporal.plot,
                                     blank.plot, blank.plot,
                                     ncol = 2, heights = c(0.5, 1, 0.5), widths = c(0.2, 1))
  
  # Produce combined spatial and temporal panel
  combined.panel.trio <- arrangeGrob(spatial.panel.trio, temporal.panel.trio, ncol = 2, widths = c(1.3, 0.9))
  
  combined.panel.five <- arrangeGrob(spatial.plot.list[[1]], blank.plot, spatial.plot.list[[2]], blank.plot,
                                     spatial.plot.list[[3]], blank.plot, spatial.plot.list[[4]], blank.plot, 
                                     temporal.plot, ncol = 9, widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1, 0.3, 1.3))
  
  # # Create a figure caption object
  # caption.panel <- grid::textGrob("\n Figure _ | PLACEHOLDER CAPTION",
  #                                 x = 0, y = 0.5, just = "left", gp = grid::gpar(fontsize = 18))
  # 
  # # Create combined panel with caption
  # combined.panel.with.caption <- arrangeGrob(combined.panel, caption.panel, ncol = 1, heights = c(8, 1))
  
  # Export plot
  ggsave(combined.panel.trio,
         width = 21, height = 12,
         filename = paste0("outputs/figures/manuscript/fig_", BetaMetric,
                           "_", SpectralMetric, "_trio", filepath.37, filepath.top.hits, ".png"))
  
  ggsave(combined.panel.five,
         width = 23, height = 5.4,
         filename = paste0("outputs/figures/manuscript/fig_", BetaMetric,
                           "_", SpectralMetric, "_five", filepath.37, filepath.top.hits, ".png"))

  
} # End of function


# FUNCTION: NDVI Plots ----

# Create function to produce plots
create.spatial.temporal.plots.NDVI <- function(SpectralMetric, BetaMetric, SpatialData, TemporalData, SpatialStatisticsData, TemporalStatisticsData){
  
  
  # Generate colour palettes for plots
  colours.Functional <- c("#FAC5ED", "#F56CB0", "#E64A90", "#CD1076")
  colours.Taxonomic <- c("#F29191", "#FF6961", "#ED2D2D", "#C20404")
  colours.Biomass <- c("#FADF96", "#F5D062", "#FFB30F", "#FFAE00")
  
  # colours.Functional <- c("#FFB8B8", "#FA7878", "#D62B2B", "#911313")
  # colours.Taxonomic <- c("#DAFFD4", "#8FC980", "#46A32A", "#167000")
  # colours.Biomass <- c("#D6F5FF", "#86C6DB", "#229EC7", "#1A4E91")
  
  colour.palette <- data.frame(cbind(colours.Functional, colours.Taxonomic, colours.Biomass)) %>% 
    rename(Functional = colours.Functional, Taxonomic = colours.Taxonomic, Biomass = colours.Biomass)
  
  # Create a list to put the spatial plot outputs into
  spatial.plot.list <- list()
  
  # Create dataframe for determining maximum axes values
  axes.max.values <- SpatialData %>% 
    dplyr::select(BetaMetric, SpectralMetric) %>% 
    rename(BetaMetric = paste0(BetaMetric),
           SpectralMetric = paste0(SpectralMetric))
  
  # Determine spatial axes maximums
  x.max <- max(axes.max.values$BetaMetric, na.rm = TRUE) * 1.1
  y.max <- max(axes.max.values$SpectralMetric, na.rm = TRUE) * 1.1
  
  # Create loop counter
  spatial.plot.counter <- 0
  
  
  # Run loop to produce plot for each year
  for (i in c(2017, 2018, 2019, 2020)){
    
    
    # Increase loop counter by 1
    spatial.plot.counter <- spatial.plot.counter + 1
    
    # Assign year variable
    Year <- i
    
    # Create input dataframe
    input.data.spatial <- filter(SpatialData, Year == i) %>% 
      dplyr::select(Year, BetaMetric, SpectralMetric) %>% 
      rename(BetaMetric = paste0(BetaMetric),
             SpectralMetric = paste0(SpectralMetric))
    
    # Create input dataframe with all years for axis
    input.data.spatial.all.years <- SpatialData %>% 
      dplyr::select(Year, BetaMetric, SpectralMetric) %>% 
      rename(BetaMetric = paste0(BetaMetric),
             SpectralMetric = paste0(SpectralMetric))
    
    # Produce yearly spatial plot
    spatial.plot <- ggplot() +
      geom_point(data = input.data.spatial, aes(x = BetaMetric, y = SpectralMetric),
                 shape = 21, size = 2, alpha = 0.75, colour = "#000000", fill = dplyr::select(colour.palette, BetaMetric)[spatial.plot.counter, ]) +
      stat_smooth(data = input.data.spatial, aes(x = BetaMetric, y = SpectralMetric),
                  method = lm, colour = "#000000", alpha = 1, se = FALSE, linetype = "solid") +
      scale_x_continuous(limits = c(0, x.max), expand = c(0,0)) +
      scale_y_continuous(limits = c(0.0000000001, y.max), expand = c(0,0)) +
      labs(title = ifelse(Year == 2017, "         a) Spatial \n", "\n"),
           subtitle = paste0(Year),
           caption = paste0("\n R: ",
                            round(filter(SpatialStatisticsData,
                                         Matrix_2 == BetaMetric,
                                         Matrix_1 == SpectralMetric,
                                         Year == i)$R_value, digits = 3),
                            "; p: ",
                            ifelse(filter(SpatialStatisticsData,
                                          Matrix_2 == BetaMetric,
                                          Matrix_1 == SpectralMetric,
                                          Year == i)$p_value < 0.001,
                                   "< 0.001",
                                   round(filter(SpatialStatisticsData,
                                                Matrix_2 == BetaMetric,
                                                Matrix_1 == SpectralMetric,
                                                Year == i)$p_value, digits = 3)),
                            ifelse(round(filter(SpatialStatisticsData,
                                                Matrix_2 == BetaMetric,
                                                Matrix_1 == SpectralMetric,
                                                Year == i)$p_value, digits = 3) < 0.001, "***",
                                   ifelse(round(filter(SpatialStatisticsData,
                                                       Matrix_2 == BetaMetric,
                                                       Matrix_1 == SpectralMetric,
                                                       Year == i)$p_value, digits = 3) < 0.01, "**",
                                          ifelse(round(filter(SpatialStatisticsData,
                                                              Matrix_2 == BetaMetric,
                                                              Matrix_1 == SpectralMetric,
                                                              Year == i)$p_value, digits = 3) < 0.05, "*", "")))),
           x = paste0("\n", BetaMetric, ifelse(BetaMetric %in% c("Taxonomic", "Functional"), " Dissimilarity", " Difference")),
           y = paste0(SpectralMetric, ifelse(SpectralMetric == "Spectral", " Dissimilarity", " Difference"), "\n")) +
      theme_ms() +
      theme(legend.position = "none",
            plot.caption = element_text(size = 19, vjust = 0, hjust = 0.5, face = "bold"))
    
    # Save plot to list
    spatial.plot.list[[paste0(BetaMetric, "_", SpectralMetric, "_", i)]] <- spatial.plot
    
    
  } # End of loop
  
  
  # Produce grid object of the four year plots
  spatial.panel.trio <- arrangeGrob(spatial.plot.list[[1]], spatial.plot.list[[2]],
                                    spatial.plot.list[[3]], spatial.plot.list[[4]],
                                    ncol = 2)
  
  # Create dataframe of temporal data - NEW STYLE (histogram)
  input.data.temporal <- TemporalData %>% 
    filter(Matrix_2 == paste0(BetaMetric))
  
  # Determine the mean R value
  mean.temporal.R <- mean(input.data.temporal$R_value, na.rm = TRUE)
  
  # Plot temporal panel - NEW STYLE (histogram)
  (temporal.plot <- ggplot() +
      geom_histogram(data = input.data.temporal, aes(x = R_value), colour = "#000000",
                     alpha = 0.75, fill = dplyr::select(colour.palette, BetaMetric)[3,],
                     binwidth = 0.2) +
      geom_vline(aes(xintercept = mean.temporal.R), colour = "#000000", size = 0.5, linetype = "dashed") +
      labs(title = "         b) Temporal \n",
           subtitle = "2017 ~ 2020",
           caption = paste0("\nR: ", round(mean.temporal.R, digits = 3)),
           x = paste0("\n", BetaMetric, ifelse(BetaMetric %in% c("Taxonomic", "Functional"), " Dissimilarity", " Difference")),
           y = paste0(SpectralMetric, ifelse(SpectralMetric == "Spectral", " Dissimilarity", " Difference"), "\n")) +
      theme_ms()
  )
  
  # # Plot the temporal panel - OLD STYLE (scatterplot)
  # (temporal.plot <- ggplot() +
  #     geom_point(data = input.data.temporal, aes(x = BetaMetric, y = SpectralMetric, fill = timeframe),
  #                shape = 21, size = 2, alpha = 0.75, colour = "#000000") +
  #     stat_smooth(data = input.data.temporal, aes(x = BetaMetric, y = SpectralMetric),
  #                 method = lm, colour = "#000000", alpha = 1, se = FALSE, linetype = "solid") +
  #     scale_fill_manual(values = dplyr::select(colour.palette, BetaMetric)[1:3,]) +
  #     scale_x_continuous(limits = c(0, max(input.data.temporal$BetaMetric)*1.1), expand = c(0,0)) +
  #     scale_y_continuous(limits = c(0.000000000001, max(input.data.temporal$SpectralMetric)*1.1), expand = c(0,0)) +
  #     labs(title = "         b) Temporal \n",
  #          subtitle = "2017 ~ 2020",
  #          caption = paste0("\n R: ", round(filter(TemporalStatisticsData,
  #                                                  Paired_Metrics == paste0(SpectralMetric, ":", BetaMetric))$mean_R, digits = 3)),
  #          x = paste0("\n", BetaMetric, ifelse(BetaMetric %in% c("Taxonomic", "Functional"), " Dissimilarity", " Difference")),
  #          y = paste0(SpectralMetric, ifelse(SpectralMetric == "Spectral", " Dissimilarity", " Difference"), "\n"),
  #          fill = "Pariwise\nInterval\n(Years)") +
  #     theme_ms() +
  #     theme(legend.position = "right"))
  
  # Create blank plot for panel
  blank.plot <- grid::grid.rect(gp = grid::gpar(col = "white"))
  
  # Produce grid object of the single temporal plot
  temporal.panel.trio <- arrangeGrob(blank.plot, blank.plot,
                                     blank.plot, temporal.plot,
                                     blank.plot, blank.plot,
                                     ncol = 2, heights = c(0.5, 1, 0.5), widths = c(0.2, 1))
  
  # Produce combined spatial and temporal panel
  combined.panel.trio <- arrangeGrob(spatial.panel.trio, temporal.panel.trio, ncol = 2, widths = c(1.3, 0.9))
  
  combined.panel.five <- arrangeGrob(spatial.plot.list[[1]], blank.plot, spatial.plot.list[[2]], blank.plot,
                                     spatial.plot.list[[3]], blank.plot, spatial.plot.list[[4]], blank.plot, 
                                     temporal.plot, ncol = 9, widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1, 0.3, 1.3))
  
  # # Create a figure caption object
  # caption.panel <- grid::textGrob("\n Figure _ | PLACEHOLDER CAPTION",
  #                                 x = 0, y = 0.5, just = "left", gp = grid::gpar(fontsize = 18))
  # 
  # # Create combined panel with caption
  # combined.panel.with.caption <- arrangeGrob(combined.panel, caption.panel, ncol = 1, heights = c(8, 1))
  
  # Export plot
  ggsave(combined.panel.trio,
         width = 21, height = 12,
         filename = paste0("outputs/figures/manuscript/fig_", BetaMetric,
                           "_", SpectralMetric, "_trio.png"))
  
  ggsave(combined.panel.five,
         width = 23, height = 5.4,
         filename = paste0("outputs/figures/manuscript/fig_", BetaMetric,
                           "_", SpectralMetric, "_five.png"))
  
  
} # End of function
