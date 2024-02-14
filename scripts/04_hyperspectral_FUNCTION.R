# 04a - Functions for working with NEON hyperspectral data
# Joseph Everest & NEON
# February 2023, modified April 2023, May 2023, September 2023


# LOAD PACKAGES & THEMES ----

# Load packages
library(tidyverse)
library(viridis)
library(gridExtra)
library(rgdal)
library(raster)
# install.packages("BiocManager")
# BiocManager::install("rhdf5")
library(rhdf5)
library(neonhs)

# Load themes
source("hyperspectral/scripts/EX1_ggplot_themes.R")


# FUNCTION: IMPORT HYPERSPECTRAL TILE ----

  # Imports full NEON hyperspectral tile, function created in NEON tutorial (see EX2)

import.hyperspec.tile <- function(file, band, noDataValue, extent, CRS){
  
  # First, read in the raster
  tile <- h5read(file, "/NIWO/Reflectance/Reflectance_Data", index = list(band, NULL, NULL))
  
  # Convert from array to matrix
  tile <- (tile[1, , ])
  
  # Transpose data to fix flipped row and column order 
  tile <- t(tile)

  # Assign data ignore values to NA
  tile[tile == noDataValue] <- NA
  
  # Turn the object into a raster
  tile.r <- raster(tile, crs = CRS)
  
  # Assign the extent to the raster
  extent(tile.r) <- extent
  
  # Return the raster object
  return(tile.r)
  
}


# FUNCTION: TILE INFORMATION ----

  # Extracts information from metadata for each tile (filepath) listed in vector

extract.tile.metadata <- function(tile.list){

  
  # Create empty dataframe for outputting tile information
  tile.information <- data.frame()
  
  
  # Run a loop to get important information on each raster
  for (i in tile.list){
    
    # Generate simplified tile name
    tile.name <- i %>% 
      str_remove_all("data/NEON/tiles/NEON_D13_NIWO_DP3_") %>% 
      str_remove_all("_reflectance.h5") %>% 
      str_remove_all("_449000") %>% 
      str_remove_all("000")
    
    # Extract the tile reflectance metadata
    tile.metadata <- h5readAttributes(i, "/NIWO/Reflectance/Reflectance_Data")
    
    # Read the dimensions of the hyperspectral tile
    tile.nrows <- tile.metadata$Dimensions[1]
    tile.ncols <- tile.metadata$Dimensions[2]
    tile.nbands <- tile.metadata$Dimensions[3]
    
    # Define the tile no data value
    tile.NoDataValue <- as.numeric(tile.metadata$Data_Ignore_Value)
    
    # Define the tile CRS
    tile.EPSG <- h5read(i, "/NIWO/Reflectance/Metadata/Coordinate_System/EPSG Code")
    
    # Determine extent values of tile
    tile.xmin <- tile.metadata$Spatial_Extent_meters[1]
    tile.xmax <- tile.metadata$Spatial_Extent_meters[2]
    tile.ymin <- tile.metadata$Spatial_Extent_meters[3]
    tile.ymax <- tile.metadata$Spatial_Extent_meters[4]
    
    # Extract scale factor from reflectance attributes
    tile.scalefactor <- tile.metadata$Scale_Factor
    
    # Create a dataframe of these values
    tile.output <- data.frame("Tile_Name" = tile.name,
                              "Rows" = tile.nrows,
                              "Cols" = tile.ncols,
                              "Bands" = tile.nbands,
                              "NoDataValue" = tile.NoDataValue,
                              "EPSG_Code" = tile.EPSG,
                              "Scale_Factor" = tile.scalefactor,
                              "xmin" = tile.xmin,
                              "xmax" = tile.xmax,
                              "ymin" = tile.ymin,
                              "ymax" = tile.ymax) %>% 
      separate(Tile_Name, c("Year", "Tile"), sep = "_")
    
    # Append output to dataframe
    tile.information <- rbind(tile.information, tile.output)
    
    # Remove intermediate objects
    rm(i, tile.name, tile.metadata, tile.nrows, tile.ncols, tile.nbands, tile.NoDataValue,
       tile.EPSG, tile.xmin, tile.xmax, tile.ymin, tile.ymax, tile.output)
    
    
  } # End of loop
  
  
  # Write csv of output
  write.csv(tile.information, file = paste0("outputs/output_NIWOT_tile_metadata", filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Remove intermediate objects
  rm(tile.information)
  
  
} # End of function



# FUNCTION: PIXEL EXTRACTION ----

extract.pixel.spectra <- function(tile.list, year.list, buffer){
  
  
  # Load in tile information
  tile.information <- read.csv(paste0("outputs/output_NIWOT_tile_metadata", filepath.top.hits, ".csv"))
  
  # Import saddle locations
  saddle.locations <- read.csv(paste0("outputs/output_saddle_locations", filepath.top.hits, ".csv"))
  
  # Create key of saddle locations plots to ID
  saddle.key <- saddle.locations %>% 
    mutate(id = as.numeric(row.names(.))) %>% 
    relocate(id, .before = )
  
  # Get just two required columns for saddle extraction
  saddle.locations.spatial <- saddle.locations %>% 
    dplyr::select(PLOT, EASTING, NORTHING) %>% 
    sf::st_as_sf(., coords = c("EASTING", "NORTHING"), crs = 26913, remove = F) %>% 
    sf::st_transform(., crs = as.numeric(unique(tile.information$EPSG_Code)))
  
  # Import the tile wavelengths from metadata
  saddle.bands.remove.1 <- h5read(tile.list[1], "/NIWO/Reflectance/Metadata/Spectral_Data/Wavelength")
  
  # Create dataframe with band numbers
  saddle.bands.remove.2 <- data.frame(c(1:426), saddle.bands.remove.1)
  
  # Assign names to columns
  names(saddle.bands.remove.2)[1] <- "Band"
  names(saddle.bands.remove.2)[2] <- "Wavelength"
  
  # Determine which bands to remove depending on wavelengths
  saddle.bands.remove.3 <- saddle.bands.remove.2 %>% 
    mutate(REMOVE = ifelse(Wavelength <= 400, TRUE, FALSE),
           REMOVE = ifelse(Wavelength >= 2400, TRUE, REMOVE),
           REMOVE = ifelse(Wavelength >= 1340 & Wavelength <= 1445, TRUE, REMOVE),
           REMOVE = ifelse(Wavelength >= 1790 & Wavelength <= 1955, TRUE, REMOVE)) %>% 
    filter(REMOVE == TRUE) %>% 
    dplyr::select(Band) %>% 
    unique()
  
  # Create vector of those bands
  saddle.bands.remove.4 <- unique(saddle.bands.remove.3$Band)
  
  # Create output dataframe for outputting spectral outputs
  spectra.output <- data.frame()
  
  
  # Run loop to extract spectra at each plot location for each tile
  for (i in tile.list){
    
    
    # Extract tile name from the input filepath
    tile.ID <- i %>% 
      str_remove_all("data/NEON/tiles/NEON_D13_NIWO_DP3_") %>% 
      str_remove_all("_reflectance.h5") %>% 
      str_remove_all("_449000") %>% 
      str_remove_all("000")
    
    # Filter tile information to that tile
    tile.information.cut <- tile.information %>% 
      mutate(ID = paste0(Year, "_", Tile)) %>% 
      relocate(ID, .after = Tile) %>% 
      filter(ID == tile.ID)
    
    # Determine input parameters for importing hyperspectral tile
    tile.bands <- c(1:426)
    tile.NoDataValue <- unique(tile.information.cut$NoDataValue)
    tile.filename <- i
    tile.crs <- crs(paste0("+init=epsg:", unique(tile.information.cut$EPSG_Code)))
    tile.extent <- extent(unique(tile.information.cut$xmin),
                          unique(tile.information.cut$xmax),
                          unique(tile.information.cut$ymin),
                          unique(tile.information.cut$ymax))
    
    # Import the tile
    tile.raster.hyperspec <- lapply(tile.bands,
                                    FUN = import.hyperspec.tile,
                                    file = tile.filename,
                                    noDataValue = tile.NoDataValue,
                                    extent = tile.extent,
                                    CRS =  tile.crs)
    
    # Now create a raster stack from those 426 rasters
    tile.stack.hyperspec <- stack(tile.raster.hyperspec)
    
    # Create a list of band names
    tile.band.names <- paste("Band_", tile.bands, sep = "")
    
    # Set the rasterStack's names to the list of bandNames created above
    names(tile.stack.hyperspec) <- tile.band.names
    
    # Create dataframe to append each extracted band to (ID == plot, see saddle.key)
    saddle.spectra.1 <- as.data.frame(1:79) %>% 
      rename(ID = "1:79")
    
    
    # Take single raster for each band to extract reflectance from with buffer
    for (j in 1:426){
      
      
      # Create individual raster to extract reflectance for
      tile.specific.band <- tile.stack.hyperspec[[j]]
      
      # Extract spectra at the saddle point locations, with a 3m buffer
      saddle.spectra.indiv.band <- as.data.frame(raster::extract(tile.specific.band, saddle.locations.spatial, method = "simple", buffer = buffer, fun = "mean"))
      
      # Rename column to reflectance
      names(saddle.spectra.indiv.band)[1] <- paste0("Band_", j)
      
      # Modify dataframe and join to overall tile output
      saddle.spectra.1 <- saddle.spectra.indiv.band %>% 
        mutate(ID = as.numeric(row.names(.))) %>% 
        left_join(saddle.spectra.1, ., by = c("ID" = "ID"))
      
      
    } # End of loop running through each band within a single tile
    
    
    # Import the tile wavelengths from metadata
    tile.wavelengths <- h5read(i, "/NIWO/Reflectance/Metadata/Spectral_Data/Wavelength")
    
    # Create a wavelengths key
    tile.wavelengths.key <- data.frame(tile.bands, tile.wavelengths) %>% 
      rename(Band = tile.bands, Wavelength = tile.wavelengths)
    
    # Remove plots that didn't fall within that tile (all bands == NA)
    saddle.spectra.2 <- saddle.spectra.1 %>% 
      na.omit(.)
    
    # Turn saddle spectra output into usable dataframe AND unscale reflectance
    saddle.spectra.3 <- saddle.spectra.2 %>% 
      pivot_longer(cols = 2:ncol(.), names_to = "Band", values_to = "Reflectance") %>% 
      mutate(scale_Reflectance = Reflectance / as.integer(tile.information.cut$Scale_Factor),
             Band = as.integer(str_remove(Band, pattern = "Band_"))) %>% 
      left_join(., tile.wavelengths.key, by = c("Band" = "Band")) %>% 
      left_join(., saddle.key, by = c("ID" = "id")) %>% 
      relocate(PLOT, ID, EASTING, NORTHING, Band, Wavelength, .before = ) %>% 
      rename(id = ID)
    
    # Input plot information for joining back to main dataframe
    saddle.spectra.4 <- saddle.spectra.3 %>% 
      mutate(Tile = tile.ID) %>% 
      separate(Tile, into = c("Year", "Tile"), sep = "_") %>% 
      mutate(Tile = paste0(Year, "_449000_", Tile, "000")) %>% 
      relocate(Year, Tile, .after = id)
    
    # Run loop to determine whether brightness normalizing or not
    if (brightness == "Yes"){
      
      # Brightness normalise the spectra for each plot
          # Divide each band by the sqrt of the sum of squared reflectances
      saddle.spectra.5 <- saddle.spectra.4 %>% 
        group_by(PLOT, Year) %>% 
        mutate(bn_key = sqrt(sum(scale_Reflectance^2))) %>% 
        ungroup() %>% 
        mutate(bn_Reflectance = scale_Reflectance / bn_key) %>% 
        dplyr::select(-bn_key)
      
    } else {
      
      # Don't brightness normalize
      saddle.spectra.5 <- saddle.spectra.4 %>% 
        mutate(bn_Reflectance = scale_Reflectance)
      
    }

    # Run loop to determine whether smoothing or not
    if (smoothing == "Yes"){
      
      # Create dataframe of spectra for smoothing
      saddle.smooth.1 <- saddle.spectra.5 %>% 
        dplyr::select(id, Band, bn_Reflectance) %>% 
        pivot_wider(names_from = "Band", values_from = "bn_Reflectance") %>% 
        column_to_rownames(var = "id")
      
      # Create vector of remaining wavelengths and bands after filtering
      saddle.wavelengths <- sort(unique(saddle.spectra.5$Wavelength))
      
      # Turn into class object 'speclib' (DON'T LOAD "hsdar" PACKAGE (confuses tidyverse)
      saddle.smooth.2 <- hsdar::speclib(as.matrix(saddle.smooth.1), saddle.wavelengths)
      
      # Apply smoothing function
      saddle.smooth.3 <- hsdar::noiseFiltering(saddle.smooth.2, method = 'sgolay', n = 7)
      
      # Turn smoothed speclib object back into matrix
      saddle.smooth.4 <- as.data.frame(hsdar::spectra(saddle.smooth.3))
      
      # Create vector of IDs for joining back to main dataframe
      saddle.ids <- rownames(saddle.smooth.1)
      
      # Create vector of remaining wavelengths and bands after filtering
      saddle.bands <- unique(saddle.spectra.5$Band)
      
      # Assign row and column names to dataframe for rejoining
      rownames(saddle.smooth.4) <- saddle.ids
      colnames(saddle.smooth.4) <- saddle.bands
      
      # Turn back into long format for rejoining
      saddle.smooth.5 <- saddle.smooth.4 %>% 
        rownames_to_column(var = "id") %>% 
        mutate(id = as.integer(id)) %>% 
        pivot_longer(cols = 2:ncol(.), names_to = "Band", values_to = "smooth_Reflectance") %>% 
        mutate(Band = as.integer(Band))
      
      # Join smoothed reflectance data back to the main spectral output for this tile
      saddle.spectra.6 <- left_join(saddle.spectra.5, saddle.smooth.5, by = c("id" = "id", "Band" = "Band")) %>% 
        dplyr::select(-id) %>% 
        rename(Plot = PLOT)
      
    } else {
      
      # Don't smooth
      saddle.spectra.6 <- saddle.spectra.5 %>% 
        mutate(smooth_Reflectance = bn_Reflectance)
      
    }
    
    # Filter water absorption and extreme bands from the output
    saddle.spectra.7 <- saddle.spectra.6 %>% 
      mutate(REMOVE = ifelse(Band %in% saddle.bands.remove.4, TRUE, FALSE)) %>% 
      filter(REMOVE != TRUE) %>% 
      dplyr::select(-REMOVE)

    # Append output to main output dataframe
    spectra.output <- rbind(spectra.output, saddle.spectra.7)
    
    # Remove intermediate objects
    # rm(tile.ID, tile.information.cut, tile.bands, tile.NoDataValue, tile.filename, tile.crs,
    #    tile.extent, tile.raster.hyperspec, tile.stack.hyperspec, tile.band.names, tile.specific.band,
    #    saddle.spectra.indiv.band, tile.wavelengths, tile.wavelengths.key, saddle.smooth.1,
    #    saddle.wavelengths, saddle.smooth.2, saddle.smooth.3, saddle.smooth.4, saddle.ids,
    #    saddle.bands, saddle.smooth.5, saddle.spectra.1, saddle.spectra.2, saddle.spectra.3,
    #    saddle.spectra.4, saddle.spectra.5, saddle.spectra.6, saddle.spectra.7)
    
    
  } # End of tile loop
  
  
  # Tidy the spectra output
  spectra.output.tidy <- spectra.output %>%
    arrange(Year, PLOT, Wavelength)
  
  # Calculate the mean NIR for each plot year combo
  spectra.NIR <- spectra.output.tidy %>%
    dplyr::select(Year, PLOT, Wavelength, scale_Reflectance) %>% 
    filter(Wavelength >= 752 & Wavelength <= 1048) %>% 
    group_by(Year, PLOT) %>% 
    summarise(mean_NIR = mean(scale_Reflectance)) %>%
    ungroup()
  
  # Join this information by plot and year to the main output dataframe
  spectra.output.NIR <- left_join(spectra.output.tidy, spectra.NIR, by = c("Year" = "Year", "PLOT" = "PLOT"))
  
  # Generate broad-band NDVI using Sentinel-2 bands
  spectra.NDVI.broad <- spectra.output.NIR %>% 
    dplyr::select(PLOT, Year, Wavelength, scale_Reflectance) %>% 
    mutate(Sentinel_Band = NA, # Label wavelengths in Sentinel-2 bands
           Sentinel_Band = ifelse(Wavelength > 633 & Wavelength < 695, "Band_4", Sentinel_Band),
           Sentinel_Band = ifelse(Wavelength > 726 & Wavelength < 938, "Band_8", Sentinel_Band)) %>% 
    filter(Sentinel_Band %in% c("Band_4", "Band_8")) %>% # Retain only Sentintel-2 bands
    dplyr::select(-Wavelength) %>% 
    group_by(PLOT, Year, Sentinel_Band) %>% 
    summarise(mean_Reflectance = mean(scale_Reflectance)) %>% # Calculate mean reflectance for each Sentinel-2 band
    ungroup() %>% 
    pivot_wider(names_from = "Sentinel_Band", values_from = "mean_Reflectance") %>% 
    mutate(NDVI_broad = (Band_8 - Band_4) / (Band_8 + Band_4)) %>% 
    dplyr::select(-c(Band_8, Band_4))
  
  # Generate narrow-band NDVI using NEON specified bands
  spectra.NDVI.narrow <- spectra.output.NIR %>% 
    dplyr::select(PLOT, Year, Band, scale_Reflectance) %>% 
    filter(Band %in% c(58, 90)) %>% 
    mutate(Band = paste0("Band_", Band)) %>% 
    pivot_wider(names_from = "Band", values_from = "scale_Reflectance") %>% 
    mutate(NDVI_narrow = (Band_90 - Band_58) / (Band_90 + Band_58)) %>% 
    dplyr::select(-c(Band_90, Band_58))
  
  # Join two NDVI outputs together in one dataframe
  spectra.NDVI <- left_join(spectra.NDVI.broad, spectra.NDVI.narrow, by = c("Year" = "Year", "PLOT" = "PLOT"))

  # Join this information by plot and year to the main output dataframe
  spectra.output.NDVI <- left_join(spectra.output.NIR, spectra.NDVI, by = c("Year" = "Year", "PLOT" = "PLOT"))
  
  # Save the overall spectra output
  write.csv(spectra.output.NDVI, file = paste0("outputs/output_saddle_spectra_b", buffer, 
                                               filepath.brightness, filepath.smoothing,
                                               filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Remove intermediate objects
  rm(spectra.output, spectra.NIR, spectra.output.NIR, spectra.NDVI, spectra.output.NDVI)
  
  
  # Run loop to generate signatures for each plot in each year
  for (k in year.list){
    
    
    # Plot the spectral signature for each plot using ggplot() and facet_wrap()
    (pixel.signatures <- ggplot(data = filter(spectra.output.tidy, Year == k)) +
       geom_point(aes(x = Wavelength, y = smooth_Reflectance), colour = "#850101", size = 0.2) +
       labs(x = "Wavelength (nm)",
            y = "Reflectance (BN & Smooth)",
            title = "Saddle Spectra by Plot",
            subtitle = paste0("Year: ", k)) +
       facet_wrap(~ PLOT)) + 
      theme_1()
    
    # Save panel of spectra
    ggsave(pixel.signatures, filename = paste0("outputs/figures/saddle_spectra_by_plot_b", buffer,
                                               "_", k, filepath.brightness, filepath.smoothing,
                                               filepath.top.hits, ".png"), width = 12, height = 9)
    
    
  } # End of year loop
  
  
  # Remove intermediate objects
  rm(spectra.output.tidy, tile.information, saddle.spdf, saddle.key, saddle.bands.remove.1,
     saddle.bands.remove.2, saddle.bands.remove.3, saddle.bands.remove.4)
  
  
} # End of function


# FUNCTION: BUFFER COMPARISON PLOTS ----

extract.comparison.plots <- function(buffers){
  
  
  # Run if statement to ensure importing spectra of specified buffers
  if (length(buffers) == 3){
    
    
    # Import dataframes for each of the specified buffers
    buffer.input.1 <- read.csv(paste0("outputs/output_saddle_spectra_b", buffers[1],
                                      filepath.brightness, filepath.smoothing, filepath.PCA, filepath.top.hits, ".csv"))
    buffer.input.2 <- read.csv(paste0("outputs/output_saddle_spectra_b", buffers[2], 
                                      filepath.brightness, filepath.smoothing, filepath.PCA, filepath.top.hits, ".csv"))
    buffer.input.3 <- read.csv(paste0("outputs/output_saddle_spectra_b", buffers[3], 
                                      filepath.brightness, filepath.smoothing, filepath.PCA, filepath.top.hits, ".csv"))
    
    # Determine what the buffers are
    buffer.1 <- paste0(buffers[1])
    buffer.2 <- paste0(buffers[2])
    buffer.3 <- paste0(buffers[3])
    
    # Modify the reflectance for if shaded (e.g. NIR < 0.2)#
    buffer.input.1 <- mutate(buffer.input.1, smooth_Reflectance = ifelse(mean_NIR < 0.2, NA, smooth_Reflectance))
    buffer.input.2 <- mutate(buffer.input.2, smooth_Reflectance = ifelse(mean_NIR < 0.2, NA, smooth_Reflectance))
    buffer.input.3 <- mutate(buffer.input.3, smooth_Reflectance = ifelse(mean_NIR < 0.2, NA, smooth_Reflectance))
    
    # Rename the smooth_Reflectance column for the dataframe with the buffer distance
    colnames(buffer.input.1)[colnames(buffer.input.1) == "smooth_Reflectance"] <- "Buffer_1"
    colnames(buffer.input.2)[colnames(buffer.input.2) == "smooth_Reflectance"] <- "Buffer_2"
    colnames(buffer.input.3)[colnames(buffer.input.3) == "smooth_Reflectance"] <- "Buffer_3"
    
  } # End of importing data if statement
  
  
  # Join together the input dataframes for plotting
  buffer.input.full <- buffer.input.1 %>% 
    dplyr::select(-c(Tile, EASTING, NORTHING, Reflectance, Wavelength, scale_Reflectance, bn_Reflectance, mean_NIR)) %>% 
    left_join(., buffer.input.2, by = c("PLOT" = "PLOT", "Year" = "Year", "Band" = "Band")) %>% 
    dplyr::select(-c(Tile, EASTING, NORTHING, Reflectance, Wavelength, scale_Reflectance, bn_Reflectance, mean_NIR)) %>% 
    left_join(., buffer.input.3, by = c("PLOT" = "PLOT", "Year" = "Year", "Band" = "Band")) %>% 
    dplyr::select(-c(Tile, EASTING, NORTHING, Reflectance, Wavelength, scale_Reflectance, bn_Reflectance, mean_NIR))
  
  # Plot two of the buffers against each other to observe the relationship
  (buffer.plot.1.2 <- ggplot(data = buffer.input.full, aes(x = Buffer_1, y = Buffer_2, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 0.5, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +

      ggpmisc::stat_fit_glance(method = 'lm', # Don't load this package (breaks tidyverse)
                               method.args = list(formula = y ~ x),  geom = 'text',
                               aes(label = paste0("\nR2 = ", signif(..r.squared.., digits = 3)))) +

      scale_fill_viridis(option = "magma", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(buffer.input.full$Buffer_1), max(buffer.input.full$Buffer_1)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(buffer.input.full$Buffer_2), max(buffer.input.full$Buffer_2)), expand = c(0,0)) +
      labs(#title = "Smoothed Reflectance",
           #subtitle = paste0("Buffers: ", buffer.1, " m & ", buffer.2, " m"),
           x = paste0("\nReflectance - Buffer ", buffer.1, " m"),
           y = paste0("Reflectance - Buffer ", buffer.2, " m\n"),
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  (buffer.plot.1.3 <- ggplot(data = buffer.input.full, aes(x = Buffer_1, y = Buffer_3, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 0.5, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      
      ggpmisc::stat_fit_glance(method = 'lm', # Don't load this package (breaks tidyverse)
                               method.args = list(formula = y ~ x),  geom = 'text',
                               aes(label = paste0("\nR2 = ", signif(..r.squared.., digits = 3)))) +
      
      scale_fill_viridis(option = "viridis", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(buffer.input.full$Buffer_1), max(buffer.input.full$Buffer_1)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(buffer.input.full$Buffer_3), max(buffer.input.full$Buffer_3)), expand = c(0,0)) +
      labs(#title = "Smoothed Reflectance",
           #subtitle = paste0("Buffers: ", buffer.1, " m & ", buffer.3, " m"),
           x = paste0("\nReflectance - Buffer ", buffer.1, " m"),
           y = paste0("Reflectance - Buffer ", buffer.3, " m\n"),
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  (buffer.plot.2.3 <- ggplot(data = buffer.input.full, aes(x = Buffer_2, y = Buffer_3, fill = as.factor(Year))) +
      geom_point(shape = 21, size = 2, alpha = 0.5, colour = "#000000") +
      stat_smooth(method = lm, fill = "#D3D3D3", colour = "#000000", alpha = 0.5, se = TRUE) +
      
      ggpmisc::stat_fit_glance(method = 'lm', # Don't load this package (breaks tidyverse)
                               method.args = list(formula = y ~ x),  geom = 'text',
                               aes(label = paste0("\nR2 = ", signif(..r.squared.., digits = 3)))) +
      
      scale_fill_viridis(option = "mako", begin = 0, end = 1, direction = 1, discrete = TRUE) +
      scale_x_continuous(limits = c(min(buffer.input.full$Buffer_2), max(buffer.input.full$Buffer_2)), expand = c(0,0)) +
      scale_y_continuous(limits = c(min(buffer.input.full$Buffer_3), max(buffer.input.full$Buffer_3)), expand = c(0,0)) +
      labs(#title = "Smoothed Reflectance",
           #subtitle = paste0("Buffers: ", buffer.2, " m & ", buffer.3, " m"),
           x = paste0("\nReflectance - Buffer ", buffer.2, " m"),
           y = paste0("Reflectance - Buffer ", buffer.3, " m\n"),
           fill = "Year") +
      theme_1() +
      theme(legend.position = "right"))
  
  # Create a panel of the three plots
  buffer.panel <- grid.arrange(buffer.plot.1.2, buffer.plot.1.3, buffer.plot.2.3, ncol = 3)
  
  # Output plot to .png
  ggsave(buffer.panel, filename = paste0("outputs/figures/beta_spectral_comparison_buffer_b",
                                         buffer.1, "_", buffer.2, "_", buffer.3, filepath.brightness,
                                         filepath.smoothing, filepath.top.hits, ".png"), width = 18, height = 5)
    
  
} # End of function
