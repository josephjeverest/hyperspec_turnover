# Script 0 - FUNCTIONS to download data required for project
# June 2023
# Joseph Everest

# FUNCTION: DOWNLOAD DATA FROM EDI PORTAL ----

download.edi <- function(edi.identifier, entityId, filename){
  
  # View citation for data
  data.citation <- read_data_package_citation(packageId  = edi.identifier)
  
  # View full data package
  data.package <- read_data_entity_names(packageId = edi.identifier)
  
  # Extract the raw data
  data.raw <- read_data_entity(packageId = edi.identifier, entityId = data.package$entityId[entityId])
  
  # Read in the data as .csv
  data <- read_csv(file = data.raw)
  
  # Export .csv to data folder
  write.csv(data, file = paste0("data/", filename), row.names = FALSE)
  
}


# FUNCTION: DOWNLOAD NEON HYPERSPECTRAL TILES ----

download.neon.hyperspec <- function(){
  
  # Create directory to move files to
  dir.create("data/NEON/tiles")
  
  # Import saddle locations
  saddle.locations.raw <- read.csv("data/nwt_lter_locations_utm_II_v10.csv")
  
  # Clean up locations dataframe to retain only required information
  saddle.locations <- saddle.locations.raw %>% 
    dplyr::select(ALT_SITECODE, UTM_E, UTM_N) %>% 
    filter(str_detect(ALT_SITECODE, pattern = "Pt Quad")) %>% 
    rename(PLOT = ALT_SITECODE,
           EASTING = UTM_E,
           NORTHING = UTM_N) %>% 
    mutate(PLOT = str_remove(PLOT, pattern = "Saddle Pt Quad "),
           PLOT = str_remove(PLOT, pattern = "^0+"),
           PLOT = str_replace(PLOT, pattern = "A", replacement = "1"),
           PLOT = as.integer(PLOT)) %>% 
    arrange(PLOT)
  
  # Define Eastings and Northings
  saddle.eastings <- c(unique(saddle.locations$EASTING))
  saddle.northings <- c(unique(saddle.locations$NORTHING))
  
  
  # Create counter of Niwot revisit for determining file renaming (2017 was revisit 1, 2018 was 2...)
  niwot.revisit <- 0
  
  # Run for loop for each year
  for (year in c("2017", "2018", "2019", "2020")){
    
    # Incresae Niwot revisit counter by 1
    niwot.revisit <- niwot.revisit + 1
    
    # Run download function (has to be run on separate years, can't vector them)
    byTileAOP("DP3.30006.001",
              site = "NIWO",
              year = year, # Data available from 2017 to 2020
              easting = saddle.eastings,
              northing = saddle.northings,
              buffer = 1,
              check.size = TRUE,
              savepath = "data/NEON/BOB",
              token = NA_character_)
    
    # Run loop to rename (move) files to easier location with year in filename
    for (tile.extent in c("4433", "4434")){
      
      file.copy(paste0("data/NEON/DP3.30006.001/neon-aop-products/", year, "/FullSite/D13/", year,
                       "_NIWO_", niwot.revisit, "/L3/Spectrometer/Reflectance/NEON_D13_NIWO_DP3_",
                       "449000_", tile.extent, "000_reflectance.h5"),
                paste0("data/NEON/tiles/NEON_D13_NIWO_DP3_", year, "_449000_", tile.extent,
                       "000_reflectance.h5"))
      
    } # End of extent loop
    
  } # End of year loop
  
}


# FUNCTION: DOWNLOAD NEON RGB TILES ----

download.neon.rgb <- function(){
  
  # Create directory to move files to
  dir.create("data/NEON/tiles")
  
  # Import saddle locations
  saddle.locations.raw <- read.csv("data/nwt_lter_locations_utm_II_v10.csv")
  
  # Clean up locations dataframe to retain only required information
  saddle.locations <- saddle.locations.raw %>% 
    dplyr::select(ALT_SITECODE, UTM_E, UTM_N) %>% 
    filter(str_detect(ALT_SITECODE, pattern = "Pt Quad")) %>% 
    rename(PLOT = ALT_SITECODE,
           EASTING = UTM_E,
           NORTHING = UTM_N) %>% 
    mutate(PLOT = str_remove(PLOT, pattern = "Saddle Pt Quad "),
           PLOT = str_remove(PLOT, pattern = "^0+"),
           PLOT = str_replace(PLOT, pattern = "A", replacement = "1"),
           PLOT = as.integer(PLOT)) %>% 
    arrange(PLOT)
  
  # Define Eastings and Northings
  saddle.eastings <- c(unique(saddle.locations$EASTING))
  saddle.northings <- c(unique(saddle.locations$NORTHING))
  
  # Run download function
  byTileAOP("DP3.30010.001",
            site = "NIWO",
            year = "2020", # Data available from 2017 to 2020
            easting = saddle.eastings,
            northing = saddle.northings,
            buffer = 1,
            check.size = TRUE,
            savepath = "data/NEON/",
            token = NA_character_)
  
  # Run loop to rename (move) files to easier location with year in filename
  for (tile.extent in c("4433", "4434")) {
    
    file.copy(paste0("data/NEON/DP3.30010.001/neon-aop-products/2020/",
                     "FullSite/D13/2020_NIWO_4/L3/Camera/Mosaic/2020_NIWO_4_449000_",
                     tile.extent, "000_image.tif"),
              paste0("data/NEON/tiles/2020_NIWO_4_449000_", tile.extent,
                     "000_image.tif"))
    
  } # End of extent loop
  
} # End of years loop
