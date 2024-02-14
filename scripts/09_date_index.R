# 09 - Script to check dates that were flown by NEON at Niwot 2017-2020
# Joseph Everest & Claire Lunch (NEON)
# January 2024


# LIBRARIES ----

# Load required packages
library(rhdf5)
library(terra)


# **[CHANGE]** - DECIDE ON PRE-PROCESSING DECISIONS & WHETHER TO USE TOP HITS ONLY OR NOT ----

# Decisions
brightness <- "Yes" # Default = "Yes"
smoothing <- "No" # Default = "No"
top.hits.only <- "No" # Default = "No"
buffer <- 1

# Generate output folder paths
if (brightness == "No"){ filepath.brightness <- "" } else { filepath.brightness <- "_brightness_normalized" }
if (smoothing == "No"){ filepath.smoothing <- "" } else {filepath.smoothing <- "_smoothed"}
if (top.hits.only == "No"){ filepath.top.hits <- "" } else { filepath.top.hits <- "_top_hits_only" }


# IMPORT FUNCTION FOR DATE RASTER ----

# Function to extract date rasters from .h5 files
  # Written by Claire Lunch of the NEON Data Skills team
getDateRaster <- function(filepath, site) {
  
  dataSelectIndex <- h5read(paste(site,
                                  '/Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index',
                                  sep=''), 
                            file=filepath, read.attributes=T)
  dataSelectLookup <- attributes(dataSelectIndex)$Data_Files
  dataSelectLookupSplit <- unlist(strsplit(dataSelectLookup, split=','))
  
  dateLookup <- regmatches(dataSelectLookupSplit,
                           regexpr('20[0-9]{2}[0-9]{4}', 
                                   dataSelectLookupSplit))
  
  epsg <- h5read(paste(site, '/Reflectance/Metadata/Coordinate_System/EPSG Code',
                       sep=''), file=filepath)
  xmin <- as.numeric(regmatches(filepath, regexpr('[0-9]{3}000', filepath)))
  ymin <- as.numeric(regmatches(filepath, regexpr('[0-9]{4}000', filepath)))
  
  dateRaster <- dataSelectIndex
  for(i in 1:ncol(dataSelectIndex)) {
    for(j in 1:nrow(dataSelectIndex)) {
      dateRaster[i,j] <- dateLookup[dataSelectIndex[i,j]]
    }
  }
  
  dateRaster <- rast(dateRaster, crs=paste('EPSG:', epsg, sep=''),
                     extent=ext(c(xmin, xmin+1000, ymin, ymin+1000)))
  
  return(dateRaster)
  
}


# EXRACT NIWOT SAMPLE DATES ----

# Import saddle locations
saddle.locations <- read.csv(paste0("outputs/output_saddle_locations", filepath.top.hits, ".csv"))

# Create key of saddle locations plots to ID
saddle.key <- saddle.locations %>% 
  mutate(id = as.numeric(row.names(.))) %>% 
  relocate(id, .before = )

# Create a vector of all the HDF5 files to work with
tiles.folderpath <- "data/NEON/tiles/"
tiles.filename <- list.files(path = tiles.folderpath, pattern = "*.h5")
tiles.filepath <- paste0(tiles.folderpath, tiles.filename)

# Create blank output dataframe
date.output <- data.frame()

# Write a loop to extract the rasters in turn and save them
for (i in tiles.filepath){

  # Load in the date raster from .h5 file
  date.raster <- getDateRaster(i, "NIWO")
  
  # Trim saddle locations to required columns
  saddle.locations.cut <- dplyr::select(saddle.locations, EASTING, NORTHING)
  
  # Extract dates at each saddle location without specified buffer
  date.extract <- terra::extract(date.raster, saddle.locations.cut, method = "simple", buffer = buffer)
  
  # Generate information for exporting
  raster.year <- str_remove(i, pattern = "data/NEON/tiles/NEON_D13_NIWO_DP3_")
  raster.year <- str_split(raster.year, "\\_", simplify=T)[,1] # Remove everything after "_"
  raster.name <- str_remove(i, "data/NEON/tiles/NEON_D13_NIWO_DP3_") %>% 
    str_remove(., "_reflectance.h5")

  # Add information to output 
  date.tidy <- date.extract %>% 
    mutate(Year = raster.year,
           Raster = raster.name) %>% 
    rename(Extraction_Date = lyr.1) %>% 
    relocate(Year, Raster, Extraction_Date, .after = ID)
  
  # Add to overall output
  date.output <- rbind(date.output, date.tidy)
  
  }

# Tidy up date output
date.output.clean <- date.output %>% 
  filter(!is.na(Extraction_Date)) %>% 
  left_join(., saddle.key, by = c("ID" = "id")) %>% 
  relocate(Year, Raster, PLOT, EASTING, NORTHING, Extraction_Date, .after = ) %>% 
  dplyr::select(-ID) %>% 
  arrange(Year, PLOT)

# Add statistics columns
date.output.stats <- date.output.clean %>% 
  group_by(Year) %>% 
  mutate(Identical = ifelse(length(unique(Extraction_Date)) == 1, TRUE, FALSE)) %>% 
  ungroup() %>% 
  group_by(Year, Extraction_Date) %>% 
  mutate(NumberOfObs = length(PLOT)) %>% 
  ungroup() %>% 
  mutate(PercentageOfObs = (NumberOfObs / length(unique(PLOT)))*100)


# EXPORT EXTRACTED DATES ----

# Export to .csv
write.csv(date.output.stats, file = paste0("outputs/output_saddle_extraction_dates", filepath.top.hits, ".csv"))