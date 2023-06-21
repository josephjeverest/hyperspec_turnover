# Script 0 - Download data required for project
# June 2023
# Joseph Everest


  # ** NOTE: ** Only run this script if not previously downloaded the data 


# PACKAGES & FUNCTIONS ----

# Download the packages required
library(tidyverse)
library(EDIutils)
library(neonUtilities)

# Call functions to download EDI and NEON data
source("scripts/00_data_download_FUNCTION.R")


# DATA DOWNLOAD: COMPOSITION ----

# Run data download for composition data
download.edi(edi.identifier = "knb-lter-nwt.93.6", entityId = 1, filename = "nwt_composition.csv")

# View citation
read_data_package_citation(packageId  = "knb-lter-nwt.93.6")


# DATA DOWNLOAD: TRAITS ----

# Run data download for traits data
download.edi(edi.identifier = "knb-lter-nwt.500.3", entityId = 1, filename = "nwt_traits.csv")

# View citation
read_data_package_citation(packageId  = "knb-lter-nwt.500.3")


# DATA DOWNLOAD: BIOMASS ----

# Run data download for biomass data
download.edi(edi.identifier = "knb-lter-nwt.16.7", entityId = 1, filename = "nwt_biomass.csv")

# View citation
read_data_package_citation(packageId  = "knb-lter-nwt.16.7")


# DATA DOWNLOAD: METEOROLOGICAL DATA ----

# Run data download for meteorological data
download.edi(edi.identifier = "knb-lter-nwt.314.1", entityId = 1, filename = "nwt_met-data.csv")

# View citation
read_data_package_citation(packageId  = "knb-lter-nwt.314.1")


# DATA DOWNLOAD: NEON HYPERSPECTRAL TILES ----

# Run data download for NEON hyperspectral tiles
download.neon.hyperspec() # Note: must enter 'y' followed by 'Enter' to download each year's tiles


# DATA DOWNLOAD: NEON RGB TILES ----

# Run data download for NEON rgb tiles
download.neon.rgb() # Note: must enter 'y' followed by 'Enter' to download tiles
