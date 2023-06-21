# 01 - Saddle Composition Data
# Joseph Everest
# February 2023, modified March 2023


# LOAD PACKAGES ----

# Load packages
library(tidyverse)
library(Taxonstand)
library(sp)


# **[CHANGE]** - DECIDE WHETHER TO USE ALL HITS OR TOP HITS ONLY ----

# Decision
top.hits.only <- "No" # Default = "No"

# Generate output folder path
if (top.hits.only == "No"){ filepath.top.hits <- "" } else { filepath.top.hits <- "_top_hits_only" }


# IMPORT DATA (saddle.1) ----

# Load saddle composition data
saddle.1 <- read.csv("data/nwt_composition.csv")


# DATA MANIPULATION (saddle.2) ----

# Modify the saddle data to get into most useful format
saddle.2 <- saddle.1 %>% 
  dplyr::select(-c(LTER_site, local_site, USDA_code)) %>% # Unneeded columns
  filter(year %in% c(2013:2021)) %>% # Match up with hyperspectral data
  mutate(PlotYear = paste0(plot, ":", year)) %>% 
  relocate(PlotYear, plot, .before = year)


# ADD PLOT COORDINATES (saddle.3) ----

# Import and manipulate old coordinates
coords.old <- read.csv("data/nwt_lter_locations_scelmendorf.csv") %>% 
  dplyr::select(SITECODE, LATITUDE, LONGITUDE, ALT_SITECODE) %>% # Retain only relevant columns
  filter(str_detect(ALT_SITECODE, pattern = "addle")) %>% # Retain locations with mention of s"addle"
  mutate(SITECODE = str_remove_all(SITECODE, "PTQUAD"), # Remove "PTQUAD" to leave just the plot number
         SITECODE = str_replace_all(SITECODE, "A", "1"), # Replace the 'A' in sitecodes with '1' to get matching sitecode to saddle data
         SITECODE = as.numeric(SITECODE)) %>% 
  dplyr::select(-ALT_SITECODE) %>% # Remove additional site code column
  rename(LAT_old = LATITUDE, LONG_old = LONGITUDE) %>% 
  arrange(SITECODE)

# Cut to just the LAT and LONG columns
latlong.old <- coords.old %>% 
  dplyr::select(LAT_old, LONG_old)

# Transform into spatial points object
coordinates(latlong.old) <- c("LONG_old", "LAT_old")

# Apply CRS to spatial points object
proj4string(latlong.old) <- CRS("+proj=longlat +datum=WGS84")

# Transform spatial points object to UTM
utm.old <- spTransform(latlong.old, CRS("+proj=utm +zone=13 +datum=WGS84"))

# Change back to dataframe to rejoin with saddle data
utm.old.df <- as.data.frame(utm.old) %>%
  rename(EASTING_old = LONG_old, NORTHING_old = LAT_old)

# Join easting and northing data to saddle dataframe
coords.old <- cbind(coords.old, utm.old.df) %>%
  relocate(NORTHING_old, EASTING_old, .after = LONG_old)

# Import new coordinates
coords.new <- read.csv("data/nwt_lter_locations_utm_II_v10.csv") %>% 
  dplyr::select(ALT_SITECODE, UTM_E, UTM_N) %>% 
  filter(str_detect(ALT_SITECODE, pattern = "Pt Quad")) %>% 
  rename(PLOT = ALT_SITECODE,
         EASTING_new = UTM_E,
         NORTHING_new = UTM_N) %>% 
  mutate(PLOT = str_remove(PLOT, pattern = "Saddle Pt Quad "),
         PLOT = str_remove(PLOT, pattern = "^0+"),
         PLOT = str_replace(PLOT, pattern = "A", replacement = "1"),
         PLOT = as.integer(PLOT)) %>% 
  arrange(PLOT)

# Join together and determine differences
coords.full <- left_join(coords.old, coords.new, by = c("SITECODE" = "PLOT")) %>% 
  relocate(LONG_old, LAT_old, NORTHING_old, NORTHING_new, EASTING_old, EASTING_new, .after = SITECODE) %>% 
  mutate(NORTHING_dif = abs(NORTHING_new - NORTHING_old),
         EASTING_dif = abs(EASTING_new - EASTING_old),
         TOTAL_dif = NORTHING_dif + EASTING_dif) %>% 
  arrange(desc(TOTAL_dif))

# Export the locations and differences between old and new measurements
write.csv(coords.full, file = "outputs/output_saddle_locations_full.csv", row.names = FALSE)

# Trim down coordinates for joining with composition data
locations <- coords.full %>% 
  dplyr::select(SITECODE, EASTING_new, NORTHING_new) %>% 
  rename(EASTING = EASTING_new, NORTHING = NORTHING_new) %>% 
  arrange(SITECODE)

# Combine coordinates with saddle data
saddle.2a <- left_join(saddle.2, locations, by = c("plot" = "SITECODE")) %>%
  relocate(EASTING, NORTHING, .after = year)

# Remove plot 3 as large discrepancy in coordinate location and outlier in later analyses
saddle.2b <- saddle.2a %>% 
  filter(plot != 3)

# Modify point coordinates from (5-95 x 5-95) to (1-10 x 1-10)
saddle.3 <- saddle.2b %>% 
  mutate(x = (x + 5)/10, # Converts 5 -> 1, 15 -> 2,... 95 -> 10 etc. 
         y = (y + 5)/10) # As above...


# PLOT TYPE (saddle.4) ----

# Load in saddle NPP dataset which has veg class per plot
veg.class.1 <- read.csv("data/nwt_biomass.csv")

  # NOTE: 81 of 88 plots present
    # MISSING: 1, 13, 41, 42, 63, 71, 301
    # These seven plots are of type 'barren', so not recorded in NPP dataset
    # Added manually based on additional information from Sarah Elmendorf

# Trim down the NPP dataset to useful veg type data
veg.class.2 <- veg.class.1 %>% 
  dplyr::select(grid_pt, year, veg_class) %>% 
  rename(plot = grid_pt)

# Check that vegetation type is the same for each plot across years
veg.class.test <- veg.class.2 %>% 
  group_by(plot) %>% 
  summarise(Types = length(unique(veg_class)))
    # CHECKED: veg type of each plot is consistent across years

# Create dataframe with one veg class per plot
veg.class.3 <- veg.class.2 %>% 
  filter(year == 2020) %>%  # All 81 recorded plots assessed in 2020
  dplyr::select(-year) %>%
  mutate(veg_class = case_when(veg_class %in% c("WM") ~ "wet_meadow",
                               veg_class %in% c("SF") ~ "snowfence",
                               veg_class %in% c("FF") ~ "fellfield",
                               veg_class %in% c("DM") ~ "dry_meadow",
                               veg_class %in% c("MM") ~ "moist_meadow",
                               veg_class %in% c("ST") ~ "shrub_tundra",
                               veg_class %in% c("SB") ~ "snowbed"))

# Join plot vegetation type to main saddle data frame and modify NAs and snowfence
saddle.4 <- left_join(saddle.3, veg.class.3, by = c("plot" = "plot")) %>% 
  mutate(veg_class = ifelse(is.na(veg_class), "barren", veg_class)) %>% # Replace the 7 NAs with 'barren'
  filter(!veg_class %in% c("snowfence")) %>% # Remove the two plots that are snowfence plots (28 and 38)
  relocate(veg_class, .after = NORTHING)


# TAXONOMY CHECK (saddle.5) ----

# Create vector of unique species in saddle dataset
species.saddle <- unique(saddle.4$USDA_name)

# # Run taxonomy checker on full dataset to check species names
# species.checked <- TPL(species.saddle)
# 
# # Write csv of taxonomy checker results
# write.csv(species.checked, "outputs/output_tpl_saddle.csv")

# Read csv of taxonomy checker results and trim dataframe
species.checked <- read.csv("outputs/output_tpl_saddle.csv") %>% 
  dplyr::select(Taxon, New.Genus, New.Species, New.Taxonomic.status, Family) %>%
  mutate(Name_TPL = paste(New.Genus, New.Species, sep = " ")) %>%
  relocate(Name_TPL, New.Genus, Family, .before = New.Taxonomic.status) %>%
  dplyr::select(-New.Species) %>% 
  rename(New.Species = Name_TPL, New.Family = Family) %>% 
  arrange(Taxon)

# Create a vector of non-vascular plant material 'Taxons'
non.vasc.plant <- c("Bare ground", "Bare Ground", "Elk scat", "Hole", "Lichen", "LIchen", "Litter", "Moss",
                    "Rock fragments", "Rock, fragments", "Unknown soil crust", "Unknown composite")
  
# Modify the taxonomy checked dataframe - SPECIES
species.checked.s <- species.checked %>% 
  mutate(New.Species = case_when(Taxon == "Antennaria" ~ "XXXANTENNARIA:SADDLE",
                                 Taxon %in% non.vasc.plant ~ "NON-VASC. PLANT MATERIAL",
                                 Taxon == "Carex" ~ "XXXCAREX:SADDLE",
                                 Taxon == "Cerastium" ~ "XXXCERASTIUM:SADDLE",
                                 Taxon == "Draba" ~ "XXXDRABA:SADDLE",
                                 Taxon == "Erigeron" ~ "XXXERIGERON:SADDLE",
                                 Taxon %in% c("Forb", "Forb (herbaceous, not grass nor grasslike)") ~ "XXXFORB:SADDLE",
                                 Taxon %in% c("Graminoid", "Graminoid (grass or grasslike)") ~ "XXXGRAMINOID:SADDLE",
                                 Taxon == "Juncus" ~ "XXXJUNCUS:SADDLE",
                                 Taxon == "Poa" ~ "XXXPOA:SADDLE",
                                 str_detect(Taxon, "Unknown Carex spp") ~ "XXXCAREX:SADDLE",
                                 TRUE ~ New.Species))

# Modify the taxonomy checked dataframe - GENUS
species.checked.g <- species.checked.s %>% 
  mutate(New.Genus = case_when(Taxon == "Antennaria" ~ "Antennaria",
                               Taxon %in% non.vasc.plant ~ "NA",
                               Taxon == "Carex" ~ "Carex",
                               Taxon == "Cerastium" ~ "Cerastium",
                               Taxon == "Draba" ~ "Draba",
                               Taxon == "Erigeron" ~ "Erigeron",
                               Taxon %in% c("Forb", "Forb (herbaceous, not grass nor grasslike)") ~ "NA",
                               Taxon %in% c("Graminoid", "Graminoid (grass or grasslike)") ~ "NA",
                               Taxon == "Juncus" ~ "Juncus",
                               Taxon == "Poa" ~ "Poa",
                               str_detect(Taxon, "Unknown Carex spp") ~ "Carex",
                               TRUE ~ New.Genus))

# Modify the taxonomy checked dataframe - FAMILY
species.checked.f <- species.checked.g %>% 
  mutate(New.Family = case_when(Taxon == "Antennaria" ~ "Compositae",
                               Taxon %in% non.vasc.plant ~ "NA",
                               Taxon == "Carex" ~ "Cyperaceae",
                               Taxon == "Cerastium" ~ "Caryophyllaceae",
                               Taxon == "Draba" ~ "Brassicaceae",
                               Taxon == "Erigeron" ~ "Compositae",
                               Taxon %in% c("Forb", "Forb (herbaceous, not grass nor grasslike)") ~ "NA",
                               Taxon %in% c("Graminoid", "Graminoid (grass or grasslike)") ~ "NA",
                               Taxon == "Juncus" ~ "Juncaceae",
                               Taxon == "Poa" ~ "Poaceae",
                               str_detect(Taxon, "Unknown Carex spp") ~ "Cyperaceae",
                               TRUE ~ New.Family))

# Modify the taxonomy checked dataframe - TAXONOMIC STATUS
species.checked.t <- species.checked.f %>% 
  mutate(New.Taxonomic.status = case_when(Taxon == "Antennaria" ~ "Modified",
                                          Taxon %in% non.vasc.plant ~ "Modified",
                                          Taxon == "Carex" ~ "Modified",
                                          Taxon == "Cerastium" ~ "Modified",
                                          Taxon == "Draba" ~ "Modified",
                                          Taxon == "Erigeron" ~ "Modified",
                                          Taxon %in% c("Forb", "Forb (herbaceous, not grass nor grasslike)") ~ "Modified",
                                          Taxon %in% c("Graminoid", "Graminoid (grass or grasslike)") ~ "Modified",
                                          Taxon == "Juncus" ~ "Modified",
                                          Taxon == "Poa" ~ "Modified",
                                          str_detect(Taxon, "Unknown Carex spp") ~ "Modified",
                                          TRUE ~ New.Taxonomic.status)) %>% 
  mutate(New.Genus = ifelse(New.Genus == "NA", NA, New.Genus),
         New.Family = ifelse(New.Family == "NA", NA, New.Family))

# Join new species names to the saddle dataframe
saddle.5 <- left_join(saddle.4, species.checked.t, by = c("USDA_name" = "Taxon")) %>% 
  dplyr::select(-New.Taxonomic.status) %>% 
  rename(SPECIES = New.Species, GENUS = New.Genus, FAMILY = New.Family)


# ADD GROWTH FORM (saddle.6) ----

# Add dataset of growth forms for species
growth.forms.1 <- read.csv("data/itex_growth-forms.csv", stringsAsFactors = FALSE) %>% 
  dplyr::select(-X)

# Create dataframe containing ITEX site genus, species, and growth forms
growth.forms.2 <- growth.forms.1 %>% 
  mutate(Name = paste(GENUS, SPECIES, sep = " ")) %>% 
  dplyr::select(Name, GFNARROWwalker)

# Check which species have multiple growth forms
growth.forms.multi <- growth.forms.2 %>% 
  unique() %>% 
  group_by(Name) %>% 
  summarise(GFs = length(unique(GFNARROWwalker))) %>% 
  ungroup() %>% 
  filter(GFs > 1) # Species that are listed as having more than one growth form

# Export list of species with multiple growth forms
growth.forms.multi.sp <- unique(growth.forms.multi$Name)

# See if any of the saddle species have multiple growth forms (4 species do)
growth.forms.multi.saddle <- saddle.5 %>% 
  filter(SPECIES %in% growth.forms.multi.sp) %>% 
  dplyr::select(SPECIES) %>% 
  unique()
    # TWO SPECIES WITH MULTIPLE GF: Sibbaldia procumbens & Kobresia myosuroides

# Export list of those two saddle species with multiple growth forms
growth.forms.multi.saddle.sp <- unique(growth.forms.multi.saddle$SPECIES)

# Determine whichs growth forms these 2 species already have:
growth.forms.3 <- growth.forms.2 %>% 
  filter(Name %in% growth.forms.multi.saddle.sp) %>% 
  dplyr::select(Name, GFNARROWwalker) %>% 
  unique()
    # Kobresia myosuroides - different types of graminoid, will get resolved with simplification
    # Sibbaldia procumbens - forb / subshrub, appears to be mostly accepted as a forb

# Modify growth forms of these species in the growth forms dataset
growth.forms.4 <- growth.forms.2 %>% 
  mutate(GFNARROWwalker = ifelse(Name %in% c("Sibbaldia procumbens"), "FORB", GFNARROWwalker),
         GFNARROWwalker = ifelse(Name %in% c("Kobresia myosuroides"), "GRAMINOID", GFNARROWwalker)) %>% 
  distinct()

# Join growth forms to saddle.5 data and simplify growth form names
saddle.5a <- left_join(saddle.5, growth.forms.4, by = c("SPECIES" = "Name")) %>% 
  rename(GROWTH_FORM = GFNARROWwalker) %>% 
  mutate(GROWTH_FORM = case_when(GROWTH_FORM %in% c("GRASS", "RUSH", "SEDGE") ~ "GRAMINOID",
                                 GROWTH_FORM == "SDECI" ~ "SHRUB",
                                 TRUE ~ GROWTH_FORM))

# See which species are missing growth form data
growth.form.missing.saddle <- saddle.5a %>% 
  dplyr::select(SPECIES, GROWTH_FORM) %>% 
  unique() %>% 
  filter(is.na(GROWTH_FORM)) %>%
  arrange(SPECIES)

# Modify missing growth form data for saddle species
saddle.6 <- saddle.5a %>% 
  mutate(GROWTH_FORM = case_when(str_detect(SPECIES, "Carex") ~ "GRAMINOID",
                                 SPECIES == "Comastoma tenellum" ~ "FORB",
                                 SPECIES == "Cymopterus alpinus" ~ "FORB",
                                 SPECIES == "Deschampsia cespitosa" ~ "GRAMINOID",
                                 SPECIES == "Elymus trachycaulus" ~ "GRAMINOID",
                                 SPECIES == "Eremogone fendleri" ~ "FORB",
                                 SPECIES == "Eritrichium nanum" ~ "FORB",
                                 SPECIES == "Erysimum capitatum" ~ "FORB",
                                 SPECIES == "Gagea serotina" ~ "FORB",
                                 SPECIES == "Gentiana algida" ~ "FORB",
                                 SPECIES == "Geum rossii" ~ "FORB",
                                 SPECIES == "Hieracium froelichianum" ~ "FORB",
                                 SPECIES == "Packera cana" ~ "FORB",
                                 SPECIES == "Packera crocata" ~ "FORB",
                                 SPECIES == "Packera werneriifolia" ~ "FORB",
                                 SPECIES == "Persicaria vivipara" ~ "FORB",
                                 SPECIES == "Poa rupicola" ~ "GRAMINOID",
                                 SPECIES == "Ranunculus eschscholtzii" ~ "FORB",
                                 SPECIES == "Salix petrophila" ~ "SHRUB",
                                 SPECIES == "Selaginella scopulorum" ~ "FORB",
                                 SPECIES == "Solidago simplex" ~ "FORB",
                                 SPECIES == "Tetraneuris acaulis" ~ "FORB",
                                 SPECIES == "XXXANTENNARIA:SADDLE" ~ "FORB",
                                 SPECIES == "XXXCAREX:SADDLE" ~ "GRAMINOID",
                                 SPECIES == "XXXCERASTIUM:SADDLE" ~ "FORB",
                                 SPECIES == "XXXDRABA:SADDLE" ~ "FORB",
                                 SPECIES == "XXXFORB:SADDLE" ~ "FORB",
                                 SPECIES == "XXXGRAMINOID:SADDLE" ~ "GRAMINOID",
                                 SPECIES == "XXXJUNCUS:SADDLE" ~ "GRAMINOID",
                                 SPECIES == "XXXPOA:SADDLE" ~ "GRAMINOID",
                                 TRUE ~ GROWTH_FORM))


# PIN EDITS (saddle.7) ----

# Correct for pins that have 'extra' hits
saddle.6a <- saddle.6 %>% 
  filter(hit_type != "extra")

# Check to make sure no species occurs twice on same pin (one duplicate)
pin.duplicate <- saddle.6a %>% 
  filter(SPECIES != "NON-VASC. PLANT MATERIAL") %>% 
  group_by(PlotYear, x, y) %>% 
  summarise(total.pin = length(SPECIES),
            sp.pin = length(unique(SPECIES))) %>% 
  ungroup() %>%
  mutate(SAME = ifelse(total.pin == sp.pin, TRUE, FALSE)) %>% 
  filter(SAME == FALSE) %>% 
  mutate(PlotYearXY = paste0(PlotYear, ":", x, ":", y))

# Create vector of IDs of duplicate pins
pin.duplicate.ID <- pin.duplicate$PlotYearXY

# Create dataframe of those records
pin.duplicate.2 <- saddle.6a %>% 
  mutate(PlotYearXY = paste0(PlotYear, ":", x, ":", y)) %>% 
  filter(PlotYearXY %in% pin.duplicate.ID)
    # One species: occurred as a top and middle hit so okay to leave in

# Check each plot has 100 hits
pin.hits <- saddle.6a %>% 
  group_by(PlotYear) %>% 
  mutate(hit_count = length(unique(point))) %>% 
  ungroup() %>% 
  filter(hit_count != 100)
    # All have 100 hits

# No more edits required
saddle.7 <- saddle.6a


# PLOT CHECKS (saddle.8) ----

# Check that we have the same number of plots in every year
top.hits.checks <- saddle.7 %>% 
  group_by(year) %>% 
  summarise(num_plots = length(unique(plot))) %>% # 86 in all years
  ungroup()

# Carry over saddle.7 as saddle.8
saddle.8 <- saddle.7


# FINAL TIDY (saddle.9) ----

# Re-add in columns for site and subsite etc. and cut to 2017-2020
saddle.9 <- saddle.8 %>% 
  mutate(SITE = "NIWOT",
         SUBSITE = "SADDLE_GRID") %>% 
  rename(PLOT = plot, YEAR = year, CLASS = veg_class, FuncGroup = GROWTH_FORM) %>% 
  relocate(SITE, SUBSITE, PLOT, YEAR, .before = PlotYear) %>% 
  filter(YEAR %in% c(2017, 2018, 2019, 2020))


# DATAFRAME FOR TOP HITS ONLY ----

# Run loop for if top.hits.only is "YES"
if (top.hits.only == "Yes"){
  
  # Convert all pins with only a single hit to 'top hits'
  saddle.top.hits.only <- saddle.9 %>%
    group_by(PlotYear, x, y) %>%
    mutate(num_hits = length(hit_type)) %>%
    ungroup() %>%
    mutate(hit_type = ifelse(num_hits == 1, "top", hit_type)) %>%
    dplyr::select(-num_hits) %>% 
    filter(hit_type == "top")
  
  # Check that each plot has exactly 100 hits
  saddle.top.hits.checks <- saddle.top.hits.only %>% 
    group_by(PlotYear) %>% 
    summarise(count = length(hit_type)) %>% 
    ungroup()
  
  # All correct so assign object forward
  saddle.10 <- saddle.top.hits.only
  
} else {
  
  # Carry forward original saddle.9
  saddle.10 <- saddle.9
  
} # End of top.hits.only if statement


# EXPORT CLEANED SADDLE DATAFRAME ----

# Export dataframe to csv
write.csv(saddle.10, file = paste0("outputs/output_saddle_composition", filepath.top.hits, ".csv"),
          row.names = FALSE)


# CALCULATE COVER (saddle.9) ----

# Remove any rows that are classes as 'non-vasc. plant material'
cover.vasc <- saddle.10 %>% 
  filter(SPECIES != "NON-VASC. PLANT MATERIAL")

# Calculate cover for each plot
cover.calc <- cover.vasc %>% 
  group_by(PlotYear) %>% 
  mutate(total_hits = length(point)) %>% 
  ungroup() %>% 
  group_by(PlotYear, SPECIES) %>% 
  mutate(species_hits = length(point)) %>% 
  ungroup() %>% 
  mutate(RelativeCover = (species_hits / total_hits) * 100)

# Transform dataframe so one row per species/plot/year
cover.cut <- cover.calc %>% 
  dplyr::select(-c(point, x, y, hit_type, USDA_name)) %>% 
  distinct(SITE, SUBSITE, PLOT, YEAR, NORTHING, EASTING, CLASS, SPECIES,
           GENUS, FAMILY, FuncGroup, RelativeCover, .keep_all = TRUE)

# Check that relative cover in every plot/year sums to 100
cover.sum <- cover.cut %>% 
  group_by(PlotYear) %>% 
  summarise(TotalCover = sum(RelativeCover)) %>% 
  ungroup()
    # All sum to 100

# Determine the cover per functional group
cover.fg <- cover.cut %>% 
  group_by(PlotYear, FuncGroup) %>% 
  mutate(FGCover = sum(RelativeCover)) %>% 
  ungroup()
  
# Create a key for joining plot-dominating FG information to saddle dataframe
cover.fg.key <- cover.fg %>% 
  dplyr::select(PlotYear, FuncGroup, FGCover) %>% # Only columns we need for calculating plot-dominating FG
  distinct(PlotYear, FuncGroup, .keep_all = TRUE) %>% # Keeps one record of each type of FG in each plot:year combo
  pivot_wider(names_from = FuncGroup, values_from = FGCover, values_fill = list(FGCover = 0)) %>% # Fill missing values with zero
  mutate(PlotDominatingFG = case_when(SHRUB > 50 ~ "Shrub-Dominated",
                                      GRAMINOID > 50 ~ "Graminoid-Dominated",
                                      FORB > 50 ~ "Forb-Dominated",
                                      TRUE ~ "None")) %>%
  rename(ShrubCover = SHRUB, GraminoidCover = GRAMINOID, ForbCover = FORB) %>%
  dplyr::select(PlotYear, GraminoidCover, ShrubCover, ForbCover, PlotDominatingFG)

# Rejoin FG cover and dominating information to main saddle dataframe
cover.final <- left_join(cover.cut, cover.fg.key, by = c("PlotYear" = "PlotYear"))

# Tidy up dataframe
saddle.11 <- cover.final %>% 
  dplyr::select(-c(total_hits, species_hits)) %>% 
  group_by(PlotYear) %>% 
  mutate(SpeciesRichness = length(unique(SPECIES))) %>% 
  ungroup() %>% 
  relocate(SpeciesRichness, .before = RelativeCover)


# EXPORT COVER SADDLE DATAFRAME ----

# Export dataframe to csv
write.csv(saddle.11, file = paste0("outputs/output_saddle_cover", filepath.top.hits, ".csv"),
          row.names = FALSE)
