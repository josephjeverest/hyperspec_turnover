# 02 - Saddle Trait Data
# Joseph Everest
# February 2023, modified September 2023


# LOAD PACKAGES & THEMES ----

# Load packages
library(tidyverse)
library(viridis)
library(Taxonstand)

# Load themes
source("scripts/EX1_ggplot_themes.R")


# **[CHANGE]** - DECIDE WHETHER TO USE TOP HITS ONLY OR NOT ----

# Decision
top.hits.only <- "No" # Default = "No"

# Generate output folder path
if (top.hits.only == "No"){ filepath.top.hits <- "" } else { filepath.top.hits <- "_top_hits_only" }


# IMPORT DATA (traits.1) ----

# Load Niwot trait data
traits.1 <- read.csv("data/plant_trait.ms.data.csv")


# DATA MANIPULATION (traits.2) ----

# Create vectors of columns to remove
remove.experiment <- c("Species", "Family", "USDA.Code", "Code", "SoilMoisture")
remove.traits <- c("VegHeight", "CC1", "CC2", "CC3", "CC.Relative", "Stomatal.Conductance",
                   "LeafArea", "WetWeight", "DryWeight", "LWC", "LMA", "CN_ratio")

# Modify the trait data to get into most useful format
traits.2 <- traits.1 %>% 
  dplyr::select(-remove.experiment, -remove.traits) %>% # Remove processing
  mutate(Site_Rep = paste0(Exp, "_", Rep)) %>% 
  dplyr::select(-Rep) %>% 
  rename(Original_Name = Latin.name,
         Original_FuncGroup = Life.Form,
         TREATMENT = TRT,
         LAT = Lat, LONG = Long,
         SITE = Exp,
         YEAR = YearCollected,
         Height = OHeight,
         Chlorophyll = CCAvg,
         LeafN = Percent_N,
         LeafC = Percent_C) %>% 
  relocate(SITE, Site_Rep, .after = YEAR)


# REMOVE UNWANTED TREATMENTS & VALUES (traits.3) ----

# SITES: "MUH", "SAD", "ITEX", "EASTKNOLL", "NaN" & ""
  # For "SAD", "EASTKNOLL", "NaN" & "": retain all treatments
  # For "ITEX": retain only XXX as everything else is treated
  # For "MUH": ...

# Filter out unwanted treatments & early chlorphyll values (inconsistent methods)
traits.3 <- traits.2 %>% 
  mutate(REMOVE = "NO",
         REMOVE = ifelse(SITE == "ITEX" & TREATMENT != "XXX", "YES", "NO")) %>% 
  filter(REMOVE == "NO") %>% # Retain only records to keep
  mutate(CLASS = case_when(TREATMENT %in% c("DRY") ~ "dry_meadow",
                               TREATMENT %in% c("FF") ~ "fellfield",
                               TREATMENT %in% c("XXX", "MOIST") ~ "moist_meadow",
                               TREATMENT %in% c("SNOWBED") ~ "snowbed",
                               TREATMENT %in% c("WET") ~ "wet_meadow")) %>% # Rename treatments for matching to classifications in composition data
  dplyr::select(-c(REMOVE, SITE, TREATMENT)) %>% 
  relocate(CLASS, .after = LONG) %>% 
  mutate(Chlorophyll = ifelse(YEAR == "2008_2009", NaN, Chlorophyll))


  # TAXONOMY CHECKS (traits.4) ----

# Create vector of unique species in trait dataset
species.traits <- unique(traits.3$Original_Name)

# # Run taxonomy checker on full dataset to check species names
# species.checked <- TPL(species.traits)
# 
# # Write csv of taxonomy checker results
# write.csv(species.checked, "outputs/output_tpl_traits.csv")

# Read csv of taxonomy checker results and trim dataframe
species.checked <- read.csv(paste0("outputs/output_tpl_traits.csv")) %>% 
  dplyr::select(Taxon, New.Genus, New.Species, New.Taxonomic.status, Family) %>%
  mutate(Name_TPL = paste(New.Genus, New.Species, sep = " ")) %>%
  relocate(Name_TPL, New.Genus, Family, .before = New.Taxonomic.status) %>%
  dplyr::select(-New.Species) %>% 
  rename(New.Species = Name_TPL, New.Family = Family) %>% 
  arrange(Taxon)

# Modify the one unresolved name in the dataset
species.resolved <- species.checked %>% 
  mutate(New.Species = ifelse(Taxon == "Poa spp. tall", "XXXPOA:NIWOT", New.Species),
         New.Genus = ifelse(Taxon == "Poa spp. tall", "Poa", New.Genus),
         New.Family = ifelse(Taxon == "Poa spp. tall", "Poaceae", New.Family),
         New.Taxonomic.status = ifelse(Taxon == "Poa spp. tall", "Modified", New.Taxonomic.status)) %>% 
  arrange(New.Species) %>% 
  dplyr::select(-New.Taxonomic.status) %>% 
  rename(SPECIES = New.Species, GENUS = New.Genus, FAMILY = New.Family)

# Join back together with trait data
traits.4 <- left_join(traits.3, species.resolved, by = c("Original_Name" = "Taxon")) %>% 
  relocate(SPECIES, GENUS, FAMILY, .after = )


# FUNCTIONAL GROUP CHECKS (traits.5) ----

# Observe funcgroups currently listed in trait dataframe
unique(traits.4$Original_FuncGroup)

# Generate vectors to simplify functional groups
fg.forb <- c("Forb", "NFixing.Forb", "Forb,Herb", "Forb,herb")
fg.graminoid <- c("Graminoid")
fg.shrub <- c("Subshrub", "Shrub,Tree", "Shrub", "Shrub,Subshrub,Tree", "Shrub,Subshrub",
              "Subshrub, Shrub")
fg.forb.subshrub <- c("Subshrub, Forb", "Subshrub, Shrub, Forb", "Forb,herb,Subshrub",
                      "Forb, Subshrub", "Nfixing.Subshrub, Forb")

# Resolve functional group names
fg.simple <- traits.4 %>% 
  mutate(FuncGroup = case_when(Original_FuncGroup %in% fg.forb ~ "FORB",
                               Original_FuncGroup %in% fg.graminoid ~ "GRAMINOID",
                               Original_FuncGroup %in% fg.shrub ~ "SHRUB",
                               Original_FuncGroup %in% fg.forb.subshrub ~ "FORB, SHRUB",
                               TRUE ~ Original_FuncGroup)) %>% 
  relocate(FuncGroup, .after = FAMILY) %>% 
  dplyr::select(-Original_FuncGroup)

# Determine the species in the "Forb, Shrub" group
fg.forbshrub <- fg.simple %>% 
  filter(FuncGroup == "FORB, SHRUB") %>% 
  dplyr::select(SPECIES) %>% 
  unique() %>% 
  arrange(SPECIES)

# Manually correct the "Forb, Shrub" species
traits.5 <- fg.simple %>% 
  mutate(FuncGroup = case_when(SPECIES == "Agoseris aurantiaca" ~ "FORB",
                               SPECIES == "Campanula uniflora" ~ "FORB",
                               SPECIES == "Castilleja occidentalis" ~ "FORB",
                               SPECIES == "Cryptantha cana" ~ "FORB",
                               SPECIES == "Dryas octopetala" ~ "SHRUB",
                               SPECIES == "Eriogonum jamesii" ~ "SHRUB",
                               SPECIES == "Erysimum capitatum" ~ "FORB",
                               SPECIES == "Minuartia obtusiloba" ~ "FORB",
                               SPECIES == "Packera cana" ~ "FORB",
                               SPECIES == "Pedicularis racemosa" ~ "FORB",
                               SPECIES == "Penstemon whippleanus" ~ "FORB",
                               SPECIES == "Potentilla concinna" ~ "FORB",
                               SPECIES == "Potentilla nivea" ~ "FORB",
                               SPECIES == "Senecio triangularis" ~ "FORB",
                               SPECIES == "Trifolium nanum" ~ "FORB",
                               TRUE ~ FuncGroup))
                                  # Each species has one unique functional group, checked!


# DUPLICATES CHECK (traits.6) ----

# Remove duplicate rows
traits.6 <- traits.5 %>% 
  distinct(SPECIES, YEAR, Site_Rep, LAT, LONG, Height, Chlorophyll, LDMC, SLA, D15N, D13C, LeafN, LeafC, .keep_all = TRUE) %>% 
  dplyr::select(-Site_Rep)


# EXPORT CLEANED TRAIT DATAFRAME ----

# Export dataframe to csv
write.csv(traits.6, file = paste0("outputs/output_traits_full.csv"), row.names = FALSE)


# SIMPLIFY TRAIT DATAFRAME FOR CALCULATING AVERAGES (traits.7) ----

# Generate a dataframe with only required columns
traits.7 <- traits.6 %>% 
  dplyr::select(-c(Original_Name, YEAR, LAT, LONG)) %>% 
  pivot_longer(names_to = "Trait_Name", values_to = "Trait_Value",
               cols = c("Height", "SLA", "LDMC", "Chlorophyll", "D15N", "D13C", "LeafN", "LeafC")) %>% 
  mutate(Trait_Value = ifelse(is.nan(Trait_Value), NA, Trait_Value)) %>% 
  mutate(ID = row.names(.)) %>% 
  relocate(ID, .before = )


# LOAD IN COMPOSITION DATA & CUT (saddle.1) ----

# Load in composition data
saddle.1 <- read.csv(paste0("outputs/output_saddle_cover", filepath.top.hits, ".csv")) %>% 
  mutate(ID = row.names(.)) %>% 
  relocate(ID, .before = ) 


# CUT COMPOSITION DATA (saddle.2) ----

saddle.2 <- saddle.1 %>% 
  dplyr::select(ID, SPECIES, GENUS, FAMILY, FuncGroup, CLASS)


# PREPARE DATAFRAMES FOR GAP-FILLING (saddle.2 / traits.8) ----

# Generate unique identifiers for all the levels we want to gap-fill by for composition
saddle.3 <- saddle.2 %>% 
  mutate(sp_class = paste0(SPECIES, "_", CLASS),
         sp = SPECIES,
         g_class = paste0(GENUS, "_", CLASS),
         g = GENUS,
         f_class = paste0(FAMILY, "_", CLASS),
         f = FAMILY,
         fun_class = paste0(FuncGroup, "_", CLASS),
         fun = FuncGroup) %>% 
  dplyr::select(-CLASS) %>% 
  pivot_longer(names_to = "Unit_Type", values_to = "Unit", cols = 6:ncol(.)) %>% 
  mutate(Unit = ifelse(str_detect(SPECIES, pattern = "XXX") & str_detect(Unit_Type, pattern = "sp"), NA, Unit), # Species morphpspecies
         Unit = ifelse(str_detect(SPECIES, pattern = "XXX") & is.na(GENUS) & str_detect(Unit_Type, pattern = "g_"), NA, Unit), # Genus morphospecies
         Unit = ifelse(str_detect(SPECIES, pattern = "XXX") & is.na(FAMILY) & str_detect(Unit_Type, pattern = "f_"), NA, Unit)) %>% # Family morphospecies
  dplyr::select(-c(GENUS, FAMILY, FuncGroup))

# Generate unique identifiers for all the levels we want to gap-fill by for traits
traits.8 <- traits.7 %>% 
  mutate(sp_class = paste0(SPECIES, "_", CLASS),
         sp = SPECIES,
         g_class = paste0(GENUS, "_", CLASS),
         g = GENUS,
         f_class = paste0(FAMILY, "_", CLASS),
         f = FAMILY,
         fun_class = paste0(FuncGroup, "_", CLASS),
         fun = FuncGroup) %>% 
  dplyr::select(-CLASS) %>% 
  pivot_longer(names_to = "Unit_Type", values_to = "Unit", cols = 8:ncol(.)) %>% 
  mutate(Trait_Value = ifelse(str_detect(Unit, pattern = "_NA"), NA, Trait_Value)) %>% 
  # mutate(Unit = ifelse(is.na(GENUS) | is.na(FAMILY) | is.na(FuncGroup), NA, Unit)) %>% # If no genus/family/fg, assigning unit as NA
  mutate(Trait_Value = ifelse(str_detect(SPECIES, pattern = "XXX") & str_detect(Unit_Type, pattern = "sp"), NA, Trait_Value), # Species morphpspecies
         Trait_Value = ifelse(str_detect(SPECIES, pattern = "XXX") & is.na(GENUS) & str_detect(Unit_Type, pattern = "g_"), NA, Trait_Value), # Genus morphospecies
         Trait_Value = ifelse(str_detect(SPECIES, pattern = "XXX") & is.na(FAMILY) & str_detect(Unit_Type, pattern = "f_"), NA, Trait_Value)) %>% # Family morphospecies
  group_by(Unit, Trait_Name) %>% # Add columns counting number of records for that value type
  mutate(Records = length(Trait_Value)) %>% 
  ungroup() %>% 
  mutate(Records = ifelse(is.na(Unit), NA, Records)) # if unit is NA, assigning Records as NA
  
# Generate error risks per Unit for each trait
traits.error.risk <- traits.8 %>% 
  dplyr::select(ID, Unit, Unit_Type, Trait_Name, Trait_Value, Records) %>% 
  group_by(Unit, Trait_Name) %>% 
  mutate(sd_unit = sd(Trait_Value, na.rm = TRUE),
         mean_unit = mean(Trait_Value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(Error_Risk = (abs(Trait_Value - mean_unit) / sd_unit)) %>% # abs() is absolute value, makes all differences positive
  dplyr::select(-c(sd_unit, mean_unit))

# Generate average trait values per Unit for each trait
traits.medians <- traits.error.risk %>% 
  filter(Records >= 4) %>% # Don't use a unit if 1-3 records
  mutate(Retain = "Remove") %>% # Set default for column to remove
  mutate(Retain = case_when(Records >= 4 & Records < 10 & Error_Risk < 2.25 ~ "Keep", # 4-9 records
                            Records >= 10 & Records < 20 & Error_Risk < 2.75 ~ "Keep", # 10-19 records
                            Records >= 20 & Records < 30 & Error_Risk < 3.25 ~ "Keep", # 20-29 records
                            Records >= 30 & Error_Risk < 4 ~ "Keep", # 30+ records
                            TRUE ~ Retain)) %>%
  filter(Retain == "Keep") %>% # Only keep records that fall within Error Risk stipulations
  dplyr::select(-c(Retain, Records)) %>%
  group_by(Unit, Trait_Name) %>% # Create new column counting number of records for that Unit type now that values outside error risk band have been removed
  mutate(Records = length(Trait_Value)) %>%
  ungroup() %>%
  filter(Records >= 4) %>% # Ensure that Unit still has minimum four units for calculating average trait value
  group_by(Unit, Trait_Name) %>% # Calculating averages by unit (e.g species/genus etc. in certain region etc.)
  summarise(Median_Trait_Value = median(Trait_Value, na.rm = TRUE)) %>% 
  ungroup()


# GAP-FILLING (combo.1) ----

# Create gap-filling function
gap_fill <- function(df){
  mutate(df, accepted_value = ifelse(!is.na(sp_class), sp_class,
                                     ifelse(!is.na(sp), sp,
                                            ifelse(!is.na(g_class), g_class,
                                                   ifelse(!is.na(g), g,
                                                          ifelse(!is.na(f_class), f_class,
                                                                 ifelse(!is.na(f), f,
                                                                        ifelse(!is.na(fun_class), fun_class,
                                                                               ifelse(!is.na(fun), fun, NA)))))))))}

# Now create individual composition dataframes for each trait
saddle.Height <- left_join(saddle.3, filter(traits.medians, Trait_Name == "Height"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(Height = accepted_value) %>% 
  dplyr::select(ID, Height)

saddle.SLA <- left_join(saddle.3, filter(traits.medians, Trait_Name == "SLA"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(SLA = accepted_value) %>% 
  dplyr::select(ID, SLA)

saddle.LDMC <- left_join(saddle.3, filter(traits.medians, Trait_Name == "LDMC"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(LDMC = accepted_value) %>% 
  dplyr::select(ID, LDMC)

saddle.Chlorophyll <- left_join(saddle.3, filter(traits.medians, Trait_Name == "Chlorophyll"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(Chlorophyll = accepted_value) %>% 
  dplyr::select(ID, Chlorophyll)

saddle.D15N <- left_join(saddle.3, filter(traits.medians, Trait_Name == "D15N"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(D15N = accepted_value) %>% 
  dplyr::select(ID, D15N)

saddle.D13C <- left_join(saddle.3, filter(traits.medians, Trait_Name == "D13C"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(D13C = accepted_value) %>% 
  dplyr::select(ID, D13C)

saddle.LeafN <- left_join(saddle.3, filter(traits.medians, Trait_Name == "LeafN"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(LeafN = accepted_value) %>% 
  dplyr::select(ID, LeafN)

saddle.LeafC <- left_join(saddle.3, filter(traits.medians, Trait_Name == "LeafC"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gap_fill(.) %>% 
  rename(LeafC = accepted_value) %>% 
  dplyr::select(ID, LeafC)

# Join the individual trait values together to generate the complete gap-filled trait dataset
combo.1 <- left_join(saddle.Height, saddle.SLA, by = c("ID" = "ID")) %>% 
  left_join(., saddle.LDMC, by = c("ID" = "ID")) %>% 
  left_join(., saddle.Chlorophyll, by = c("ID" = "ID")) %>% 
  left_join(., saddle.D15N, by = c("ID" = "ID")) %>% 
  left_join(., saddle.D13C, by = c("ID" = "ID")) %>% 
  left_join(., saddle.LeafN, by = c("ID" = "ID")) %>% 
  left_join(., saddle.LeafC, by = c("ID" = "ID"))

# JOIN TO COMPOSITION DATAFRAME ----

# Join traits to saddle composition data
combo.2 <- left_join(saddle.1, combo.1, by = c("ID" = "ID"))


# EXPORT COMBINED DATAFRAME ----

# Output as .csv
write.csv(combo.2, file = paste0("outputs/output_saddle_composition_traits", filepath.top.hits, ".csv"), row.names = FALSE)


# IDENTIFYING LEVEL OF TRAIT MEDIANS BY TAXONOMIC LEVEL (prop.1) ----

# Create gap-filling function
gapfill_level <- function(df){
  mutate(df, gapfill_level = ifelse(!is.na(sp_class), "sp_class",
                                    ifelse(!is.na(sp), "sp",
                                           ifelse(!is.na(g_class), "g_class",
                                                  ifelse(!is.na(g), "g",
                                                         ifelse(!is.na(f_class), "f_class",
                                                                ifelse(!is.na(f), "f",
                                                                       ifelse(!is.na(fun_class), "fun_class",
                                                                              ifelse(!is.na(fun), "fun", NA)))))))))}

# Determine the level from which the value was gap-filled
prop.Height <- left_join(saddle.3, filter(traits.medians, Trait_Name == "Height"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gapfill_level(.) %>% 
  rename(Height = gapfill_level) %>% 
  dplyr::select(ID, Height)

prop.SLA <- left_join(saddle.3, filter(traits.medians, Trait_Name == "SLA"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gapfill_level(.) %>% 
  rename(SLA = gapfill_level) %>% 
  dplyr::select(ID, SLA)

prop.LDMC <- left_join(saddle.3, filter(traits.medians, Trait_Name == "LDMC"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gapfill_level(.) %>% 
  rename(LDMC = gapfill_level) %>% 
  dplyr::select(ID, LDMC)

prop.Chlorophyll <- left_join(saddle.3, filter(traits.medians, Trait_Name == "Chlorophyll"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gapfill_level(.) %>% 
  rename(Chlorophyll = gapfill_level) %>% 
  dplyr::select(ID, Chlorophyll)

prop.D15N <- left_join(saddle.3, filter(traits.medians, Trait_Name == "D15N"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gapfill_level(.) %>% 
  rename(D15N = gapfill_level) %>% 
  dplyr::select(ID, D15N)

prop.D13C <- left_join(saddle.3, filter(traits.medians, Trait_Name == "D13C"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gapfill_level(.) %>% 
  rename(D13C = gapfill_level) %>% 
  dplyr::select(ID, D13C)

prop.LeafN <- left_join(saddle.3, filter(traits.medians, Trait_Name == "LeafN"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gapfill_level(.) %>% 
  rename(LeafN = gapfill_level) %>% 
  dplyr::select(ID, LeafN)

prop.LeafC <- left_join(saddle.3, filter(traits.medians, Trait_Name == "LeafC"), by = c("Unit" = "Unit")) %>% 
  dplyr::select(-c(Unit, Trait_Name)) %>% 
  pivot_wider(names_from = "Unit_Type", values_from = "Median_Trait_Value") %>% # Trait values vanishing here (using Achillea millefolium as an example)
  gapfill_level(.) %>% 
  rename(LeafC = gapfill_level) %>% 
  dplyr::select(ID, LeafC)

# Join the individual proportion dataframes together to generate the complete gap-filled trait dataset
prop.1 <- left_join(prop.Height, prop.SLA, by = c("ID" = "ID")) %>% 
  left_join(., prop.LDMC, by = c("ID" = "ID")) %>% 
  left_join(., prop.Chlorophyll, by = c("ID" = "ID")) %>% 
  left_join(., prop.D15N, by = c("ID" = "ID")) %>% 
  left_join(., prop.D13C, by = c("ID" = "ID")) %>% 
  left_join(., prop.LeafN, by = c("ID" = "ID")) %>% 
  left_join(., prop.LeafC, by = c("ID" = "ID"))


# WORK OUT PERCENTAGE OF EACH TRAIT COMPRISED BY EACH LEVEL (prop.2) ----

# Generate percentage of each trait comprised by each taxonomic level and geographic division
prop.2 <- prop.1 %>% 
  pivot_longer(2:ncol(.), names_to = "Trait", values_to = "Data_Type") %>%
  group_by(Trait) %>% 
  mutate(Total_Plots = length(unique(ID))) %>% # Total number of plot records (32,272 each time)
  ungroup() %>% 
  group_by(Trait, Data_Type, Total_Plots) %>% 
  tally() %>%
  ungroup() %>% 
  mutate(Percentage = (n / Total_Plots) * 100) %>% # Percentage of each trait comprised by each level
  dplyr::select(-n, -Total_Plots)

# Generate rows for mean proportion of trait values comprised by each level across all six traits
prop.3 <- prop.2 %>% 
  group_by(Data_Type) %>% 
  summarise(Percentage = mean(Percentage)) %>% # Calculate mean percentage derived from each level
  ungroup() %>% 
  mutate(Trait = "Mean") %>% # Generate column saying that the trait is in fact 'mean'
  rbind(prop.2, .)


# STACKED BAR CHART OF PROPORTIONS (prop.plot) ----

# Rename and reorder factor levels for the Taxonomic hierarchical column (biggest to smallest)
prop.3$Data_Type <- factor(prop.3$Data_Type,
                           levels = c("fun", "fun_class", "f", "f_class", "g", "g_class", "sp", "sp_class"),
                           labels = c("Functional Group", "Function Group by Class",
                                      "Family", "Family by Class",
                                      "Genus", "Genus by Class",
                                      "Species", "Species by Class"))

# Rename and reorder factor levels for the Traits
prop.3$Trait <- factor(prop.3$Trait,
                       levels = c("Height", "SLA", "LDMC", "LeafN", "D15N", "D13C", "LeafC", "Chlorophyll", "Mean"))

# Generate bar plot of the levels at which the data is gap-filled
(prop.plot <- ggplot(data = prop.3, aes(x = Trait, y = Percentage, fill = Data_Type), colour = "black") +
    geom_histogram(stat = "identity", colour = "black", alpha = 0.8) +
    scale_fill_viridis(option = "magma", begin = 0.55, end = 1, direction = -1, discrete = TRUE) +
    labs(title = "Proportion of Gap-filled Trait Data",
         subtitle = "Comprised by each Taxonomic Level",
         x = "\n Trait",
         y = "Proportion of Species Records (%)\n",
         fill = "Gap-fill Level") +
    theme_1() +
    theme(legend.position = "right",
          axis.text.x = element_text(face = c("plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "bold"))))

# Export plot
ggsave(prop.plot, filename = paste0("outputs/figures/saddle_trait_proportions", filepath.top.hits, ".png"), width = 15, height = 7.5)
