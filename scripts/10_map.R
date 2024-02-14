# 10 - Map of Niwot with Veg. Classes
# Joseph Everest
# June 2023


# LOAD PACKAGES & THEMES ----

# Load packages
library(tidyverse)
library(raster)
library(neonUtilities)
library(sf)
library(basemaps)
library(cowplot)


# **[CHANGE]** - WHETHER OR NOT TO DOWNLOAD NEON RGB IMAGERY? ----

# Decide whether to download
download.NEON.tiles <- FALSE # Default = FALSE


# **[CHANGE]** - DECIDE WHETHER TO RETAIN PLOT 37 OR NOT ----

# Decision
retain.37 <- "Yes" # Default = "Yes"

# Generate output folder path
if (retain.37 == "Yes") {
  filepath.37 <- ""
} else {
  filepath.37 <- "_removed_37"
}


# **[CHANGE]** - DECIDE WHETHER TO RETAIN PLOT 37 OR NOT ----

# Decision
top.hits.only <- "No" # Default = "No"

# Generate output folder path
if (top.hits.only == "No") {
  filepath.top.hits <- ""
} else {
  filepath.top.hits <- "_top_hits_only"
}


# DOWNLOAD RGB DATA ----

# Import saddle locations
saddle.locations <- read.csv(paste0("outputs/output_saddle_locations", filepath.37, filepath.top.hits, ".csv"))

# Define Eastings and Northings
saddle.eastings <- c(unique(saddle.locations$EASTING))
saddle.northings <- c(unique(saddle.locations$NORTHING))


if (download.NEON.tiles == TRUE) {
  # Create directory to move files to
  dir.create("data/NEON/tiles")

  # Run download function
  byTileAOP("DP3.30010.001",
    site = "NIWO",
    year = "2020", # Data available from 2017 to 2020
    easting = saddle.eastings,
    northing = saddle.northings,
    buffer = 1,
    check.size = TRUE,
    savepath = "data/NEON/",
    token = NA_character_
  )


  # Run loop to rename (move) files to easier location with year in filename
  for (tile.extent in c("4433", "4434")) {
    file.copy(
      paste0(
        "data/NEON/DP3.30010.001/neon-aop-products/2020/",
        "FullSite/D13/2020_NIWO_4/L3/Camera/Mosaic/2020_NIWO_4_449000_",
        tile.extent, "000_image.tif"
      ),
      paste0(
        "data/NEON/tiles/2020_NIWO_4_449000_", tile.extent,
        "000_image.tif"
      )
    )
  } # End of extent loop
} # End of if statement


# IMPORT VEG. CLASS INFORMATION ----

# Load in saddle NPP dataset which has veg class per plot
veg.class.1 <- read.csv("data/saddgrid_npp.hh.data.csv")

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
  filter(year == 2020) %>% # All 81 recorded plots assessed in 2020
  dplyr::select(-year) %>%
  mutate(veg_class = case_when(
    veg_class %in% c("WM") ~ "wet meadow",
    veg_class %in% c("SF") ~ "snowfence",
    veg_class %in% c("FF") ~ "fellfield",
    veg_class %in% c("DM") ~ "dry meadow",
    veg_class %in% c("MM") ~ "moist meadow",
    veg_class %in% c("ST") ~ "shrub tundra",
    veg_class %in% c("SB") ~ "snowbed"
  ))

# Join plot vegetation type to main saddle data frame and modify NAs and snowfence
saddle.locations.class <- left_join(saddle.locations, veg.class.3, by = c("PLOT" = "plot")) %>%
  relocate(veg_class, .after = NORTHING)

pal <- data.frame(
  colors =
    c("#FCEED3", "#FAED8E", "#F2B76A", "#F28888", "#E05C5C", "#ED71A7E7"),
  veg_class = c(unique(saddle.locations.class$veg_class))
)
saddle.locations.class <- left_join(saddle.locations.class, pal, by = c("veg_class" = "veg_class"))

# IMPORT RASTERS ----
# read metadata
terra::describe("data/NEON/tiles/2020_NIWO_4_449000_4433000_image.tif")

# Import the two required rasters
rgb.raster.1 <- stack("data/NEON/tiles/2020_NIWO_4_449000_4433000_image.tif")
rgb.raster.2 <- stack("data/NEON/tiles/2020_NIWO_4_449000_4434000_image.tif")


# # Mosaic into one full raster
rgb.raster.full <- mosaic(rgb.raster.1, rgb.raster.2, fun = mean)


# Create extent from saddle locations
saddle.extent <- extent(
  min(saddle.locations$EASTING) - 50,
  max(saddle.locations$EASTING) + 50,
  min(saddle.locations$NORTHING) - 150,
  max(saddle.locations$NORTHING) + 50
)

# Crop the raster to the correct extent
rgb.raster.cut <- stack(crop(rgb.raster.full, saddle.extent))

# @joe not sure if anything needs reprojecting at this point
# I assume the saddle locs you hav area already in the same CRS as the NEON?
xy <- saddle.locations.class[, c("EASTING", "NORTHING")]
saddle.locations.class <- SpatialPointsDataFrame(
  coords = xy, data = saddle.locations.class,
  proj4string = CRS(proj4string(rgb.raster.cut))
)

# PRODUCE MAP OF SADDLE PLOTS ----
tiff(paste0("outputs/figures/manuscript/map.tiff"),
  width = 5, height = 6, units = "in", res = 300
)

terra::plotRGB(stack(rgb.raster.cut),
  r = 1, g = 2, b = 3,
  stretch = "lin", mar = 0
)


terra::plot(terra::vect(saddle.locations.class),
  col = pal$colors, y = "veg_class",
  bg = "black", pch = 22, cex = 1.5, add = T,
  legend = "bottomleft", plg = list(
    horiz = TRUE, bg = "white",
    # cex argument doesn't seem to work here
    bty = "o", cex = 0.6
  )
)

terra::sbar(200,
  below = "meters", divs = 4, type = "bar",
  xy = "topright"
)

dev.off()

# PRODUCE MAP OF NIWOT IN CO ----


# approx location of Niwot Ridge
my_points.df <-
  data.frame(
    lon = c(-105.6169466),
    lat = c(40.0597085)
  )

# approx location of Denver for reference
my_points.denver <-
  data.frame(
    lon = c(-104.991531),
    lat = c(39.742043)
  )


co_extent <- extent(c(-109.0489, -102.0424, 36.9949, 41.0006))
co_extent <- as(co_extent, "SpatialPolygons")
sp::proj4string(co_extent) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

my_points.sf <- st_as_sf(my_points.df, coords = c("lon", "lat"), crs = 4326)
my_points.denver <- st_as_sf(my_points.denver, coords = c("lon", "lat"), crs = 4326)

# view all available maps one can use with basemaps
# get_maptypes()

# set defaults for the basemap
# usa_topo_maps esri's look ok
set_defaults(map_service = "esri", map_type = "usa_topo_maps")

# sort out the CRS it's using
get_the_crs <- basemap_raster(co_extent)
my_points.sf <- sf::st_transform(my_points.sf, crs(get_the_crs))
my_points.denver <- sf::st_transform(my_points.denver, crs(get_the_crs))


mymap <- basemap_ggplot(proj_co_extent) +
  geom_sf_text(
    data = my_points.denver, label = "Denver",
    nudge_x = 90000,
    nudge_y = -10000, size = 6
  ) +
  geom_sf(data = my_points.denver, size = 3, color = "black") +
  geom_sf_text(
    data = my_points.sf, label = "Niwot Ridge",
    nudge_x = 180000,
    nudge_y = 15000, size = 8
  ) +
  geom_sf(data = my_points.sf, size = 10, color = "black", shape = "\u2605") +
  theme_map() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 5),
    panel.background = element_rect(fill = "white")
  )


ggsave(mymap,
  file = paste0("outputs/figures/manuscript/regional_map.tiff"),
  width = 5, height = 4.25, units = "in", dpi = 600,
)


# COMBINE INTO ONE MAP W INSET ----
p1 <- ggdraw() +
  draw_image(paste0(
    "outputs/figures/manuscript/regional_map.tiff"
  ), scale = 0.9)
p2 <- ggdraw() + cowplot::draw_image(paste0("outputs/figures/manuscript/map.tiff"), scale = 0.9)


plot.with.inset <-
  ggdraw() +
  draw_plot(p2) +
  draw_plot(p1,
    x = 0.15, y = 0.12, width = 5 / 18, height = 4.25 / 18,
    hjust = 0, vjust = 0
  )

tiff(paste0("outputs/figures/manuscript/inset_map.tiff"),
  width = 5, height = 4.25, units = "in", res = 600
)

plot(plot.with.inset)
dev.off()
