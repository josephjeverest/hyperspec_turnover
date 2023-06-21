# TESTS

# 1 - COVER ----

comp.new <- read.csv("outputs/output_saddle_cover.csv")

comp.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_saddle_cover.csv")

arsenal::comparedf(comp.new, comp.old)


# 2 - TRAITS ----

traits.new <- read.csv("outputs/output_saddle_composition_traits.csv")

traits.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_saddle_composition_traits.csv")

arsenal::comparedf(traits.new, traits.old)


# 3 - BIOMASS ----

biomass.new <- read.csv("outputs/output_biomass.csv")

biomass.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_biomass.csv")

arsenal::comparedf(biomass.new, biomass.old)


# 4 - PCA SPECTRA ----

pca.new <- read.csv("outputs/output_saddle_spectra_b1_PCA.csv")

pca.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_saddle_spectra_b1_PCA.csv")

arsenal::comparedf(pca.new, pca.old)


# 5 - DISTANCE ----

dist.s.new <- read.csv("outputs/output_beta_FULL_spatial.csv")

dist.s.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_beta_FULL_spatial.csv")

arsenal::comparedf(dist.s.new, dist.s.old)


dist.t.new <- read.csv("outputs/output_beta_FULL_temporal.csv")

dist.t.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_beta_FULL_temporal.csv")

arsenal::comparedf(dist.t.new, dist.t.old)


# 7 - MANTEL ----

mantel.s.new <- read.csv("outputs/output_statistics_PCA_spatial_Euclidean.csv")

mantel.s.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_statistics_PCA_spatial_Euclidean.csv")

arsenal::comparedf(mantel.s.new, mantel.s.old)


mantel.t.new <- read.csv("outputs/output_statistics_PCA_temporal_Euclidean.csv")

mantel.t.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_statistics_PCA_temporal_Euclidean.csv")

arsenal::comparedf(mantel.t.new, mantel.t.old)


# 8 - NDVI ----

dist.ndvi.s.new <- read.csv("outputs/output_beta_ndvi_spatial.csv")

dist.ndvi.s.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_beta_ndvi_spatial.csv")

arsenal::comparedf(dist.ndvi.s.new, dist.ndvi.s.old)


dist.ndvi.t.new <- read.csv("outputs/output_beta_ndvi_temporal.csv")

dist.ndvi.t.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_beta_ndvi_temporal.csv")

arsenal::comparedf(dist.ndvi.t.new, dist.ndvi.t.old)


mantel.ndvi.s.new <- read.csv("outputs/output_statistics_spatial_ndvi.csv")

mantel.ndvi.s.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_statistics_spatial_ndvi.csv")

arsenal::comparedf(mantel.ndvi.s.new, mantel.ndvi.s.old)


mantel.ndvi.t.new <- read.csv("outputs/output_statistics_temporal_ndvi.csv")

mantel.ndvi.t.old <- read.csv("C:/Users/s2125726/Documents/CUBoulder/hyperspectral/outputs_remove_no/output_statistics_temporal_ndvi.csv")

arsenal::comparedf(mantel.ndvi.t.new, mantel.ndvi.t.old)
