#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#               09_Predictors_climate_niches.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#


source(here::here("Code/00_Configuration.R"))
lapply(package_list, require, character = TRUE)

#if(!"rasterSp" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/rasterSp", build_vignettes = T)
#if(!"climateNiche" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/climateNiche", build_vignettes = T)

library(rasterSp)
library(climateNiche)
library(factoextra)

tax_path <-
  here("Data/input/Tax_lookup.csv")

# planar CRS
Mollweide_CRS <- 'PROJCS["ProjWiz_Custom_Mollweide",
 GEOGCS["GCS_WGS_1984",
  DATUM["D_WGS_1984",
   SPHEROID["WGS_1984",6378137.0,298.257223563]],
  PRIMEM["Greenwich",0.0],
  UNIT["Degree",0.0174532925199433]],
 PROJECTION["Mollweide"],
 PARAMETER["False_Easting",0.0],
 PARAMETER["False_Northing",0.0],
 PARAMETER["Central_Meridian",0],
 UNIT["Meter",1.0]]'

sf_use_s2(TRUE)

#----------------------------------------------------------#
# Load data  -----
#----------------------------------------------------------#

# Global ranges
BirdLife <- st_read("Data/input/shp_global/BirdLife.shp")

# Taxonomic match table
Tax <-
  read.csv(tax_path) %>%
  select(verbatimIdentification, scientificName)


#----------------------------------------------------------#
# Calculate global range size    -----
#----------------------------------------------------------#

BirdLife_list <- BirdLife %>%
  group_by(sci_name) %>%   # Group by species
  group_split()


#--------------------------------------------------#

# Set up parallel session:
library(furrr)
library(purrr)
#--------------------------------------------------#

plan(multisession, workers = 4)  # Set up parallel processing
sf_use_s2(FALSE)
RangeSizes_l <- furrr::future_map_dfr(BirdLife_list, function(sp_data) {
  sp_data %>%
    group_by(sci_name) %>%
    st_transform(crs = Mollweide_CRS) %>%
    st_make_valid() %>%
    summarise() %>%
    mutate(GlobRangeSize_km2 = as.numeric(st_area(.))/ 1e6) %>%  # Convert to km²
    st_drop_geometry() %>%
    group_by(sci_name) %>%
    summarise(
      GlobRangeSize_km2 = sum(GlobRangeSize_km2, na.rm = TRUE),  # Sum range sizes
      .groups = "keep"
    )
}, .options = furrr_options(seed = TRUE))  # Ensure reproducibility

# reset options
sf_use_s2(TRUE)
plan(sequential)

# Combine lists
RangeSizes <- RangeSizes_l %>%
  rename(scientificName = sci_name) %>%
  right_join(Tax, by = "scientificName") %>%
  distinct(scientificName, verbatimIdentification, .keep_all = TRUE) %>%
  mutate(GlobRangeSize_km2 = ifelse(GlobRangeSize_km2==0, NA, GlobRangeSize_km2)) %>%
  ungroup()


colSums(is.na(RangeSizes))
RangeSizes %>% filter(is.na(GlobRangeSize_km2))
#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#


skim(RangeSizes)
saveRDS(RangeSizes, here::here("Data/output/A_predictors/RangeSizes.rds"))
