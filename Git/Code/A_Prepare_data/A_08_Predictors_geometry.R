#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                   10_Predictors_geometry.R
#                
#
#              Friederike Wölke, Gabriel Ortega-Solís  
#                          2025
#
#----------------------------------------------------------#


# Start with clean environment
rm(list = ls())
gc()


#----------------------------------------------------------#
# Install and load libraries -----
#----------------------------------------------------------#

# Source 00_Configuration.R
source(here::here("Code/00_Configuration.R"))

# subset of packages to load to avoid conflicts:
list_p <- c("here", "sf" ,"tidyverse", "purrr", "tictoc", "terra", "broom", "geosphere")
lapply(list_p, require, character = TRUE)


# Source poly-attribute function
source(here("Code/Functions/poly_attr.R"))

# run function and skip errors:
poly_attr <- possibly(poly_attr, otherwise = NA_real_, quiet = FALSE)
#----------------------------------------------------------#

#----------------------------------------------------------#

sf_use_s2(FALSE) # we use planar projections so this can be turned off

#----------------------------------------------------------#
# Load data -----
#----------------------------------------------------------#

atlas_sf <-
  readRDS("Data/output/1_data_sf.rds") %>%
  filter(scalingID == 1 & cell_sampling_repeats == 2) %>%
  group_split(datasetID)

#----------------------------------------------------------#
# Atlas geometries  -----
#----------------------------------------------------------#

# Define function for processing a single atlas
process_atlas <- function(atlas, crs_value) {

  # Filter and transform dataset
  sf_current_atlas <- atlas %>%
    group_by(datasetID) %>%
    select(datasetID, geometry) %>%
    st_as_sf() %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    st_transform(crs = crs_value) %>%
    st_make_valid()

  atlasID <- unique(sf_current_atlas$datasetID)

  message("Processing atlas ID = ", atlasID)

  # Convert to terra SpatVector
  terra_current_atlas <- vect(sf_current_atlas)

  # Store border lines
  atlas_border_lines <- terra::as.lines(terra_current_atlas)

  # Bounding box metrics
  bbox <- terra::ext(terra_current_atlas)

  # Calculate predictors
  data_atlas_res <- data.frame(
    datasetID = atlasID,
    atlas_xmin = bbox[1],
    atlas_xmax = bbox[2],
    atlas_xhalf = bbox[1] + (bbox[2] - bbox[1]) / 2,
    atlas_ymin = bbox[3],
    atlas_ymax = bbox[4],
    atlas_yhalf = bbox[3] + (bbox[4] - bbox[3]) / 2,
    #atlas_nsDist = poly_attr(terra_current_atlas, "nsDist"),
    #atlas_ewDist = poly_attr(terra_current_atlas, "ewDist"),
    #atlas_maxDist = poly_attr(terra_current_atlas, "maxDist"),
    #atlas_lengthMinRect = poly_attr(terra_current_atlas, "lengthMinRect"),
    #atlas_widthMinRect = poly_attr(terra_current_atlas, "widthMinRect"),
    #atlas_elonMinRect = poly_attr(terra_current_atlas, "elonMinRect"),
    #atlas_elonRatio = poly_attr(terra_current_atlas, "elonRatio"),
    #atlas_circ = poly_attr(terra_current_atlas, "circ"),
    atlas_circNorm = poly_attr(terra_current_atlas, "circNorm")
    #atlas_relCirc = poly_attr(terra_current_atlas, "relCirc"),
    #atlas_lin = poly_attr(terra_current_atlas, "lin"),
    #atlas_bearingMinRect = poly_attr(terra_current_atlas, "bearingMinRect"),
    #atlas_bearing = poly_attr(terra_current_atlas, "bearing"),
    #atlas_perimeter = terra_current_atlas %>% project("epsg:4326") %>% perim()
  )

  return(
    list(
      geometry_metrics = data_atlas_res,
      borders = atlas_border_lines))
}

#----------------------------------------------------------#

# Apply function to each atlas using map2()
results <- map2(atlas_sf, vars$crs, process_atlas)

#----------------------------------------------------------#

# Extract results into separate lists
data_atlas_geometry <- map_dfr(results, "geometry_metrics")

# Combine results into a single data frame
row.names(data_atlas_geometry) <- NULL

# Extract borders
list_borders <- map(results, "borders")


#----------------------------------------------------------#
# Single species range geometries  -----
#----------------------------------------------------------#

list_range_geometries <-
  replicate(4, list())

#----------------------------------------------------------#

tictoc::tic()
for (atlas_i in seq_along(data_atlas_geometry$datasetID)){

  atlas <- data_atlas_geometry[[atlas_i]]
  border_lines <- list_borders[[atlas_i]]
  atlas_centroid <- centroids(border_lines)

  grouped_list <- atlas_sf[[atlas_i]] %>%
    st_as_sf() %>%
    group_by(samplingPeriodID, verbatimIdentification) %>%
    group_split()

  #----------------------------------------#
    range_geometries <- map_dfr(grouped_list, function(dta) {

    #message("Processing atlas nr = ", atlas_i)
    print(paste0("processing = ", unique(dta$verbatimIdentification),
                 ", datasetID = ", unique(dta$datasetID),
                 ", tp = ", unique(dta$samplingPeriodID) ))


    sp_range <- dta %>%
      dplyr::select(geometry, siteID, datasetID) %>%
      st_transform(crs = st_crs(border_lines)) %>%
      summarise() %>%
      terra::vect()

    sp_range_centroid <- centroids(sp_range)

    #----------------------------------------#

    # Compute range attributes
    tibble(
      datasetID = unique(dta$datasetID),
      samplingPeriodID = unique(dta$samplingPeriodID),
      verbatimIdentification = unique(dta$verbatimIdentification),
      circNorm = poly_attr(sp_range, "circNorm"), # based on area-perimeter ratio (Cn=P^2/4πA) --> scale independent ( = 1 if it is a circle,  > 1 if it deviates from a circle)
      minDist_toBorder_centr = min(distance(crds(sp_range_centroid), crds(border_lines), lonlat = F))
    )
  })


  list_range_geometries[[atlas_i]] <- range_geometries

}
tictoc::toc() # 297.64 sec elapsed

#----------------------------------------------------------#

# Bind back together:
data_range_geometries <- list_range_geometries %>%
  bind_rows()

#----------------------------------------------------------#

data_final <- data_range_geometries %>%
  left_join(data_atlas_geometry, by = "datasetID") %>%
  ungroup() %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification, .keep_all = TRUE)

#----------------------------------------------------------#

colSums(is.na(data_final)) #no NAs

#----------------------------------------------------------#
# Save data to .rds  -----
#----------------------------------------------------------#


saveRDS(data_final, here("Data/output/A_predictors/Range_geometries.rds"))
