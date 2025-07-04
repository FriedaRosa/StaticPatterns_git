#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                   11_Predictors_spatial_autocorrelation.R
#                
#
#              Friederike Wölke, Carmen Soria  
#                          2025
#
#----------------------------------------------------------#


# Start with clean environment
rm(list = ls())
gc()

library(here)
#----------------------------------------------------------#
# Install and load libraries -----
#----------------------------------------------------------#

# Source 00_Configuration.R
source(here::here("Code/00_Configuration.R"))

# subset of packages to load to avoid conflicts:
list_p <- c("here", "sf" ,"tidyverse", "purrr", "tictoc", "spdep", "rstatix", "furrr", "futures")
lapply(list_p, require, character = TRUE)


# set variable:
# Define a mapping of datasetID to CRS
dataset_crs_map <- c("5" = "epsg:5514",  # Czechia
                     "6" = "epsg:32118", # NewYork
                     "13" = "epsg:6684", # Japan
                     "26" = "epsg:3035") # Europe

#----------------------------------------------------------#
# Load data -----
#----------------------------------------------------------#

# Load spatial dataset
atlas_sf <- readRDS(here::here("Data/output/1_data_sf.rds")) %>%
  filter(
    scalingID == 1 &
    cell_sampling_repeats == 2 &
    sp_sampling_repeats == 2 &
    species_keep == 1
  )

#----------------------------------------------------------#
# Custom function to process one dataset -----
#----------------------------------------------------------#

library(tidyverse)
library(sf)
library(spdep)
library(purrr)

## Test:
dataset_group <- atlas_sf %>%
  group_by(datasetID) %>%
  group_split() %>% .[[4]]



process_dataset <- function(dataset_group) {
  dataset_id <- unique(dataset_group$datasetID)

  # Transform to correct CRS
  dataset_group <- dataset_group %>% st_as_sf() %>%
    st_transform(crs = dataset_crs_map[[as.character(dataset_id)]])

  # Get all unique sampling periods
  sampling_periods <- unique(dataset_group$samplingPeriodID)

  # Safe versions of spatial autocorrelation functions
  safe_joincount_test <- possibly(function(x, nb) joincount.test(as.factor(x), nb2listw(nb, style = "B", zero.policy = TRUE)), otherwise = NA)

  # Process each sampling period
  res <- map_dfr(sampling_periods, function(sampling_period) {

    # Subset data for the time period
    tp_subset <- dataset_group %>% filter(samplingPeriodID == sampling_period)

    # Get all unique species names
    sp_names <- tp_subset %>%
      st_drop_geometry() %>%
      distinct(verbatimIdentification) %>%
      pull() %>%
      unique()

    # Get all sites
    all_sites <- dataset_group %>% select(siteID, geometry) %>% unique()

    # Process each species
    map_dfr(sp_names, function(sp_name) {

      message(paste0("ID = ", dataset_id, ", TP = ", sampling_period, ", SP = ", sp_name))
      #print(paste0("ID = ", dataset_id, ", TP = ", sampling_period, ", SP = ", sp_name))

      # Filter for a single species and join with all sites
      sp_data <- tp_subset %>%
        filter(verbatimIdentification == sp_name) %>%
        st_drop_geometry() %>%
        full_join(all_sites, by = "siteID") %>%
        mutate(
          presence = if_else(is.na(verbatimIdentification), 0, 1),
          verbatimIdentification = sp_name
        ) %>%
        st_as_sf()

      # Count the number of presences
      num_pres <- sum(sp_data$presence, na.rm = TRUE)

      # Creating neighbors (queen)
      nb_q <- poly2nb(sp_data, queen = TRUE)

      # Join count test
      jc_res_q <- safe_joincount_test(sp_data$presence, nb_q)

      print(unique(paste0("n: ", num_pres,
                   ", joincount: ", unique(ifelse(is.na(jc_res_q), NA, jc_res_q[[2]]$estimate[1])[[1]]))))
      # Save results in a dataframe
      tibble(
        datasetID = dataset_id,
        samplingPeriodID = sampling_period,
        verbatimIdentification = sp_name,
        presence_n = num_pres,
        joincount_statistic = ifelse(is.na(jc_res_q), NA, jc_res_q[[2]]$estimate[1])[[1]],
        joincount_expectation = ifelse(is.na(jc_res_q), NA, jc_res_q[[2]]$estimate[2])[[1]],
        joincount_p_val = ifelse(is.na(jc_res_q), NA, jc_res_q[[2]]$p.value)[[1]],
        joincount_delta = ifelse(is.na(jc_res_q), NA,
                                 ((jc_res_q[[2]]$estimate[1] /
                                     num_pres) -
                                    (jc_res_q[[2]]$estimate[2] /
                                       num_pres)))[[1]]
      ) %>%
        unique()

      }) # map_dfr for species
  }) # map_dfr for sampling periods
gc()
return(res)
}

#----------------------------------------------------------#
# Run the function  -----
#----------------------------------------------------------#
tictoc::tic()
results_CZ <-
  atlas_sf %>%
  filter(datasetID == 5) %>%
  process_dataset()
tictoc::toc()
#----------------------------------------------------------#

tictoc::tic()
results_NY <-
  atlas_sf %>%
  filter(datasetID == 6) %>%
  process_dataset()
tictoc::toc()
#----------------------------------------------------------#

tictoc::tic()
results_JP <-
  atlas_sf %>%
  filter(datasetID == 13) %>%
  process_dataset()
tictoc::toc()
#----------------------------------------------------------#

tictoc::tic()
results_EU <-
  atlas_sf %>%
  filter(datasetID == 26) %>%
  process_dataset()
tictoc::toc()

#----------------------------------------------------------#
# Convert list to a single data frame ----
#----------------------------------------------------------#

results <-
  bind_rows(
  results_EU,
  results_CZ,
  results_JP,
  results_NY)

results %>% skimr::skim()

#----------------------------------------------------------#
# save results to .rds -----
#----------------------------------------------------------#
saveRDS(results, here("Data/output/A_predictors/Spatial_auto.rds"))

