#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                      13_Merge_predictors.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#

rm(list=ls())
gc()

#----------------------------------------------------------#
# Install and load libraries -----
#----------------------------------------------------------#
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(skimr)
library(sf)

#fix_windows_histograms()
#----------------------------------------------------------#
# Source config  -----
#----------------------------------------------------------#
source(here::here("Code/00_Configuration.R"))


#----------------------------------------------------------#
# Get raw species lists with filtering columns -----
#----------------------------------------------------------#
# Raw species data (including columns to filter just to check if other species were lost)

dta0 <-
  readRDS(here("Data/output/1_data_sf.rds")) %>%
  st_drop_geometry() %>%
  filter(scalingID == 1 & cells_keep == 1) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification, scientificName,
    introduced, sp_remove_expert, sp_sampling_repeats, species_keep) %>%
  filter(!is.na(verbatimIdentification)) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification, species_keep, .keep_all = TRUE) %>%
  mutate(sp_sampling_repeats = case_when(is.na(sp_sampling_repeats) ~ 0,
                                        !is.na(sp_sampling_repeats) ~ sp_sampling_repeats,
                                        .default = sp_sampling_repeats)) %>%
  distinct()


#----------------------------------------------------------#
# Load predictors -----
#----------------------------------------------------------#

## let's match first those with the same (approx) taxonomy

range_size <-
   readRDS(here("Data/output/A_predictors/RangeSizes.rds"))

avonet <-
   readRDS(here("Data/output/A_predictors/Avonet.rds"))

iucn <-
   readRDS(here("Data/output/A_predictors/IUCN_2025_03_25.rds")) %>%
   # manually fix non matches:
   mutate(code = case_when(scientificName %in% c("Morus bassanus",  "Agropsar philippensis", "Tadorna tadorna") ~ "LC",
                           .default = code))

 phylo_dist <-
   readRDS(here("Data/output/A_predictors/Phylo_distinct.rds")) # Phylogenetic distinctiveness

 species_predictors <-
   full_join(range_size, avonet, relationship = "many-to-many") %>%
   full_join(iucn, relationship = "many-to-many") %>%
   full_join(phylo_dist, relationship = "many-to-many") %>%
   distinct()

# Check NAs
colSums(is.na(species_predictors))

# Merge species with traits
species_predictors2 <- dta0 %>%
  filter() %>%
  left_join(species_predictors, relationship = "many-to-many") %>%
  rename("IUCN" = "code")

#-----------------------------------------------------#

## now let's merge those with a similar taxonomy

big_table <-
  readRDS(here("Data/output/A_predictors/Big_table.rds"))

geometry <-
  readRDS(here("Data/output/A_predictors/Range_geometries.rds")) %>% # Species ranges, Atlas geometry
  select(datasetID, samplingPeriodID, verbatimIdentification, circNorm, minDist_toBorder_centr)

sac_metrics <-
  readRDS(here("Data/output/A_predictors/Spatial_auto.rds")) %>%
  select(datasetID, samplingPeriodID, verbatimIdentification,joincount_delta)


# -----------------------------------------------------#
# Handle Lacunarity data
# -----------------------------------------------------#

lacunarity <-
  readRDS(here("Data/output/A_predictors/Lacunarity.rds")) %>%
  ungroup() %>%
  select(-name) %>%
  mutate(samplingPeriodID = as.numeric(as.character(samplingPeriodID)),
         datasetID = as.numeric(as.character(datasetID))) %>%
  mutate(verbatimIdentification = gsub("_", " ", verbatimIdentification)) %>%
  as.data.frame()
names(lacunarity) <- c("r", "ln(r)", "lac", "ln(lac)", "datasetID", "samplingPeriodID", "verbatimIdentification")


# Calculate mean across increasing window sizes
mean_lac <- lacunarity %>%
  group_by(datasetID, samplingPeriodID, verbatimIdentification) %>%
  summarize(mean_lnLac = mean(`ln(lac)`, na.rm = TRUE)) %>%
  mutate(
    verbatimIdentification = case_when(verbatimIdentification == "Fringilla moringilla" ~ "Fringilla montifringilla",
                                       verbatimIdentification == "Moringilla nivalis" ~ "Montifringilla nivalis",
                                       .default = verbatimIdentification)
  )

# quick check on lacunarity data
mean_lac %>%
  group_by(datasetID, samplingPeriodID) %>%
  summarize(n_sp = n_distinct(verbatimIdentification))


# Not matched:
setdiff(mean_lac$verbatimIdentification, dta0$verbatimIdentification) #7
setdiff(dta0$verbatimIdentification,mean_lac$verbatimIdentification) #0

#----------------------------------------------------------#
# Merge predictors -----
#----------------------------------------------------------#

predictors <-
  species_predictors2 %>%
  full_join(big_table, relationship =
              "many-to-many") %>%
  full_join(sac_metrics, relationship =
              "many-to-many") %>%
  full_join(geometry, relationship =
              "many-to-many") %>%
  left_join(mean_lac, relationship =
              "many-to-many") %>%
  distinct(datasetID, verbatimIdentification, samplingPeriodID,
           .keep_all = TRUE) %>%
  mutate(
    across(
      where(is.character) & !matches("verbatimIdentification") & !matches("scientificName"),
      as.factor)) %>%
  mutate(
    across(c("datasetID","samplingPeriodID",
             "Habitat", "IUCN",
             "Migration", "Primary.Lifestyle",
             "introduced", "sp_remove_expert", "sp_sampling_repeats", "species_keep"),
           as.factor)) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification, species_keep,
           .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(Habitat_5 = as.factor(case_when(Habitat %in% c("Grassland", "Shrubland", "Desert", "Rock") ~ "open",
                                         Habitat %in% c("Woodland", "Forest") ~ "closed",
                                         Habitat %in% c("Coastal", "Marine") ~ "marine",
                                         Habitat %in% c("Wetland", "Riverine") ~ "freshwater",
                                         Habitat == "Human Modified" ~ "human",
                                         .default = NA_character_
  ))) %>%
  mutate(Generalism = as.factor(case_when(Primary.Lifestyle == "Generalist" ~ 1,
                                          .default = 0
  ))) %>%
  mutate(Threatened = as.factor(case_when(IUCN %in% c("LC") ~ 0,
                                          is.na(IUCN) ~ NA,
                                          .default = 1
  ))) %>%
  select(-Habitat, -Primary.Lifestyle, -IUCN)

#----------------------------------------------------------#
# Check predictors -----
#----------------------------------------------------------#

predictors %>%
  glimpse()

str(predictors)
predictors %>%
  filter(species_keep == 1 & samplingPeriodID == 1) %>%
  is.na() %>%
  colSums()

predictors %>%
  filter(samplingPeriodID == 1) %>%
  skim() %>%
  to_long() %>%
  setNames(c("variable_class", "variable", "metric", "value"))

predictors %>% filter(is.na(species_keep))

predictors %>%
  filter(samplingPeriodID == 1) %>%
  group_by(datasetID) %>%
  skim() %>%
  as_tibble() %>%
  write.csv(here("Documentation/META_predictors_skim_summary_v2_reduced_version.csv"))

names(predictors$D_AOO_a) <- NULL
names(predictors$morans_I) <- NULL
names(predictors$morans_I_p) <- NULL
names(predictors$Lac) <- NULL
str(predictors)


final_predictors <- predictors %>%
  filter(species_keep == 1) %>%
  select(-sp_remove_expert, -sp_sampling_repeats, -species_keep, -introduced) %>%
  distinct()

saveRDS(final_predictors, here("Data/output/1_all_predictors_merged.rds"))


final_predictors %>% is.na() %>% colSums()

skimr::skim(final_predictors)
