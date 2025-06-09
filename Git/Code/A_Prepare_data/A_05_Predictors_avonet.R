#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#               06_Predictors_avonet.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#

source(here::here("Code/00_Configuration.R"))
lapply(package_list, require, character = TRUE)

tax_path <-
  here("Data/input/Tax_lookup.csv")

#----------------------------------------------------------#
# Load data  -----
#----------------------------------------------------------#

sci_name <-
  readRDS(here::here("Data/output/1_data_filtered.rds")) %>%
  distinct(verbatimIdentification, scientificName) # BirdLife 2024 taxonomy

Tax <- read.csv(tax_path)

traits <-
  read_excel(here::here("Data/input/AVONET Supplementary dataset 1.xlsx"),
              sheet = "AVONET1_BirdLife")  %>%
  dplyr::select(
    Species1, Mass,
    Habitat, Migration, Primary.Lifestyle) %>%
  mutate(
    across(c(Habitat, Migration, Primary.Lifestyle),
           as.factor))

traits_merged <- traits %>%
  rename("ScientificName2018" = "Species1") %>%
  right_join(Tax %>%
               select(verbatimIdentification, scientificName, ScientificName2018)) %>%
  select(-ScientificName2018)

traits_merged %>% skimr::skim()

#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#

saveRDS(traits_merged, here::here("Data/output/A_predictors/Avonet.rds"))

