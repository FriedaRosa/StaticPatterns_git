#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#               08_Predictors_IUCN.R
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

Tax <-
  read_csv(tax_path) %>%
  distinct(verbatimIdentification, scientificName)

#----------------------------------------------------------#
# Get IUCN status  -----
#----------------------------------------------------------#


IUCN.list <-
  iucn_summary(unique(Tax$scientificName),
               distr.detail = F,
               key = Sys.getenv('IUCN_REDLIST_KEY'))

# Extract codes
IUCN <- lapply(names(IUCN.list), function(name) {
  item <- IUCN.list[[name]]

  # Check if the element is a list and contains `red_list_category`
  if (is.list(item) && !is.null(item$red_list_category) && !is.null(item$red_list_category$code)) {
    return(data.frame(name = name, code = item$red_list_category$code, stringsAsFactors = FALSE))
  } else {
    return(data.frame(name = name, code = NA, stringsAsFactors = FALSE)) # Fill missing with NA
  }
}) %>% do.call(rbind, .)

# Check NA species
species_with_NA <-
  IUCN %>% filter(is.na(code))


IUCN.list2 <-
  iucn_summary(unique(species_with_NA$name),
               distr.detail = F,
               key = Sys.getenv('IUCN_REDLIST_KEY'))

IUCN2 <- lapply(names(IUCN.list2), function(name) {
  item <- IUCN.list2[[name]]

  # Check if the element is a list and contains `red_list_category`
  if (is.list(item) && !is.null(item$red_list_category) && !is.null(item$red_list_category$code)) {
    return(data.frame(name = name, code = item$red_list_category$code, stringsAsFactors = FALSE))
  } else {
    return(data.frame(name = name, code = NA, stringsAsFactors = FALSE)) # Fill missing with NA
  }
}) %>% do.call(rbind, .)

IUCN_merged <- rbind(IUCN %>% na.omit(), IUCN2) %>%
  unique() %>%
  mutate(code = case_when(name == "Nannopterum auritus" ~ "LC",
                          .default = code))


IUCN_merged %>% filter(is.na(code))

# Reshape
IUCN_df <-
  as.data.frame(IUCN_merged,
                row.names = c(IUCN_merged$name)) %>%
  rownames_to_column(var = "scientificName") %>%
  right_join(Tax) %>%
  distinct(verbatimIdentification, scientificName, code)

#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#

saveRDS(IUCN_df,
        here::here("Data/output/A_predictors/IUCN_2025_03_25.rds"))
