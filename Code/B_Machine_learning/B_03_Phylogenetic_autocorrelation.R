#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#              B_03_Phylogenetic_autocorrelation.R
#                
#
#                    Friederike Wölke 
#                           2025
#
#----------------------------------------------------------#

source(here::here("Code/00_Configuration.R"))
package_list <- 
  c("here","readr","ggplot2","dplyr", "phylosignal", "ape", "phylobase", "tidymodels", "ranger", "purrr", "tibble", "ggplot2")
x <- lapply(package_list, require, character = TRUE)
rm(list = ls())


# Data paths:
tax_path <-
  here("Data/input/Tax_lookup.csv")

dta_path <- 
  here("Data/output/1_all_predictors_merged.rds")

tree_path <- 
  here("Data/input/MRC_consensus_BirdTree.tre")

ranger_res_pooled_path <-
  here("Data/output/B_models/B_01_list_all_results_rf.rds")

ranger_res_regional_path <- 
  here("Data/output/B_models/B_02_rf_res_atlas_split.rds")

#----------------------------------------------------------#
# Load data  -----
#----------------------------------------------------------#

dta <- 
  readRDS(dta_path) %>%
  distinct(scientificName, verbatimIdentification, 
           datasetID, 
           Jaccard_dissim, log_R2_1, log_R2_1_per_year )

Tax <-
  read_csv(tax_path) %>%
  select(verbatimIdentification, scientificName, BirdTree) %>%
  full_join(dta) %>%
  filter(!c(is.na(datasetID))) %>%
  mutate(BirdTree = gsub(" ", "_", BirdTree))

tree <-
  ape::read.tree(tree_path)

# BirdTree in data per atlas
sp_in_dta <-
  Tax %>%
  group_by(datasetID) %>%
  group_split()

res_full <- 
  readRDS(ranger_res_pooled_path)
res_split <- 
  readRDS(ranger_res_regional_path)

# ──────────────────────────────────────────────────────────────#
# I. Phylogenetic autocorrelation in the raw data
# --------------------------------------------------------------#

# test subset:
# dat <- sp_in_dta[[1]] %>% head()
# tree <- tree
# region <- "czechia"
# response <- "Jaccard_dissim"


## Helper function for raw data ##
phylo_corr_df <- function(dat, tree, region) {
  
  ## 1 ─ mean value of each trait per tip
  dd_new <- dat %>%
    filter(!is.na(BirdTree)) %>%
    group_by(BirdTree) %>%
    summarise(across(
      c(Jaccard_dissim, log_R2_1, log_R2_1_per_year),
      ~ mean(.x, na.rm = TRUE)
    )) %>%
    column_to_rownames("BirdTree") %>%
    na.omit()
  
  ## 2 ─ prune tree to those tips
  tre_pruned <- tree %>%
    drop.tip(setdiff(tree$tip.label, rownames(dd_new))) %>%
    ladderize()
  
  ## 3 ─ build phylo4d object once
  p4d <- phylo4d(tre_pruned, dd_new)
  
  ## 4 ─ helper to reshape one correlogram
  corr_to_tbl <- function(response) {
    phyloCorrelogram(p4d, trait = response, ci.bs = 100)$res %>%
      as_tibble() %>%
      mutate(
        dist      = .[[1]], # mid‑point of the distance class
        I         = .[[4]], # Moran’s I 
        lower_CI  = .[[2]], # lower CI
        upper_CI  = .[[3]], # upper CI
        model     = "none",
        type      = "raw data",
        response  = response,
        region    = region,
        N         = nrow(dd_new),
        .keep = "used"
      )
  }
  
  ## 5 ─ run for the three traits & bind rows 
  map_dfr(
    #c("Jaccard_dissim"),
    c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year"),
    corr_to_tbl
  )
}

# ---------------------------------------------------- #
# Run the function for raw data phylo autocorr 
# ---------------------------------------------------- #

# Across regions:
region_names <- c("Czechia", "New York", "Japan", "Europe")

all_corr <- map2_dfr(
  sp_in_dta,                        # list of four data frames
  region_names,                     # parallel vector of region labels
  ~ phylo_corr_df(.x, tree, .y)     # .x = data, .y = region
)

# Check results
head(all_corr)
saveRDS(all_corr, here("Data/output/temp/B_03_phylo_autocorr_raw_data.rds"))
all_corr <- readRDS(here("Data/output/temp/B_03_phylo_autocorr_raw_data.rds"))
# Plot results
ggplot(all_corr %>% filter(response != "log_R2_1_per_year"),
       aes(dist, I, colour = response, fill = response)) +
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI),
              alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ region) +
  labs(x = "Phylogenetic distance", y = "Moran’s I") +
  theme_classic()



# ──────────────────────────────────────────────────────────────#
# II. Phylogenetic autocorrelation in the model residuals
# --------------------------------------------------------------#

## Helper function for residuals ##

resid_correlogram_df <- function(pred_df,
                                 tree,
                                 Tax,
                                 response = c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year"), 
                                 model   = c("pooled", "regional"),
                                 region  = c("all", "cz", "ny", "jp", "eu")) {
  
  model  <- match.arg(model)
  region <- match.arg(region)
  
  ## 1 ─ mean residual per tip -------------------------------------------
  tip_resid <- pred_df %>%
    select(verbatimIdentification, resid) %>%
    left_join(
      Tax %>% select(verbatimIdentification, BirdTree),
      by = "verbatimIdentification"
    ) %>%
    distinct() %>%
    filter(!is.na(BirdTree)) %>%
    group_by(BirdTree) %>%
    summarise(resid = mean(resid, na.rm = TRUE), .groups = "drop") %>%
    column_to_rownames("BirdTree")
  
  ## 2 ─ prune the tree ---------------------------------------------------
  tre_pruned <- tree %>%
    drop.tip(setdiff(tree$tip.label, rownames(tip_resid))) %>%
    ladderize()
  
  ## 3 ─ correlogram ------------------------------------------------------
  cor_obj <- phylo4d(tre_pruned, tip_resid) %>%
    phyloCorrelogram(trait = "resid", ci.bs = 100)
  
  ## 4 ─ tidy output ------------------------------------------------------
  as_tibble(cor_obj$res) %>%
    mutate(
      dist      = .[[1]], # mid‑point of the distance class
      I         = .[[4]], # Moran’s I 
      lower_CI  = .[[2]], # lower CI
      upper_CI  = .[[3]], # upper CI
      model     = model,
      type      = "residuals",
      response  = response,
      region    = region,
      N         = nrow(tip_resid),
      .keep = "used"
    )
}



# ---------------------------------------------------- #
# Run the function for model residuals phylo autocorr 
# ---------------------------------------------------- #

# a) Pooled model
responses <- c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year")

corr_pooled <- purrr::imap_dfr(
  res_full$predictions,                           # list length = 3
  ~ resid_correlogram_df(
    pred_df  = .x,
    tree     = tree,
    Tax   = Tax,
    response = responses[.y],                   # tag by position
    model    = "pooled",                        
    region = "all"
  )
)

# b) models split by region
regions   <- c("cz", "ny", "jp", "eu")           
responses <- c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year")

corr_regional <- purrr::imap_dfr(res_split, \(region_list, i_reg) {
  region <- regions[i_reg]
  
  purrr::imap_dfr(region_list, \(resp_list, i_resp) {
    resid_correlogram_df(
      pred_df  = resp_list$predictions,
      tree     = tree,
      Tax   = Tax,
      response = responses[i_resp],
      model    = "regional",
      region   = region
    )
  })
})


# Merge split and pooled results
corr_residuals_all <- bind_rows(corr_pooled, corr_regional) 
saveRDS(corr_residuals_all, here("Data/output/temp/B_03_phylo_autocorr_model_residuals.rds"))

# Check results
head(corr_residuals_all)

# Plot results
corr_residuals_all %>%
  mutate(region = case_when(
    region == "all" ~ "all",
    region == "cz" ~ "Czechia",
    region == "ny" ~ "New York",
    region == "jp" ~ "Japan",
    region == "eu" ~ "Europe"
  )) %>%
  ggplot(aes(dist, I, colour = response, fill = response)) +
  geom_hline(yintercept = 0) +
  # geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI),
  #             alpha = 0.25, colour = NA) +
  geom_line(size = 1.1) +
  facet_wrap(~ region + model) +
  ylim(-0.4, 0.43) +
  labs(x = "Phylogenetic distance", y = "Moran’s I of residuals") +
  ggthemes::theme_few() +
  geom_line(
    data = all_corr %>% mutate(model = "regional"), 
    aes(dist, I, linetype = response)
  )
  




ggplot(all_corr %>% filter(response != "log_R2_1_per_year"),
       aes(dist, I, colour = response, fill = response)) +
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI),
              alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ region) +
  labs(x = "Phylogenetic distance", y = "Moran’s I") +
  theme_classic()