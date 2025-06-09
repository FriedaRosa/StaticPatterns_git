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
  c(package_list, "phylosignal", "ape", "phylobase", "tidymodels", "ranger")
lapply(package_list, require, character = TRUE)


tax_path <-
  here("Data/input/Tax_lookup.csv")

dta_path <- here("Data/output/1_all_predictors_merged.rds")

#----------------------------------------------------------#
# Load data  -----
#----------------------------------------------------------#

dta <- readRDS(dta_path) %>%
  distinct(scientificName, verbatimIdentification, datasetID, Jaccard_dissim, log_R2_1, log_R2_1_per_year )

Tax <-
  read_csv(tax_path) %>%
  select(verbatimIdentification, scientificName, BirdTree) %>%
  full_join(dta) %>%
  filter(!c(is.na(datasetID))) %>%
  mutate(BirdTree = gsub(" ", "_", BirdTree))

tree <-
  ape::read.tree(
    here::here("Data/input/MRC_consensus_BirdTree.tre"))

# BirdTree in data per atlas
sp_in_dta <-
  Tax %>%
  group_by(datasetID) %>%
  group_split()

#----------------------------------------------------------#
# Calculate phylogenetic autocorrelation in the model residuals ----
#----------------------------------------------------------#
tictoc::tic()
res_full <- readRDS(here::here("Data/output/B_models/B_01_list_all_results_rf.rds"))
res_split <- readRDS(here::here("Data/output/B_models/B_02_rf_res_atlas_split.rds"))
Tax

# Overall
cor_list <- replicate(3, list)
for (resp_i in seq_along(res_full$predictions)){
  cor_list[[resp_i]] <- list()
  predictions_i <-
    res_full$predictions[[resp_i]]

  plot(predictions_i[[2]], predictions_i$resid)

  sp_residuals <-
    predictions_i %>%
    select(verbatimIdentification, resid) %>%
    left_join(Tax %>% select(verbatimIdentification, BirdTree))

  tre <- drop.tip(tree, setdiff(tree$tip.label, sp_residuals$BirdTree)) %>%
    ladderize()

  dta_i <- sp_residuals %>%
    aggregate(resid~BirdTree, FUN = mean)

  p4d_i <- phylo4d(tre, dta_i)
  cor_i <- phyloCorrelogram(p4d_i, trait = "resid")
  cor_list[[resp_i]] <- list(cor_i, p4d_i)
}

saveRDS(cor_list, here::here("Data/output/temp/B_03_p4d_residuals_all_data.rds"))
tictoc::toc()


# Per atlas:
tictoc::tic()
cor_list <- replicate(3, list)
for (resp_i in seq_along(res_full$predictions)){
  cor_list[[resp_i]] <- list()
  predictions_i <-
    res_full$predictions[[resp_i]]
  cor_list[[resp_i]] <- list()
  atlas_split_list <- predictions_i %>%
    group_by(datasetID) %>%
    group_split()

  for(atlas_i in seq_along(atlas_split_list)){
    this_atlas <- atlas_split_list[[atlas_i]]
    sp_residuals <-
      this_atlas %>%
      select(verbatimIdentification, resid) %>%
      left_join(Tax %>% select(verbatimIdentification, BirdTree))

    tre <- drop.tip(tree, setdiff(tree$tip.label, sp_residuals$BirdTree)) %>%
      ladderize()

    dta_i <- sp_residuals %>%
      aggregate(resid~BirdTree, FUN = mean)

    p4d_i <- phylo4d(tre, dta_i)
    cor_i <- phyloCorrelogram(p4d_i, trait = "resid")
    cor_list[[resp_i]][[atlas_i]] <- list(p4d_i, cor_i)

    }

  }

saveRDS(cor_list, here::here("Data/output/temp/B_03_p4d_residuals_full_model_split_by_atlas.rds"))
tictoc::toc()

# Save Plots -----

cor_list <- readRDS(here::here("Data/output/temp/B_03_p4d_residuals_full_model_split_by_atlas.rds"))


svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_Jaccard_CZ.svg", width = 6, height = 4)
cor_list[[1]][[1]][[2]] %>% plot(main = "Jaccard - Czechia") 
dev.off()

svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_Jaccard_NY.svg",  width = 6, height = 4)
cor_list[[1]][[2]][[2]] %>% plot(main = "Jaccard - New York State")
dev.off()

svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_Jaccard_JP.svg",  width = 6, height = 4)
cor_list[[1]][[3]][[2]] %>% plot(main = "Jaccard - Japan")
dev.off()

svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_Jaccard_EU.svg",  width = 6, height = 4)
cor_list[[1]][[4]][[2]] %>% plot(main = "Jaccard - Europe")
dev.off()



svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_logRatio_CZ.svg", width = 6, height = 4)
cor_list[[2]][[1]][[2]] %>% plot(main = "log ratio - Czechia")
dev.off()

svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_logRatio_NY.svg", width = 6, height = 4)
cor_list[[2]][[2]][[2]] %>% plot(main = "log ratio - New York State")
dev.off()

svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_logRatio_JP.svg", width = 6, height = 4)
cor_list[[2]][[3]][[2]] %>% plot(main = "log ratio - Japan")
dev.off()

svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_logRatio_EU.svg", width = 6, height = 4)
cor_list[[2]][[4]][[2]] %>% plot(main = "log ratio - Europe")
dev.off()



svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_logRatio_y_CZ.svg", width = 6, height = 4)
cor_list[[3]][[1]][[2]] %>% plot(main = "log ratio per year - Czechia")
dev.off()

svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_logRatio_y_NY.svg", width = 6, height = 4)
cor_list[[3]][[2]][[2]] %>% plot(main = "log ratio per year - New York State")
dev.off()

svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_logRatio_y_JP.svg", width = 6, height = 4)
cor_list[[3]][[3]][[2]] %>% plot(main = "log ratio per year - Japan")
dev.off()

svg(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_logRatio_y_EU.svg", width = 6, height = 4)
cor_list[[3]][[4]][[2]] %>% plot(main = "log ratio per year - Europe")
dev.off()




# pdf(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_Jaccard_plots.pdf", paper = "a4")
# cor_list[[1]][[1]][[2]] %>% plot() 
# cor_list[[1]][[2]][[2]] %>% plot()
# cor_list[[1]][[3]][[2]] %>% plot()
# cor_list[[1]][[4]][[2]] %>% plot()
# dev.off()
# 
# pdf(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_log_R2_1_plots.pdf", paper = "a4")
# cor_list[[2]][[1]][[2]] %>% plot()
# cor_list[[2]][[2]][[2]] %>% plot()
# cor_list[[2]][[3]][[2]] %>% plot()
# cor_list[[2]][[4]][[2]] %>% plot()
# dev.off()
# 
# pdf(file = "Figures/B_models/phylo_autocorr/B_03_phylo_autocorr_log_R2_1_y_plots.pdf", paper = "a4")
# cor_list[[3]][[1]][[2]] %>% plot()
# cor_list[[3]][[2]][[2]] %>% plot()
# cor_list[[3]][[3]][[2]] %>% plot()
# cor_list[[3]][[4]][[2]] %>% plot()
# dev.off()





#----------------------------------------------------------#
# Calculate phylogenetic autocorrelation in the raw data----
#----------------------------------------------------------#
# Based on: https://www.francoiskeck.fr/phylosignal/demo_general.html
tictoc::tic()
res <- list()

dd_i <- 1
for (dd_i in seq_along(sp_in_dta)){
  res[[dd_i]] <- list()
  
  tre <- drop.tip(tree, setdiff(tree$tip.label, sp_in_dta[[dd_i]]$BirdTree)) %>%
    ladderize()
  
  dd_new <- aggregate(cbind(Jaccard_dissim, log_R2_1, log_R2_1_per_year) ~ BirdTree, data = sp_in_dta[[dd_i]], FUN = mean)
  
  rownames(dd_new) <- dd_new$BirdTree
  dd_new$BirdTree <- NULL
  dd_new$datasetID <- NULL
  
  p4d <- phylo4d(tre, dd_new)
  
  sig <- phyloSignal(p4d = p4d, method = "K")
  #cor_LR <- phyloCorrelogram(p4d, trait = "log_R2_1")
  #cor_J <- phyloCorrelogram(p4d, trait = "Jaccard_dissim")
  #cor_LR_y <- phyloCorrelogram(p4d, trait = "log_R2_1_per_year")
  
  res[[dd_i]] <- sig %>% bind_cols() %>% setNames(c("K", "p"))
  #res[[dd_i]] <- list(sig, cor_LR, cor_J, cor_LR_y)
}
tictoc::toc()

# Visualize
#plot(res[[1]][[1]], stacked.methods = FALSE, quantiles = c(0.05, 0.95))

res




# Visualize
#plot(res[[1]][[1]], stacked.methods = FALSE, quantiles = c(0.05, 0.95))

res


