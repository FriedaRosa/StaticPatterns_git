#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                    D_05_Figure_S12a
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#

# Libraries
library(hstats)
library(here)
library(dplyr)
library(tidymodels)
library(ranger)
library(ggplot2)
library(ggthemes)

# Start with clean environment
rm(list = ls())
gc()

#----------------------------------------------------------#


## Interactions & their partial dependencies ##

data <- 
  readRDS(here::here("Data/output/B_models/B_01_list_all_results_rf.rds"))

fit_1 <- 
  data$res[[1]] %>% 
  extract_fit_parsnip()

train <- 
  data$list_split_data[[1]]$train

pdp_J_interactions <- 
  list()
#----------------------------------------------------------#


# All 20 interactions from Jaccard: #

#----------------------------------------------------------#
# Pairwise Interactions ----
#----------------------------------------------------------#

# 1
pdp_J_interactions[[1]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("mean_lnLac", "joincount_delta"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("lac ~ jc ~ dataset")+
  ggthemes::theme_few()

# 2
pdp_J_interactions[[2]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("AOO", "joincount_delta"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("AOO ~ jc ~ dataset")+
  ggthemes::theme_few()

# 3
pdp_J_interactions[[3]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("D_AOO_a", "AOO"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("D ~ AOO ~ dataset")+
  ggthemes::theme_few()

# 4
pdp_J_interactions[[4]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("D_AOO_a", "joincount_delta"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("D ~ jc ~ dataset")+
  ggthemes::theme_few()

# 5
pdp_J_interactions[[5]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("D_AOO_a", "circNorm"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("D ~ circ ~ dataset")+
  ggthemes::theme_few()

# 6
pdp_J_interactions[[6]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("mean_lnLac", "AOO"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("lac ~ AOO ~ dataset")+
  ggthemes::theme_few()

# 7
pdp_J_interactions[[7]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("AOO", "circNorm"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("AOO ~ circ ~ dataset")+
  ggthemes::theme_few()

# 8
pdp_J_interactions[[8]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("joincount_delta", "circNorm"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("jc ~ circ ~ dataset")+
  ggthemes::theme_few()

# 9
pdp_J_interactions[[9]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("mean_lnLac", "circNorm"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("lac ~ circ ~ dataset")+
  ggthemes::theme_few()

# 10
pdp_J_interactions[[10]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("D_AOO_a", "mean_lnLac"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("D ~ lac ~ dataset")+
  ggthemes::theme_few()

#----------------------------------------------------------#
# Threeway Interactions ----
#----------------------------------------------------------#

# 11
pdp_J_interactions[[11]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("D_AOO_a", "AOO"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("D ~ AOO ~ circ")+
  ggthemes::theme_few()

# 12
pdp_J_interactions[[12]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("D_AOO_a", "mean_lnLac"),
    BY = "joincount_delta",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("D ~ lac ~ jc")+
  ggthemes::theme_few()

# 13
pdp_J_interactions[[13]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("D_AOO_a", "AOO"),
    BY = "joincount_delta",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("D ~ AOO ~ jc")+
  ggthemes::theme_few()

# 14
pdp_J_interactions[[14]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("mean_lnLac", "AOO"),
    BY = "joincount_delta",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("lac ~ AOO ~ jc")+
  ggthemes::theme_few()

# 15
pdp_J_interactions[[15]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("AOO", "joincount_delta"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("AOO ~ jc ~ circ")+
  ggthemes::theme_few()

# 16
pdp_J_interactions[[16]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("D_AOO_a", "mean_lnLac"),
    BY = "AOO",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("D ~ lac ~ AOO")+
  ggthemes::theme_few()

# 17
pdp_J_interactions[[17]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("D_AOO_a", "joincount_delta"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("D ~ jc ~ circ")+
  ggthemes::theme_few()

# 18
pdp_J_interactions[[18]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("mean_lnLac", "joincount_delta"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("lac ~ jc ~ circ")+
  ggthemes::theme_few()

# 19
pdp_J_interactions[[19]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("mean_lnLac", "AOO"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("lac ~ AOO ~ circ")+
  ggthemes::theme_few()

# 20
pdp_J_interactions[[20]] <- 
  partial_dep(
    object = fit_1,
    X = train,
    v = c("D_AOO_a", "mean_lnLac"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("D ~ lac ~ circ")+
  ggthemes::theme_few()


saveRDS(pdp_J_interactions, here("Data/output/temp/D_05_J_pdp_all.rds"))
pdp_J_interactions <- readRDS(here("Data/output/temp/D_05_J_pdp_all.rds"))


lapply(pdp_J_interactions, plot)
