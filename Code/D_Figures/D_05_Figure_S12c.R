#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                    D_05_Figure_S12c
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


fit_3 <- 
  data$res[[3]] %>% 
  extract_fit_parsnip()


train <- 
  data$list_split_data[[3]]$train

pdp_lnR_y_interactions <- 
  list()
#----------------------------------------------------------#


# All 20 interactions from Jaccard: #

#----------------------------------------------------------#
# Pairwise Interactions ----
#----------------------------------------------------------#

# 1
pdp_lnR_y_interactions[[1]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("Mass", "AOO"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("Mass ~ AOO ~ dataset")+
  ggthemes::theme_few()

# 2
pdp_lnR_y_interactions[[2]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("D_AOO_a", "AOO"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("D ~ AOO ~ dataset")+
  ggthemes::theme_few()

# 3
pdp_lnR_y_interactions[[3]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("Mass", "joincount_delta"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("Mass ~ jc ~ dataset")+
  ggthemes::theme_few()

# 4
pdp_lnR_y_interactions[[4]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("D_AOO_a", "mean_lnLac"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("D ~ lac ~ dataset")+
  ggthemes::theme_few()

# 5
pdp_lnR_y_interactions[[5]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("AOO", "circNorm"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("AOO ~ circNorm ~ dataset")+
  ggthemes::theme_few()

# 6
pdp_lnR_y_interactions[[6]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("Mass", "D_AOO_a"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("D ~ D ~ dataset")+
  ggthemes::theme_few()

# 7
pdp_lnR_y_interactions[[7]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("joincount_delta", "circNorm"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("jc ~ circ ~ dataset")+
  ggthemes::theme_few()

# 8
pdp_lnR_y_interactions[[8]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("D_AOO_a", "circNorm"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("D ~ circ ~ dataset")+
  ggthemes::theme_few()

# 9
pdp_lnR_y_interactions[[9]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("Mass", "circNorm"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("Mass ~ circ ~ dataset")+
  ggthemes::theme_few()

# 10
pdp_lnR_y_interactions[[10]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("AOO", "joincount_delta"),
    BY = "datasetID",
    grid_size = 100L
  ) %>%
  plot() +
  ggtitle("AOO ~ jc ~ dataset")+
  ggthemes::theme_few()

#----------------------------------------------------------#
# Threeway Interactions ----
#----------------------------------------------------------#

# 11
pdp_lnR_y_interactions[[11]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("Mass", "D_AOO_a"),
    BY = "AOO",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("Mass ~ D ~ AOO")+
  ggthemes::theme_few()

# 12
pdp_lnR_y_interactions[[12]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("D_AOO_a", "AOO"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("D ~ AOO ~ circ")+
  ggthemes::theme_few()

# 13
pdp_lnR_y_interactions[[13]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("Mass", "D_AOO_a"),
    BY = "joincount_delta",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("Mass ~ D ~ jc")+
  ggthemes::theme_few()

# 14
pdp_lnR_y_interactions[[14]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("Mass", "AOO"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("Mass ~ AOO ~ circNorm")+
  ggthemes::theme_few()

# 15
pdp_lnR_y_interactions[[15]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("Mass", "AOO"),
    BY = "joincount_delta",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("Mass ~ AOO ~ jc")+
  ggthemes::theme_few()

# 16
pdp_lnR_y_interactions[[16]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("Mass", "joincount_delta"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("Mass ~ jc ~ circ")+
  ggthemes::theme_few()

# 17
pdp_lnR_y_interactions[[17]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("D_AOO_a", "AOO"),
    BY = "joincount_delta",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("D ~ AOO ~ jc")+
  ggthemes::theme_few()

# 18
pdp_lnR_y_interactions[[18]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("Mass", "D_AOO_a"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("Mass ~ D ~ circ")+
  ggthemes::theme_few()

# 19
pdp_lnR_y_interactions[[19]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("D_AOO_a", "joincount_delta"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("D ~ jc ~ circ")+
  ggthemes::theme_few()

# 20
pdp_lnR_y_interactions[[20]] <- 
  partial_dep(
    object = fit_3,
    X = train,
    v = c("AOO", "joincount_delta"),
    BY = "circNorm",
    grid_size = 100L,
    strategy = "uniform"
  ) %>%
  plot() +
  ggtitle("AOO ~ jc ~ circ")+
  ggthemes::theme_few()


saveRDS(pdp_lnR_y_interactions, here("Data/output/temp/D_05_lnR_y_pdp_all.rds"))
