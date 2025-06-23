#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#              B_03_Model_eval.R
#                
#
#                    Friederike Wölke 
#                           2025
#
#----------------------------------------------------------#


# Start with clean environment

library(dplyr)
library(ggplot2)

# Source 00_Configuration.R
source(here::here("Code/00_Configuration.R"))
package_list <- c(package_list, "ranger", "tidymodels", "caret", "skimr", "DALEXtra")
lapply(package_list, require, character = TRUE)
tidymodels_prefer()
rm(list = ls())
gc()



# Variables:
sp_id <- c("verbatimIdentification", "scientificName")
H1 <- c("Mass", "GlobRangeSize_km2", "Migration", "Habitat_5", "Generalism", "Threatened", "pd")
H2 <- c("D_AOO_a", "mean_lnLac", "AOO", "joincount_delta", "circNorm", "minDist_toBorder_centr")
H3 <- c("datasetID")
predictors <- c(H1, H2, H3)
responses <- c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year")

#----------------------------------------------------------#
#
# Part 1: Mp, Vip, Hstats ----
#
#----------------------------------------------------------#

mp_full <- read.csv2(here::here("Data/output/results/B_01_ranger_all_data.csv"))
mp_split <- read.csv2(here::here("Data/output/results/B_02_ranger_atlas_separate.csv"))

vip_full <- read.csv2(here::here("Data/output/results/B_01_ranger_vimp_all_data.csv"))
vip_split <- read.csv2(here::here("Data/output/results/B_02_ranger_vimp_split_data.csv"))

int_full <-  read.csv2(here::here("Data/output/results/B_01_Hstats_all_data.csv"))
int_split <- read.csv2(here::here("Data/output/results/B_02_Hstats_split_data.csv"))

#----------------------------------------------------------#

p_1 <-
  mp_full %>%
  pivot_longer(cols = c(rmse, rsq)) %>%
  ggplot() +
  geom_point(aes(x = response, y = value, shape = name), cex = 2, col = "black", bg = "black")+
  ggthemes::theme_few()+
  scale_shape_manual(values = c(rmse = 21, rsq = 1))


p_1

p_2 <-
  mp_split %>%
  pivot_longer(cols = c(rmse, rsq)) %>%
  ggplot() +
  geom_point(aes(x = factor(datasetID), y = value, shape = name), cex = 2, col = "black", bg = "black")+
  facet_wrap(~response)+
  ggthemes::theme_few()+
  scale_shape_manual(values = c(rmse = 21, rsq = 1))+
  xlab("DatasetID")


p_2

#----------------------------------------------------------#
library(tidytext)

vip_full <- vip_full %>%
  group_by(variable) %>%
  mutate(mean_importance_jaccard = mean(importance[response == "Jaccard_dissim"], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(variable = fct_reorder(variable, mean_importance_jaccard)) %>%
  mutate(hypo = case_when(variable %in% H1 ~ "H1",
                          variable %in% H2 ~ "H2",
                          variable %in% H3 ~ "H3"))

p_3 <-
  vip_full %>%
  ggplot() +
  geom_col(aes(x = importance, y = variable, fill = hypo)) +
  facet_wrap(~response, scales = "free_x") +  # free_x so y-axis stays fixed
  labs(y = "Variable", x = "Importance")+
  ggthemes::theme_few()


p_3


vip_split <- vip_split %>%
  group_by(variable, datasetID, response) %>%
  mutate(mean_importance_jaccard = mean(importance[response == "Jaccard_dissim"], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(variable = fct_reorder(variable, mean_importance_jaccard)) %>%
  mutate(hypo = case_when(variable %in% H1 ~ "H1",
                          variable %in% H2 ~ "H2",
                          variable %in% H3 ~ "H3"))

p_4 <-
  vip_split %>%
  ggplot() +
  geom_col(aes(x = importance, y = variable, fill = hypo)) +
  facet_wrap(~response+datasetID, scales = "free_x") +  # free_x so y-axis stays fixed
  labs(y = "Variable", x = "Importance")+
  ggthemes::theme_few()


p_4 # saved as 1200x800
#----------------------------------------------------------#
p_5 <-
  int_full %>%
  filter(!c(test == "total_interaction_strength")) %>%
  ggplot()+
  geom_col(aes(x = H2, y = reorder_within(variable, H2, response)))+
  facet_wrap(~test+response,scales = "free")+
  scale_y_reordered()+
  ggthemes::theme_few()

p_5 # saved as 1200x800

p_6 <-
  int_split %>%
  filter(!c(test == "total_interaction_strength")) %>%
  ggplot()+
  geom_col(aes(x = H2, y = reorder_within(variable, H2, response), fill = factor(datasetID)))+
  facet_wrap(~test+response,scales = "free")+
  scale_y_reordered()+
  ggthemes::theme_few()
p_6 # saved as 1500x1500


# save as ppt for postediting
esquisse::ggplot_to_ppt()


#----------------------------------------------------------#
#
# Part 2: Partial dependence plots ----
#
#----------------------------------------------------------#

# Start fresh: 
rm(list = ls())
gc()

# source config file
source(here::here("Code/00_Configuration.R"))

# load libraries
package_list <- 
  c(package_list, 
    "dplyr", "pdp", 
    "ggplot2","ranger", 
    "tidymodels")
lapply(package_list, require, character = TRUE)

# enable preference for tidymodels functions (since there are conflicts)
tidymodels_prefer()

#----------------------------------------------------------#

# Set variables:
H1 <- c("Mass", "GlobRangeSize_km2", "Migration", "Habitat_5", "Generalism", "Threatened", "pd")
H2 <- c("D_AOO_a", "mean_lnLac", "AOO", "joincount_delta", "circNorm", "minDist_toBorder_centr")
H3 <- c("datasetID")
predictors <- c(H1, H2, H3)
temp <- readRDS(here::here("Data/output/1_all_predictors_merged.rds"))

# Create lookup labels for clean plotting:
lookup <- 
  c("Mass" = "Body mass",
    "GlobRangeSize_km2" = "Global range size",
    "Migration" = "Migration",
    "Habitat_5" = "Habitat type",
    "Generalism" = "Generalism",
    "Threatened" = "Threatened",
    "pd" = "Phylogenetic distinctiveness",
    "D_AOO_a" = "Fractal dimension",
    "mean_lnLac" = "Mean log lacunarity",
    "AOO" = "Area of occupancy",
    "joincount_delta" = "Spatial autocorrelation",
    "circNorm" = "Normalized circularity",
    "minDist_toBorder_centr" = "Smallest distance centroid-border") %>% 
  as.data.frame()
lookup$curr_label <- row.names(lookup)
lookup$old_label <- lookup$.
row.names(lookup) <- NULL


# Categorical plot:
cat_vars <- c("Migration",
              "Habitat_5",
              "Generalism",
              "Threatened",
              "datasetID")

#----------------------------------------------------------#

# Model results:
all <- readRDS(here::here("Data/output/B_models/B_01_list_all_results_rf.rds"))
split <- readRDS(here::here("Data/output/B_models/B_02_rf_res_atlas_split.rds"))

#----------------------------------------------------------#

# Extract all fit objects into a named list:
all_models_list <- 
  list(
    # full models:
    Jaccard = all$res[[1]] %>% extract_fit_parsnip(),
    log_Ratio = all$res[[2]] %>% extract_fit_parsnip(),
    log_Ratio_y = all$res[[3]] %>% extract_fit_parsnip(),
    
    # split by datasetID:
    Jaccard_5 = split[[1]][[1]]$final_fit %>% extract_fit_parsnip(),
    Jaccard_6 = split[[2]][[1]]$final_fit %>% extract_fit_parsnip(),
    Jaccard_13 = split[[3]][[1]]$final_fit %>% extract_fit_parsnip(),
    Jaccard_26 = split[[4]][[1]]$final_fit %>% extract_fit_parsnip(),
    
    log_Ratio_5 = split[[1]][[2]]$final_fit %>% extract_fit_parsnip(),
    log_Ratio_6 = split[[2]][[2]]$final_fit %>% extract_fit_parsnip(),
    log_Ratio_13 = split[[3]][[2]]$final_fit %>% extract_fit_parsnip(),
    log_Ratio_26 = split[[4]][[2]]$final_fit %>% extract_fit_parsnip(),
    
    log_Ratio_y_5 = split[[1]][[3]]$final_fit %>% extract_fit_parsnip(),
    log_Ratio_y_6 = split[[2]][[3]]$final_fit %>% extract_fit_parsnip(),
    log_Ratio_y_13 = split[[3]][[3]]$final_fit %>% extract_fit_parsnip(),
    log_Ratio_y_26 =  split[[4]][[3]]$final_fit %>% extract_fit_parsnip()
    )

# Extract all training data to a list (note: these are identical for different responses)
train_data <- 
  list(
    # full data:
    all = all$list_split_data[[1]]$train,
    # split by dataset:
    split_5 = split[[1]][[1]]$final_fit$splits[[1]] %>% training(),
    split_6 = split[[2]][[1]]$final_fit$splits[[1]] %>% training(),
    split_13 = split[[3]][[1]]$final_fit$splits[[1]] %>% training(),
    split_26 = split[[4]][[1]]$final_fit$splits[[1]] %>% training())

#----------------------------------------------------------#


# create empty list:
pdp_list <- replicate(15, list)

# testing loop:
model_i <- 1
var_i <- 1
#----------------------------------------------------------#
for (model_i in seq_along(all_models_list)){
  
  if (model_i %in% c(1,4,5,6,7)) response <- "Jaccard_dissim" else
    if (model_i %in% c(2,8,9,10,11)) response <- "log_R2_1" else
      if (model_i %in% c(3,12,13,14,15)) response <- "log_R2_1_per_year"
      
  # select model
  this_model <- all_models_list[[model_i]]
  
  # dynamically remove datasetID as predictor depending on model
  if(model_i %in% c(4:15)) predictors <- predictors[1:13] else predictors

  # create empty list for each model_i
  pdp_list[[model_i]] <- list()
  
  # match the model to its train_data 
  # (it was previously saved in the model results list)
  if (model_i %in% c(1:3)) train_i <- 1 else
    if (model_i %in% c(4:7)) train_i <- 2 else
      if (model_i %in% c(8:11)) train_i <- 3 else
        if(model_i %in% c(12:15)) train_i <- 4
  
  if (model_i %in% c(1,2,3)) obj_ident <- "full" else
    if (model_i %in% c(4:7)) obj_ident <- "Jaccard_split" else
      if (model_i %in% c(8:11)) obj_ident <- "log_R2_1_split" else
        if (model_i %in% c(12:15)) obj_ident <- "log_R2_1_y_ split"
  
  # set current training data
  this_train <- train_data[[train_i]]

  #----------------------------------------------------------#
  # loop through all variables in the data (all predictors)
  #----------------------------------------------------------#
  
  for (var_i in seq_along(predictors)){

    # select current variable to plot
    this_variable <- predictors[var_i]
    print(paste0("object = ", obj_ident, ", Response = ", response, ", Var = ", this_variable))
  
  
    pdp_list[[model_i]][[var_i]] <- 
      pdp::partial(
        object = this_model,
        train = this_train,
        pred.var = this_variable,
        plot = FALSE)
  }
}

# save to temporary folder
saveRDS(pdp_list, here::here("Data/output/temp/B_04_partial_dependencies.rds"))

# read back in
pdp_list <- readRDS( here::here("Data/output/temp/B_04_partial_dependencies.rds"))


# Merge to plotting df:
pdp_df_list <- replicate(15, list)
pdp_df_list_cat <- replicate(15, list)

for (model_i in seq_along (pdp_list)){
  this_model <- pdp_list[[model_i]]
  pdp_df_list[[model_i]] <- list()
  pdp_df_list_cat[[model_i]] <-list()
  
  if (model_i %in% c(1,4,5,6,7)) response <- "Jaccard_dissim" else
    if (model_i %in% c(2,8,9,10,11)) response <- "log_R2_1" else
      if (model_i %in% c(3,12,13,14,15)) response <- "log_R2_1_per_year"
  
  if (model_i %in% c(1,2,3)) obj_ident <- "full" else
    if (model_i %in% c(4:7)) obj_ident <- "Jaccard_split" else
      if (model_i %in% c(8:11)) obj_ident <- "log_R2_1_split" else
        if (model_i %in% c(12:15)) obj_ident <- "log_R2_1_y_split"
  
  if (model_i %in% c(1,2,3)) datasetID <- "all" else
    if (model_i %in% c(4,8,12)) datasetID <- "5" else
      if(model_i %in% c(5,9,13)) datasetID <- "6" else
        if(model_i %in% c(6,10,14)) datasetID <- "13" else
          if(model_i %in% c(7,11,15)) datasetID <- "26"
  
  for(var_i in seq_along(this_model)){
    
    
    pdp_df_i <- this_model[[var_i]]
    pdp_df_i$variable <- names(pdp_df_i)[1]
    
    pdp_df_i$value <- pdp_df_i[[1]]
    pdp_df_i$response <- response
    pdp_df_i$datasetID <- datasetID
    pdp_df_i$model <- obj_ident
    pdp_df_i[[1]] <- NULL
    
    if (unique(pdp_df_i$variable) %in% cat_vars) {pdp_df_list_cat[[model_i]][[var_i]] <- pdp_df_i} else {pdp_df_list[[model_i]][[var_i]] <- pdp_df_i}
    }
  }

pdp_num <- pdp_df_list %>% bind_rows()
pdp_cat <- pdp_df_list_cat %>% bind_rows()

pdp_res <- list(numeric = pdp_num,
                categorical = pdp_cat)



#--------------------------------------------------#
# Plots (ggplot2) ----
#--------------------------------------------------#

num_cat_i <- 1
model_i <- 1
resp_i <- 1

lookup <- 
  c("Mass" = "Body mass",
    "GlobRangeSize_km2" = "Global range size",
    "Migration" = "Migration",
    "Habitat_5" = "Habitat type",
    "Generalism" = "Generalism",
    "Threatened" = "Threatened",
    "pd" = "Phylogenetic distinctiveness",
    "D_AOO_a" = "Fractal dimension",
    "mean_lnLac" = "Mean log lacunarity",
    "AOO" = "Area of occupancy",
    "joincount_delta" = "Spatial autocorrelation",
    "circNorm" = "Normalized circularity",
    "minDist_toBorder_centr" = "Smallest distance to the border") 




for (num_cat_i in seq_along(pdp_res)){ # there are 2 sub-lists (num, cat)
  
  if (num_cat_i == 1) { # numeric
    
    # chose sublist from results
    pdp_subset <- pdp_res[[num_cat_i]]
    
    # extract results from models dynamically
    for (model_i in seq_along(unique(pdp_subset$model))){ # there are 4 models
      
      this_model <- unique(pdp_subset$model)[model_i]
      model_dd <- pdp_subset %>% filter(model == this_model)
      
      for (resp_i in seq_along(unique(model_dd$response))){
        
        this_resp <- unique(model_dd$response)[resp_i]
        
        print(paste0("Num/Cat = ",num_cat_i,"; model = ",this_model,"; response = ", this_resp))
        
        response_dd <- model_dd %>% 
          filter(response == this_resp)
        
        p0 <- response_dd %>%
          ggplot()+
          geom_line(aes(
            y = yhat,
            x = value,
            col = datasetID,
            group = interaction(variable, datasetID)),
            linewidth = 0.8, 
            alpha = 0.8)+
          labs(x = NULL,
               y = bquote(italic(.(this_resp)) ~ "(Partial)"))+
          facet_wrap(~ variable , scales = "free_x", labeller = labeller(variable = lookup))+
          ggthemes::theme_par()+
          ggtitle(paste(this_resp), paste0(this_model, "-model"))+
          scale_color_brewer(palette = "PuOr")
        
        
        # adjust axes and plot hline
      
        if (this_resp == "Jaccard_dissim") {
          p1 <-
            p0 +
            ylim(0, 1) +
            geom_hline(
              yintercept = 0.5,
              linewidth = 0.3,
              col = "darkgrey"
            )
        } else if (this_resp == "log_R2_1") {
          p1 <-
            p0 +
            #ylim(min(temp$log_R2_1), max(temp$log_R2_1)) +
            geom_hline(
              yintercept = 0,
              linewidth = 0.3,
              col = "darkgrey"
            )
          
        } else if (this_resp == "log_R2_1_per_year") {
          p1 <-
            p0 +
            #ylim(min(temp$log_R2_1_per_year), max(temp$log_R2_1_per_year)) +
            geom_hline(
              yintercept = 0,
              linewidth = 0.3,
              col = "darkgrey"
            )
        }
        
        print(p1)
        
        } # responses loop closing
      
      } # models loop closing
    
    } else { # Categorical variables 
  
    # chose sublist from results
    pdp_subset <- pdp_res[[num_cat_i]]
    
    # extract results from models dynamically
    for (model_i in seq_along(unique(pdp_subset$model))){ # there are 4 models
      
      this_model <- unique(pdp_subset$model)[model_i]
      model_dd <- pdp_subset %>% filter(model == this_model)
      
      for (resp_i in seq_along(unique(model_dd$response))){
        
        this_resp <- unique(model_dd$response)[resp_i]
        
        print(paste0("Num/Cat = ",num_cat_i,"; model = ",this_model,"; response = ", this_resp))
        
        response_dd <- model_dd %>% 
          filter(response == this_resp)
        
        p0 <- response_dd %>%
          ggplot()+
          geom_point(
            aes(y = yhat,
                x = value,
                col = datasetID),
            alpha = 0.8)+
          ggthemes::theme_par()+
          labs(y = bquote(italic(.(this_resp)) ~ "(Partial)"))+
          facet_wrap(~ variable , scales = "free_x", labeller = labeller(variable = lookup))+
          scale_color_brewer(palette = "PuOr")
        
        
        # adjust axes and plot hline
        
        if (this_resp == "Jaccard_dissim") {
          p1 <-
            p0 +
            ylim(0, 1) +
            geom_hline(
              yintercept = 0.5,
              linewidth = 0.3,
              col = "darkgrey"
            )
        } else if (this_resp == "log_R2_1") {
          p1 <-
            p0 +
            #ylim(min(temp$log_R2_1), max(temp$log_R2_1)) +
            geom_hline(
              yintercept = 0,
              linewidth = 0.3,
              col = "darkgrey"
            )
          
        } else if (this_resp == "log_R2_1_per_year") {
          p1 <-
            p0 +
            #ylim(min(temp$log_R2_1_per_year), max(temp$log_R2_1_per_year)) +
            geom_hline(
              yintercept = 0,
              linewidth = 0.3,
              col = "darkgrey"
            )
        }
        
        print(p1)
        
        } # responses loop closing
      } # models loop closing
    } # categorical variables loop closing
  } # closing the whole thing






## For maintext:

# Jaccard per model:

plots_numeric <- list()
for (var_i in seq_along(unique(pdp_res[[1]]$variable))){
  this_var <- unique(pdp_res[[1]]$variable)[var_i]
  this_resp <- "Jaccard_dissim"
  this_model <- "Jaccard_split"
  # Numeric Jaccard plots

  p0 <-pdp_res[[1]] %>%
    filter(model == "Jaccard_split" & variable == this_var) %>%
    ggplot()+
    geom_hline(yintercept = 0.5, linewidth = 0.3, col = "darkgrey")+
    geom_line(
      aes(
        x = value,
        y = yhat,
        group = datasetID,
        col = datasetID
      ),
      linewidth = 0.8
    )+
    ggthemes::theme_few()+
    scale_color_manual(values = c("all" = "black",
                                  "5" = "#E66101", 
                                  "6" = "#FDB863" ,
                                  "13" = "#B2ABD2", 
                                  "26" = "#5E3C99"))+
    labs(x = lookup[this_var],
         y = bquote(italic(.(this_resp)) ~ "(Partial)"),
         col = "Dataset")+
    ggtitle(NULL,paste0(this_model, "-model"))+
    ylim(0,1)
  
  p0
  plots_numeric[[var_i]] <- p0
  #ggsave(plot = last_plot(), path = here::here("Figures/B_models/performance"), filename = paste0("B_04_PartialPlot_Jaccard_split_",this_var,".svg" ))
}

plots_numeric

ggsave(plot = plots_numeric[[1]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","mass",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[2]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","globalrangesize",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[3]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","phylodist",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[4]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","D",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[5]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","Lacunarity",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[6]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","AOO",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[7]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","SAC",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[8]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","Circ",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[9]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","DistBorder",".svg"),
       width = 3.5, height = 2.8)



plots_categorical <- list()  
for (var_i in seq_along(unique(pdp_res[[2]]$variable))){
  this_var <- unique(pdp_res[[2]]$variable)[var_i]
  this_resp <- "Jaccard_dissim"
  this_model <- "Jaccard_split"
  # Categorical Jaccard plots
  
  p1 <- pdp_res[[2]] %>%
    filter(model == "Jaccard_split" & variable == this_var) %>%
    ggplot()+
    geom_hline(yintercept = 0.5, linewidth = 0.3, col = "darkgrey")+
    geom_point(
      aes(
        x = value,
        y = yhat,
        col = datasetID
      ),
      size = 3,
      alpha = 0.75
    )+
    ggthemes::theme_few()+
    scale_color_brewer(palette = "PuOr")+
    labs(x = lookup[this_var],
         y = bquote(italic(.(this_resp)) ~ "(Partial)"),
         col = "Dataset")+
    ggtitle(NULL,paste0(this_model, "-model"))+
    ylim(0,1)
  
  plots_categorical[[var_i]] <- p1
  #ggsave(plot = last_plot(), device = "svg", path = here::here("Figures/B_models/performance"), filename = paste0("B_04_PartialPlot_Jaccard_split_",this_var,".svg" ))
}

plots_categorical <- plots_categorical[c(1:4)]


ggsave(plot = plots_categorical[[1]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","Migration",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_categorical[[2]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","Habitat",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_categorical[[3]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","Generalism",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_categorical[[4]], 
       path = here::here("Figures/B_models/performance/J/"), 
       filename = paste0("B_04_PartialPlot_Jaccard_split_","Threat",".svg"),
       width = 3.5, height = 2.8)




patchwork::wrap_plots(plots_numeric, ncol = 3) + 
  patchwork::plot_layout(tag_level = 'keep', guides = 'collect', axes = "collect_y") +
  patchwork::plot_annotation(tag_levels = 'a', tag_suffix = ')')

export::graph2ppt(
  file=here::here("Figures/B_models/performance/B_04_PartialSep_Jaccard.ppt"), width=8, height=10)

patchwork::wrap_plots(plots_categorical, ncol = 3) + 
  patchwork::plot_layout(tag_level = 'keep', guides = 'collect', axes = "collect_y") +
  patchwork::plot_annotation(tag_levels = 'a', tag_suffix = ')')

export::graph2ppt(
  file=here::here("Figures/B_models/performance/B_04_PartialSep_Categorical_Jaccard.ppt"), width=8, height=(10/3)*2)





# Log Ratio per model:
plots_numeric <- list()
for (var_i in seq_along(unique(pdp_res[[1]]$variable))){
  this_var <- unique(pdp_res[[1]]$variable)[var_i]
  this_resp <- "log_R2_1"
  this_model <- "log_R2_1_split"
  # Numeric Jaccard plots
  
  p0 <-pdp_res[[1]] %>%
    filter(model == "log_R2_1_split" & variable == this_var) %>%
    ggplot()+
    geom_hline(yintercept = 0, linewidth = 0.3, col = "darkgrey")+
    geom_line(
      aes(
        x = value,
        y = yhat,
        group = datasetID,
        col = datasetID
      ),
      linewidth = 0.8
    )+
    ggthemes::theme_few()+
    scale_color_manual(values = c("all" = "black",
                                  "5" = "#E66101", 
                                  "6" = "#FDB863" ,
                                  "13" = "#B2ABD2", 
                                  "26" = "#5E3C99"))+
    labs(x = lookup[this_var],
         y = bquote(italic(.(this_resp)) ~ "(Partial)"),
         col = "Dataset")+
    ggtitle(NULL,paste0(this_model, "-model"))
  
  p0
  plots_numeric[[var_i]] <- p0
  #ggsave(plot = last_plot(), path = here::here("Figures/B_models/performance"), filename = paste0("B_04_PartialPlot_Jaccard_split_",this_var,".svg" ))
}

plots_numeric


plots_categorical <- list()  
for (var_i in seq_along(unique(pdp_res[[2]]$variable))){
  this_var <- unique(pdp_res[[2]]$variable)[var_i]
  this_resp <- "log_R2_1"
  this_model <- "log_R2_1_split"
  # Categorical Jaccard plots
  
  p1 <- pdp_res[[2]] %>%
    filter(model == "log_R2_1_split" & variable == this_var) %>%
    ggplot()+
    geom_hline(yintercept = 0, linewidth = 0.3, col = "darkgrey")+
    geom_point(
      aes(
        x = value,
        y = yhat,
        col = datasetID
      ),
      size = 3,
      alpha = 0.75
    )+
    ggthemes::theme_few()+
    scale_color_brewer(palette = "PuOr")+
    labs(x = lookup[this_var],
         y = bquote(italic(.(this_resp)) ~ "(Partial)"),
         col = "Dataset")+
    ggtitle(NULL,paste0(this_model, "-model"))
  plots_categorical[[var_i]] <- p1
  #ggsave(plot = last_plot(), device = "svg", path = here::here("Figures/B_models/performance"), filename = paste0("B_04_PartialPlot_Jaccard_split_",this_var,".svg" ))
}

plots_categorical <- plots_categorical[c(1:4)]




ggsave(plot = plots_numeric[[1]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","mass",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[2]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","globalrangesize",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[3]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","phylodist",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[4]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","D",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[5]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","Lacunarity",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[6]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","AOO",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[7]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","SAC",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[8]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","Circ",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[9]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","DistBorder",".svg"),
       width = 3.5, height = 2.8)

ggsave(plot = plots_categorical[[1]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","Migration",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_categorical[[2]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","Habitat",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_categorical[[3]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","Generalism",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_categorical[[4]], 
       path = here::here("Figures/B_models/performance/lnR/"), 
       filename = paste0("B_04_PartialPlot_lnR_split_","Threat",".svg"),
       width = 3.5, height = 2.8)






patchwork::wrap_plots(plots_numeric, ncol = 3) + 
  patchwork::plot_layout(tag_level = 'keep', guides = 'collect', axes = "collect_y") +
  patchwork::plot_annotation(tag_levels = 'a', tag_suffix = ')')

export::graph2ppt(
  file=here::here("Figures/B_models/performance/B_04_PartialSep_logR.ppt"), width=8, height=10)

patchwork::wrap_plots(plots_categorical, ncol = 3) + 
  patchwork::plot_layout(tag_level = 'keep', guides = 'collect', axes = "collect_y") +
  patchwork::plot_annotation(tag_levels = 'a', tag_suffix = ')')

export::graph2ppt(
  file=here::here("Figures/B_models/performance/B_04_PartialSep_Categorical_logR.ppt"), width=8, height=(10/3)*2)








# log ratio per year per model:

plots_numeric <- list()
for (var_i in seq_along(unique(pdp_res[[1]]$variable))){
  this_var <- unique(pdp_res[[1]]$variable)[var_i]
  this_resp <- "log_R2_1_per_year"
  this_model <- "log_R2_1_y_split"
  # Numeric Jaccard plots
  
  p0 <-pdp_res[[1]] %>%
    filter(model == "log_R2_1_y_split" & variable == this_var) %>%
    ggplot()+
    geom_hline(yintercept = 0, linewidth = 0.3, col = "darkgrey")+
    geom_line(
      aes(
        x = value,
        y = yhat,
        group = datasetID,
        col = datasetID
      ),
      linewidth = 0.8
    )+
    ggthemes::theme_few()+
    scale_color_manual(values = c("all" = "black",
                                  "5" = "#E66101", 
                                  "6" = "#FDB863" ,
                                  "13" = "#B2ABD2", 
                                  "26" = "#5E3C99"))+
    labs(x = lookup[this_var],
         y = bquote(italic(.(this_resp)) ~ "(Partial)"),
         col = "Dataset")+
    ggtitle(NULL,paste0(this_model, "-model"))
  
  p0
  plots_numeric[[var_i]] <- p0
  #ggsave(plot = last_plot(), path = here::here("Figures/B_models/performance"), filename = paste0("B_04_PartialPlot_Jaccard_split_",this_var,".svg" ))
}

plots_numeric

plots_categorical <- list()  
for (var_i in seq_along(unique(pdp_res[[2]]$variable))){
  this_var <- unique(pdp_res[[2]]$variable)[var_i]
  this_resp <- "log_R2_1_per_year"
  this_model <- "log_R2_1_y_split"
  # Categorical Jaccard plots
  
  p1 <- pdp_res[[2]] %>%
    filter(model == "log_R2_1_y_split" & variable == this_var) %>%
    ggplot()+
    geom_hline(yintercept = 0, linewidth = 0.3, col = "darkgrey")+
    geom_point(
      aes(
        x = value,
        y = yhat,
        col = datasetID
      ),
      size = 3,
      alpha = 0.75
    )+
    ggthemes::theme_few()+
    scale_color_brewer(palette = "PuOr")+
    labs(x = lookup[this_var],
         y = bquote(italic(.(this_resp)) ~ "(Partial)"),
         col = "Dataset")+
    ggtitle(NULL,paste0(this_model, "-model"))
  
  plots_categorical[[var_i]] <- p1
  #ggsave(plot = last_plot(), device = "svg", path = here::here("Figures/B_models/performance"), filename = paste0("B_04_PartialPlot_Jaccard_split_",this_var,".svg" ))
}

plots_categorical <- plots_categorical[c(1:4)]






ggsave(plot = plots_numeric[[1]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","mass",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[2]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","globalrangesize",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[3]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","phylodist",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[4]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","D",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[5]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","Lacunarity",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[6]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","AOO",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[7]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","SAC",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[8]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","Circ",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_numeric[[9]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","DistBorder",".svg"),
       width = 3.5, height = 2.8)

ggsave(plot = plots_categorical[[1]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","Migration",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_categorical[[2]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","Habitat",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_categorical[[3]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","Generalism",".svg"),
       width = 3.5, height = 2.8)
ggsave(plot = plots_categorical[[4]], 
       path = here::here("Figures/B_models/performance/lnR_y/"), 
       filename = paste0("B_04_PartialPlot_lnR_y_split_","Threat",".svg"),
       width = 3.5, height = 2.8)








patchwork::wrap_plots(plots_numeric, ncol = 3) + 
  patchwork::plot_layout(tag_level = 'keep', guides = 'collect', axes = "collect_y") +
  patchwork::plot_annotation(tag_levels = 'a', tag_suffix = ')')

export::graph2ppt(
  file=here::here("Figures/B_models/performance/B_04_PartialSep_LnR_y.ppt"), width=8, height=10)

patchwork::wrap_plots(plots_categorical, ncol = 3) + 
  patchwork::plot_layout(tag_level = 'keep', guides = 'collect', axes = "collect_y") +
  patchwork::plot_annotation(tag_levels = 'a', tag_suffix = ')')

export::graph2ppt(
  file=here::here("Figures/B_models/performance/B_04_PartialSep_Categorical_LnR_y.ppt"), width=8, height=(10/3)*2)








## ## Moved also to Script D_05_Figure_S12 ###
## Interactions & their partial dependencies ##
library(hstats)
library(here)
library(dplyr)
library(tidymodels)
library(ranger)
data <- readRDS(here::here("Data/output/B_models/B_01_list_all_results_rf.rds"))

fit_1 <- data$res[[1]] %>% extract_fit_parsnip()
fit_2 <- data$res[[2]] %>% extract_fit_parsnip()
fit_3 <- data$res[[3]] %>% extract_fit_parsnip()

partial_dep(object = fit_1,
            v = "mean_lnLac",
            X = data$list_split_data[[1]]$train,
            BY = "joincount_delta") %>%
  plot(show_points = FALSE)+
  ylim(0,1)



# Heatmap
partial_dep(object = fit_1,
            v = c("mean_lnLac", "joincount_delta"),
            X = data$list_split_data[[1]]$train,
            #BY = "joincount_delta",
            grid_size = 100) %>%
  plot()


# Heatmap
partial_dep(object = fit_1,
            v = c("D_AOO_a", "AOO"),
            X = data$list_split_data[[1]]$train,
            #BY = "AOO",
            grid_size = 1000) %>%
  plot()




# Log Ratio: Body mass ~ AOO 
# Heatmap
partial_dep(object = fit_2,
            v = c("Mass", "AOO"),
            X = data$list_split_data[[1]]$train,
            grid_size = 100) %>%
  plot()

partial_dep(object = fit_3,
            v = c("Mass", "Habitat_5"),
            BY = "datasetID",
            X = data$list_split_data[[1]]$train,
            grid_size = 100) %>%
  plot()


partial_dep(object = fit_3,
            v = c("AOO", "D_AOO_a"),
            BY = "datasetID",
            X = data$list_split_data[[1]]$train,
            grid_size = 100) %>%
  plot()
