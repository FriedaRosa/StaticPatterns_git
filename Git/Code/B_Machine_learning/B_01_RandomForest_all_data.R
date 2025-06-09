#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                B_01_RandomForest_full_models.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#


# Start with clean environment
rm(list = ls())
gc()

# Source 00_Configuration.R
source(here::here("Code/00_Configuration.R"))
package_list <- c(package_list, "ranger", "tidymodels", "caret", "skimr", "DALEXtra")
lapply(package_list, require, character = TRUE)
tidymodels_prefer()
rm(list=ls())
set.seed(123)


library(ggthemes)
library(here)

#----------------------------------------------------------#
# Get final data -----
#----------------------------------------------------------#

dta <- readRDS(here::here("Data/output/1_all_predictors_merged.rds")) %>%
  filter(samplingPeriodID == 1)

# Set variables for model
sp_id <- c("verbatimIdentification", "scientificName")
H1 <- c("Mass", "GlobRangeSize_km2", "Migration", "Habitat_5", "Generalism", "Threatened", "pd")
H2 <- c("D_AOO_a", "mean_lnLac", "AOO", "joincount_delta", "circNorm", "minDist_toBorder_centr")
H3 <- c("datasetID")
predictors <- c(H1, H2, H3)
responses <- c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year")

# filter data to necessary columns
dta_new <- dta %>%
  select(all_of(c(sp_id, responses, H3, H1, H2))) %>%
  ungroup()

#----------------------------------------------------------#
# Create data subsets for both responses  -----
#----------------------------------------------------------#

dta_J <-
  dta_new %>%
  select(-log_R2_1, -log_R2_1_per_year)

dta_lnRR <-
  dta_new %>%
  select(-Jaccard_dissim, -log_R2_1_per_year)

dta_ln_RR_year <-
  dta_new %>%
  select(-Jaccard_dissim, -log_R2_1)

all_datasets <-
  list(dta_J, dta_lnRR, dta_ln_RR_year)

#----------------------------------------------------------#
# Create empty lists to fill in the loop ----
# ---------------------------------------------------------#

tuned_res <- replicate(length(all_datasets), list())
predictions <- replicate(length(all_datasets), list())
res <- replicate(length(all_datasets), list())
list_split_data <- replicate(length(all_datasets), list())


#---------------------------------------------------------#
#
# START LOOP -------------
#
# --------------------------------------------------------#
set.seed(123)
for (data_i in seq_along(all_datasets)){

  set.seed(123)
  model_dd <- all_datasets[[data_i]]

  #----------------------------------------------------------#
  # Create splits  -----
  #----------------------------------------------------------#

  set.seed(123)
  split_info <- initial_split(model_dd, strata = datasetID, prop = 0.8)
  train <- training(split_info)
  test <- testing(split_info)

  set.seed(123)
  folds <- vfold_cv(train, v=10, repeats = 3)

  list_split_data[[data_i]] <- list(split_info = split_info,
                                    train = train,
                                    test = test,
                                    folds = folds)

  #----------------------------------------------------------#
  # Model recipes  -----
  #----------------------------------------------------------#

  if (responses[data_i] == "Jaccard_dissim"){
    mod_rec <-
      recipe(Jaccard_dissim ~ . , data = train) %>%
      update_role(all_of(sp_id),
                  new_role = "speciesID",
                  old_role = "predictor")
  } else if (responses[data_i] == "log_R2_1"){
    mod_rec <-
      recipe(log_R2_1 ~ . , data = train) %>%
      update_role(all_of(sp_id),
                  new_role = "speciesID",
                  old_role = "predictor")
  } else if (responses[data_i] == "log_R2_1_per_year"){
    mod_rec <-
      recipe(log_R2_1_per_year ~ . , data = train) %>%
      update_role(all_of(sp_id),
                  new_role = "speciesID",
                  old_role = "predictor")

  }


  set.seed(123)
  mod_spec <- rand_forest(mtry = tune(),
                          trees = tune(),
                          min_n = tune()) %>%
    set_engine("ranger",
               importance = "permutation",
               respect.unordered.factors = TRUE,
               always.split.variables = "datasetID") %>%
    set_mode("regression")


  #----------------------------------------------------------#
  # Workflow ------
  #----------------------------------------------------------#

  set.seed(123)
  mod_wf <- workflow() %>%
    add_recipe(mod_rec) %>%
    add_model(mod_spec)

  #----------------------------------------------------------#
  # Bayesian hyperparameter tuning -----
  #----------------------------------------------------------#
  set.seed(123)
  rf_params <- parameters(
    mtry(range = c(2L, 10L)),
    min_n(range = c(5L, 15L)),
    trees(range = c(1000L, 5000L))
  )

  ## outcommented to save computing time and read results from file:
  # tictoc::tic()
  # tuned_bayes <- tune_bayes(mod_wf,
  #                           resamples = folds,
  #                           param_info = rf_params,
  #                           initial = 5,
  #                           iter = 50,
  #                           control = control_bayes(verbose = T,
  #                                                   no_improve = 10,
  #                                                   seed = 123))
  # tictoc::toc()
  # saveRDS(tuned_bayes, paste0(here::here("Data/output/temp/B_01_"), responses[data_i], "_tuned_bayes.rds"))

  tuned_bayes <- readRDS(paste0(here::here("Data/output/temp/B_01_"), responses[data_i], "_tuned_bayes.rds"))

  best <- tuned_bayes %>%
    select_best(metric = "rmse")

  tuned_res[[data_i]] <- list(tuned_bayes, best)

  #----------------------------------------------------------#
  # Final model fit -----
  #----------------------------------------------------------#

  final_wf <- mod_wf %>%
    finalize_workflow(best)

  final_fit <- final_wf %>%
    last_fit(split_info) # this will perform the model evaluation on the testing data directly

  saveRDS(final_fit, paste0(here::here("Data/output/temp/B_01_"), responses[data_i], "_final_fit.rds"))

  res[[data_i]] <- final_fit

  #----------------------------------------------------------#
  # Predictions on test data -----
  #----------------------------------------------------------#
  p <-
    predict(final_fit %>% tune::extract_workflow(), new_data = test) %>%
    bind_cols(test %>% select(responses[data_i], verbatimIdentification, scientificName, datasetID)) %>% # change y with response var
    setNames(c(".pred", responses[data_i],  "verbatimIdentification", "scientificName", "datasetID")) %>%
    mutate(resid = .[[responses[data_i]]] - .pred)

  predictions[[data_i]] <- p

}

# xAI:
plot_list <- replicate(3, list())
res_xAI <- replicate(length(res), list())

for (data_i in seq_along(res)){
  dd <- list_split_data[[data_i]]
  fit <- res[[data_i]] %>%
    extract_fit_parsnip()

  explainer <- DALEXtra::explain_tidymodels(
    fit,
    data = dd$train %>% select(!any_of(c(responses, "verbatimIdentification", "scientificName"))),
    y = dd$train %>% pull(any_of(responses)))

  # Partial plots / Model profile
  pd <- model_profile(explainer,
                      groups = "datasetID",
                      type = "partial",
                      N = NULL,
                      variables = explainer$data %>% names())
  # Model diagnostics
  md <- model_diagnostics(explainer)
  # model performance:
  mp <- model_performance(explainer)
  # variable importance:
  vi <- model_parts(explainer)

  # save results to list
  res_xAI[[data_i]] <- list(explainer = explainer,
                            pd = pd,
                            md = md,
                            vi = vi)

  }

plot_list <- replicate(3, list())
plot_list[[1]] <- list(pd1 = plot(res_xAI[[1]]$pd),
                  md1 = plot(res_xAI[[1]]$md %>% mutate(label = paste0(prefix))) + ggthemes::theme_few(),
                  vi1 = plot(res_xAI[[1]]$vi %>% mutate(label = paste0(prefix))) + ggthemes::theme_few())
plot_list[[2]] <- list(pd2 = plot(res_xAI[[2]]$pd),
                       md2 = plot(res_xAI[[2]]$md %>% mutate(label = paste0(prefix))) + ggthemes::theme_few(),
                       vi2 = plot(res_xAI[[2]]$vi %>% mutate(label = paste0(prefix))) + ggthemes::theme_few())
plot_list[[3]] <- list(pd3 = plot(res_xAI[[3]]$pd),
                       md3 = plot(res_xAI[[3]]$md  %>% mutate(label = paste0(prefix))) + ggthemes::theme_few(),
                       vi3 = plot(res_xAI[[3]]$vi %>% mutate(label = paste0(prefix))) + ggthemes::theme_few())



  # Define base path and file name prefix
  out_dir <- here::here("Figures/B_models/xAI/")
  prefix <- responses[data_i]


  # save plots
  plot_list[[data_i]] <- list(
    pd = plot(pd) + ggthemes::theme_few(),
    mp = plot(mp) + ggthemes::theme_few(),
    md = plot(md  %>% mutate(label = paste0(prefix))) + ggthemes::theme_few(),
    vi = plot(vi %>% mutate(label = paste0(prefix))) + ggthemes::theme_few()
  )




  pdf(onefile = TRUE,
      file = here::here(out_dir, paste0("all_data_xAI_plots.pdf")))
  plot_list
  dev.off()




# Interactions with H-stats
library(hstats)

interaction_res <- replicate(length(res), list())
for (data_i in seq_along(list_split_data)){
  # prepare data for hstats function
  dd <- list_split_data[[data_i]]
  fit <- res[[data_i]] %>%
    extract_fit_parsnip()
  preds <- dd$train %>%
    select(!any_of(
      c(responses,
        "verbatimIdentification",
        "scientificName")))

  # run hstats function
  interactions <- hstats(fit,
                         X = preds,
                         threeway_m = 5)

  # save results
  interaction_res[[data_i]] <- interactions

  library(ggthemes)
  library(here)
  library(ggplot2)

  # create the plot
  interaction_plot <- plot(interactions, which = 1:100) + theme_few()

  # save as .svg
  ggsave(filename = paste0(here::here("Figures/B_models/interactions/","B_01_", responses[data_i], "_interactions_hstats.svg")),
         plot = interaction_plot,
         device = "svg",
         width = 10, height = 10)


}


# save RDS
files_to_save <- list(
  interactions = interaction_res,
  res_xAI = res_xAI,
  predictions = predictions,
  tuned_res = tuned_res,
  res = res,
  list_split_data = list_split_data)

saveRDS(files_to_save, here::here("Data/output/B_models/B_01_list_all_results_rf.rds"))




#------------------------------------------------#
# Create Model results cvs
# -----------------------------------------------#


# read back in and save as svg.
list_res <- readRDS(here::here("Data/output/B_models/B_01_list_all_results_rf.rds"))


## save plots (outcommented to avoid overwriting)


# for (resp_i in seq_along(responses)){
# int <- list_res$interactions[[resp_i]]
# class(int) <- "hstats"
# interaction_plot <-
#   plot(int, which = 1:100) +
#   theme_few() +
#   ggtitle(paste(responses[resp_i])) +
#   xlim(c(0,0.14))
# ggsave(filename = here::here("Figures/B_models/interactions/", paste0("B_01_,responses[resp_i], "_interactions_hstats.svg")),
#       plot = interaction_plot,
#       device = "svg",
#       width = 10, height = 10)
# interaction_plot
# }


# save results to csv:

tuned_res_df <- list(list_res$tuned_res[[1]][[2]],
                     list_res$tuned_res[[2]][[2]],
                     list_res$tuned_res[[3]][[2]]) %>%
  bind_rows() %>%
  select(-.config) %>%
  mutate(
    response = c("Jaccard_dissim",
                 "log_R2_1",
                 "log_R2_1_per_year"))


mod_res <-
  list_res$res %>%
  bind_rows() %>%
  .[[3]] %>%
  bind_rows() %>%
  select(-.config, -.estimator) %>%
  mutate(
    response = c("Jaccard_dissim", "Jaccard_dissim",
                 "log_R2_1", "log_R2_1",
                 "log_R2_1_per_year", "log_R2_1_per_year")) %>%
  pivot_wider(id_cols = c("response"),
              names_from = c(".metric"),
              values_from = c(".estimate")) %>%
  full_join(tuned_res_df) %>%
  mutate(cv = "3x10folds",
         tuning_method = "bayesian",
         model_type = "ranger")

# check:
kable(mod_res)

write.csv2(mod_res, here::here("Data/output/results/B_01_ranger_all_data.csv"))

## Interactions:
list_res$interactions[[1]]


int_res_list <- list()
for(resp_i in seq_along(list_res$interactions)){

  if(resp_i == 1) this_response <- "Jaccard_dissim" else
    if(resp_i == 2) this_response <- "log_R2_1" else
      if(resp_i == 3) this_response <- "log_R2_1_per_year"

      int_res_list[[resp_i]] <- data.frame(
        H2 = c(summary(list_res$interactions[[resp_i]])$h2[[1]],
               summary(list_res$interactions[[resp_i]])$h2_overall[[1]],
               summary(list_res$interactions[[resp_i]])$h2_pairwise[[1]],
               summary(list_res$interactions[[resp_i]])$h2_threeway[[1]]),
        response = this_response,
        test = c("total_interaction_strength", rep("overall", 14), rep("pairwise", 10), rep("threeway", 10)),
        variable = c("all",
                     row.names(summary(list_res$interactions[[resp_i]])$h2_overall[[1]]),
                     row.names(summary(list_res$interactions[[resp_i]])$h2_pairwise[[1]]),
                     row.names(summary(list_res$interactions[[resp_i]])$h2_threeway[[1]])

                     ))

}

int_res <- bind_rows(int_res_list)

write.csv2(int_res, here::here("Data/output/results/B_01_Hstats_all_data.csv"))




# Variable importances:
res_all <- readRDS(here::here("Data/output/B_models/B_01_list_all_results_rf.rds"))$res

vimp_list <- list()
for(i in seq_along(res_all)){

  if(i == 1) this_response <- "Jaccard_dissim" else
    if(i == 2) this_response <- "log_R2_1" else
      if(i == 3) this_response <- "log_R2_1_per_year"
  fit <- res_all[[i]] %>% extract_fit_engine()
  vimp <- fit$variable.importance

  vimp_list[[i]] <- data.frame(variable = names(vimp),
                               importance = vimp,
                               row.names = NULL,
                               response = this_response,
                               importance_mode = fit$importance.mode
                               )
}

imp_df <-
  vimp_list %>%
  bind_rows()

write.csv2(imp_df, here::here("Data/output/results/B_01_ranger_vimp_all_data.csv"))
