#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#              B_02_RandomForest_separate_per_atlas.R
#                
#
#                    Friederike Wölke 
#                           2025
#
#----------------------------------------------------------#


# Start with clean environment
rm(list = ls())
gc()

# Source 00_Configuration.R
source(here::here("Code/00_Configuration.R"))
package_list <- c(package_list, "ranger", "tidymodels", "caret", "skimr", "DALEXtra", "hstats")
lapply(package_list, require, character = TRUE)
tidymodels_prefer()
rm(list = ls())
set.seed(123)

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


# modify data for modeling
dta_new <- dta %>%
  select(all_of(c(sp_id, responses, H3, H1, H2))) %>%
  ungroup()


# check correct coding
str(dta_new) # 14 predictors

# split by dataset
dta_atlas_split <- dta_new %>%
  group_split(datasetID)

#----------------------------------------------------------#
set.seed(123)
mod_spec <- rand_forest(mtry = tune(),
                        trees = tune(),
                        min_n = tune()) %>%
  set_engine("ranger",
             importance = "permutation",
             respect.unordered.factors = TRUE) %>%
  set_mode("regression")

rf_params <- parameters(
  mtry(range = c(2L, 10L)),
  min_n(range = c(5L, 15L)),
  trees(range = c(1000L, 5000L))
)

#----------------------------------------------------------#
tuned_res <- replicate(length(dta_atlas_split), list())
predictions <- replicate(length(dta_atlas_split), list())
res <- replicate(length(dta_atlas_split), list())
list_split_data <- replicate(length(dta_atlas_split), list())
results_response <- replicate(3, list())
results_atlas_i <- replicate(3, list())
#----------------------------------------------------------#
for (atlas_i in seq_along(dta_atlas_split)) {
  this_atlas <- dta_atlas_split[[atlas_i]] %>%
    select(-datasetID)

  # split data
  set.seed(123)
  split_info <- initial_split(this_atlas, prop = 0.8)
  train <- training(split_info)
  test <- testing(split_info)

  set.seed(123)
  folds <- vfold_cv(train, v = 10, repeats = 3)

  list_split_data[[atlas_i]] <-
    list(
      split_info = split_info,
      train = train,
      test = test,
      folds = folds
    )

  for (resp_i in seq_along(responses)) {

    resp <- responses[resp_i]

    if (responses[resp_i] == "Jaccard_dissim") {
      mod_rec <-
        recipe(Jaccard_dissim ~ ., data = train) %>%
        update_role(all_of(sp_id),
          new_role = "speciesID",
          old_role = "predictor"
        ) %>%
        update_role(all_of(responses[-resp_i]),
          old_role = "predictor",
          new_role = "alt_resp"
        )
    } else if (responses[resp_i] == "log_R2_1") {
      mod_rec <-
        recipe(log_R2_1 ~ ., data = train) %>%
        update_role(all_of(sp_id),
          new_role = "speciesID",
          old_role = "predictor"
        ) %>%
        update_role(all_of(responses[-resp_i]),
          old_role = "predictor",
          new_role = "alt_resp"
        )
    } else if (responses[resp_i] == "log_R2_1_per_year") {
      mod_rec <-
        recipe(log_R2_1_per_year ~ ., data = train) %>%
        update_role(all_of(sp_id),
          new_role = "speciesID",
          old_role = "predictor"
        ) %>%
        update_role(all_of(responses[-resp_i]),
          old_role = "predictor",
          new_role = "alt_resp"
        )
    }

  set.seed(123)
  mod_wf <- workflow() %>%
    add_recipe(mod_rec) %>%
    add_model(mod_spec)

  #----------------------------------------------------------#
  # Bayesian hyperparameter tuning -----
  #----------------------------------------------------------#


  tictoc::tic()
  set.seed(123)
  tuned_bayes <-
    tune_bayes(mod_wf,
                            resamples = folds,
                            param_info = rf_params,
                            initial = 5,
                            iter = 50,
                            control = control_bayes(verbose = T,
                                                    no_improve = 10,
                                                    seed = 123))
  tictoc::toc()

  best <-
    tuned_bayes %>%
    select_best(metric = "rmse")

  #----------------------------------------------------------#
  # Final model fit -----
  #----------------------------------------------------------#

  final_wf <- mod_wf %>%
    finalize_workflow(best)

  final_fit <- final_wf %>%
    last_fit(split_info) # this will perform the model evaluation on the testing data directly

  #----------------------------------------------------------#
  # Predictions on test data -----
  #----------------------------------------------------------#

  p <-
    predict(final_fit %>% tune::extract_workflow(), new_data = test) %>%
    bind_cols(test %>% select(responses[resp_i], verbatimIdentification, scientificName)) %>% # change y with response var
    setNames(c(".pred", responses[resp_i],  "verbatimIdentification", "scientificName")) %>%
    mutate(resid = .[[responses[resp_i]]] - .pred)


  #----------------------------------------------------------#
  # explainable AI (xAI) -----
  #----------------------------------------------------------#


  fit <- final_fit %>%
    extract_fit_parsnip()

  explainer <- DALEXtra::explain_tidymodels(
    label = paste0(responses[resp_i],"_",atlas_i),
    fit,
    data = train %>% select(!any_of(c(responses, "verbatimIdentification", "scientificName"))),
    y = train %>% pull(responses[resp_i]))

  # save results to list
  results_response[[resp_i]] <-
    list(
      tuned_bayes = tuned_bayes,
      best = best,
      final_fit = final_fit,
      predictions = p,
      explainer = explainer,
      pd = model_profile(explainer,
        type = "partial",
        N = NULL,
        variables = explainer$data %>% names()
      ),
      md = model_diagnostics(explainer),
      mp = model_performance(explainer),
      vi = model_parts(explainer)
    )    }


  results_atlas_i[[atlas_i]] <- results_response


  }

saveRDS(results_atlas_i, here::here("Data/output/B_models/B_02_rf_res_atlas_split.rds"))

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# read back in
results_atlas_i <- readRDS(here::here("Data/output/B_models/B_02_rf_res_atlas_split.rds"))
responses <- c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year")


#------------------------------------------- #
# Plot xAI ---------
# -------------------------------------------#

plot_list <- replicate(2, list())

for (atlas_i in seq_along(results_atlas_i)){
  plot_list[[atlas_i]] <- list()
  for (resp_i in seq_along(results_atlas_i[[atlas_i]])){

    this_response <- responses[resp_i]


    if (atlas_i == 1) datasetID <- 5 else
      if (atlas_i == 2) datasetID <- 6 else
        if (atlas_i == 3) datasetID <- 13 else
          if (atlas_i == 4) datasetID <- 26


      subset <- results_atlas_i[[atlas_i]][[resp_i]]

      # List with all plots for saving as pdf
      plot_list[[atlas_i]][[resp_i]] <- list(
        pd = plot(subset$pd) + ggthemes::theme_few(),
        mp = plot(subset$mp) + ggthemes::theme_few(),
        md = plot(subset$md  %>% mutate(label = paste0(datasetID,"_", this_response))) + ggthemes::theme_few(),
        vi = plot(subset$vi %>% mutate(label = paste0(datasetID,"_", this_response))) + ggthemes::theme_few()
        )

      }

}

pdf(file = paste0(here::here("Figures/B_models/xAI/B_02_xAI_separate_per_atlas.pdf")), onefile = TRUE)
plot_list
dev.off()


#------------------------------------------- #
# Save ranger res as csv ---------
# -------------------------------------------#
tuned_res_df <-
  list(
    results_atlas_i[[1]][[1]]$best,
    results_atlas_i[[1]][[2]]$best,
    results_atlas_i[[1]][[3]]$best,
    results_atlas_i[[2]][[1]]$best,
    results_atlas_i[[2]][[2]]$best,
    results_atlas_i[[2]][[3]]$best,
    results_atlas_i[[3]][[1]]$best,
    results_atlas_i[[3]][[2]]$best,
    results_atlas_i[[3]][[3]]$best,
    results_atlas_i[[4]][[1]]$best,
    results_atlas_i[[4]][[2]]$best,
    results_atlas_i[[4]][[3]]$best
  ) %>%
  bind_rows() %>%
  select(-.config) %>%
  mutate(
    response = c(
      "Jaccard_dissim",
      "log_R2_1",
      "log_R2_1_per_year",
      "Jaccard_dissim",
      "log_R2_1",
      "log_R2_1_per_year",
      "Jaccard_dissim",
      "log_R2_1",
      "log_R2_1_per_year",
      "Jaccard_dissim",
      "log_R2_1",
      "log_R2_1_per_year"
    ),
    datasetID = c(5, 5, 5, 6, 6, 6, 13, 13, 13, 26, 26, 26)
  )


mod_res <-
  list(
    results_atlas_i[[1]][[1]]$final_fit$.metrics,
    results_atlas_i[[1]][[2]]$final_fit$.metrics,
    results_atlas_i[[1]][[3]]$final_fit$.metrics,
    results_atlas_i[[2]][[1]]$final_fit$.metrics,
    results_atlas_i[[2]][[2]]$final_fit$.metrics,
    results_atlas_i[[2]][[3]]$final_fit$.metrics,
    results_atlas_i[[3]][[1]]$final_fit$.metrics,
    results_atlas_i[[3]][[2]]$final_fit$.metrics,
    results_atlas_i[[3]][[3]]$final_fit$.metrics,
    results_atlas_i[[4]][[1]]$final_fit$.metrics,
    results_atlas_i[[4]][[2]]$final_fit$.metrics,
    results_atlas_i[[4]][[3]]$final_fit$.metrics) %>%
  bind_rows() %>%
  mutate(
    response = c(
      "Jaccard_dissim",
      "Jaccard_dissim",
      "log_R2_1",
      "log_R2_1",
      "log_R2_1_per_year",
      "log_R2_1_per_year",
      "Jaccard_dissim",
      "Jaccard_dissim",
      "log_R2_1",
      "log_R2_1",
      "log_R2_1_per_year",
      "log_R2_1_per_year",
      "Jaccard_dissim",
      "Jaccard_dissim",
      "log_R2_1",
      "log_R2_1",
      "log_R2_1_per_year",
      "log_R2_1_per_year",
      "Jaccard_dissim",
      "Jaccard_dissim",
      "log_R2_1",
      "log_R2_1",
      "log_R2_1_per_year",
      "log_R2_1_per_year"),
    datasetID = c(5, 5, 5, 5, 5, 5,
                  6, 6, 6, 6, 6, 6,
                  13, 13, 13, 13, 13, 13,
                  26, 26, 26, 26, 26, 26)) %>%
  pivot_wider(
    id_cols = c(response, datasetID),
    names_from = c(".metric"),
    values_from = c(".estimate")) %>%
  full_join(tuned_res_df) %>%
  mutate(
    cv = "3x10folds",
    tuning_method = "bayesian",
    model_type = "ranger")


# check:
kable(mod_res)

write.csv2(mod_res, here::here("Data/output/B_02_ranger_atlas_separate.csv"))


# Interactions with H-stats

# Interactions with H-stats
library(hstats)

interaction_res <- replicate(length(results_atlas_i), list())

for (atlas_i in seq_along(results_atlas_i)){
  interaction_res[[atlas_i]] <- list()
  this_atlas <- results_atlas_i[[atlas_i]]

  for (response_i in seq_along(this_atlas)){

    # prepare data for hstats function
  dd <- this_atlas[[response_i]]
  fit <- dd$final_fit %>%
    extract_fit_parsnip()

  preds <- dd$explainer$data

  # run hstats function
  interactions <- hstats(fit,
                         X = preds,
                         threeway_m = 5)

  # save results
  interaction_res[[atlas_i]][[response_i]] <- interactions

  library(ggthemes)
  library(here)
  library(ggplot2)

  # create the plot
  interaction_plot <- plot(interactions, which = 1:100) + theme_few()

  if (atlas_i == 1) datasetID <- 5 else
    if (atlas_i == 2) datasetID <- 6 else
      if (atlas_i == 3) datasetID <- 13 else
        if (atlas_i == 4) datasetID <- 26

  # save as .svg
  ggsave(filename = here::here("Figures/B_models/interactions/", paste0("B_02_",responses[response_i],"_", datasetID, "_interactions_hstats.svg")),
         plot = interaction_plot,
         device = "svg",
         width = 10, height = 10)


}
}

interaction_res
saveRDS(interaction_res, here::here("Data/output/temp/B_02_Interaction_res_atlas_split.rds"))




# Save interactions in numbers (csv) for Documentation:
#Source 00_Configuration.R
source(here::here("R/00_Configuration.R"))
package_list <- c(package_list, "ranger", "tidymodels", "caret", "skimr", "DALEXtra")
lapply(package_list, require, character = TRUE)
tidymodels_prefer()
library(hstats)
interaction_res <-  readRDS(here::here("Data/output/temp/B_02_Interaction_res_atlas_split.rds"))


res_list <- replicate(4, list)
for (atlas_i in seq_along(interaction_res)){
  this_atlas <- interaction_res[[atlas_i]]
  res_list[[atlas_i]] <- list()

  if (atlas_i == 1) datasetID <- 5 else
    if (atlas_i == 2) datasetID <- 6 else
      if (atlas_i == 3) datasetID <- 13 else
        if (atlas_i == 4) datasetID <- 26

        print(datasetID)

  for (resp_i in seq_along(this_atlas)){

    if(resp_i == 1) this_response <- "Jaccard_dissim" else
      if(resp_i == 2) this_response <- "log_R2_1" else
        if(resp_i == 3) this_response <- "log_R2_1_per_year"

        print(this_response)

        response_res <- this_atlas[[resp_i]]

        H2_all_df <-
          data.frame(
            row.names = NULL,
            datasetID = datasetID,
            response = this_response,
            variable = "all",
            test = "total_interaction_strength",
            H2 = summary(response_res)$h2[[1]])

        H2_overall_df <-
          data.frame(
            row.names = NULL,
            datasetID = datasetID,
            response = this_response,
            variable = row.names(summary(response_res)$h2_overall[[1]]),
            test = "overall",
            H2 = summary(response_res)$h2_overall[[1]])

        H2_pairwise_df <-
          data.frame(
            row.names = NULL,
            datasetID = datasetID,
            response = this_response,
            variable = row.names(summary(response_res)$h2_pairwise[[1]]),
            test = "pairwise",
            H2 = summary(response_res)$h2_pairwise[[1]]
          )

        H2_threeway_df <-
          data.frame(
            row.names = NULL,
            datasetID = datasetID,
            response = this_response,
            variable = row.names(summary(response_res)$h2_threeway[[1]]),
            test = "threeway",
            H2 = summary(response_res)$h2_threeway[[1]]
          )


        temp <-
          full_join(H2_all_df, H2_overall_df) %>%
          full_join(H2_pairwise_df) %>%
          full_join(H2_threeway_df)

        res_list[[atlas_i]][[resp_i]] <- temp


  }

}


res_df <- bind_rows(res_list)

write.csv2(res_df, here::here("Data/output/results/B_02_Hstats_split_data.csv"))








# Variable importances:
res_all <- readRDS(here::here("Data/output/B_models/B_02_rf_res_atlas_split.rds"))
responses <- c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year")

vimp_list <- replicate(4, list)
for(atlas_i in seq_along(res_all)){
  vimp_list[[atlas_i]] <- list()
  this_atlas <- res_all[[atlas_i]]

  if (atlas_i == 1) datasetID <- 5 else
    if (atlas_i == 2) datasetID <- 6 else
      if (atlas_i == 3) datasetID <- 13 else
        if (atlas_i == 4) datasetID <- 26

for(i in seq_along(this_atlas)){

  if(i == 1) this_response <- "Jaccard_dissim" else
    if(i == 2) this_response <- "log_R2_1" else
      if(i == 3) this_response <- "log_R2_1_per_year"

      fit <- this_atlas[[i]]$final_fit %>% extract_fit_engine()
      vimp <- fit$variable.importance

      vimp_list[[atlas_i]][[i]] <- data.frame(datasetID = datasetID,
                                   variable = names(vimp),
                                   importance = vimp,
                                   row.names = NULL,
                                   response = this_response,
                                   importance_mode = fit$importance.mode
      )
}
}

imp_df <-
  vimp_list %>%
  bind_rows()

write.csv2(imp_df, here::here("Data/output/results/B_02_ranger_vimp_split_data.csv"))
