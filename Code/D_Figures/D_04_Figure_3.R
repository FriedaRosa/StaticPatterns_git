#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                      D_04_Figure_3.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#

source(here::here("Code/00_Configuration.R"))
lapply(package_list, require, character = TRUE)
# install.packages("tidytext")
library(tidytext)

# Start with clean environment
rm(list = ls())
gc()

#----------------------------------------------------------#

# Variables:
sp_id <- c("verbatimIdentification", "scientificName")
H1 <- c("Mass", "GlobRangeSize_km2", "Migration", "Habitat_5", "Generalism", "Threatened", "pd")
H2 <- c("D_AOO_a", "mean_lnLac", "AOO", "joincount_delta", "circNorm", "minDist_toBorder_centr")
H3 <- c("datasetID")
predictors <- c(H1, H2, H3)
responses <- c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year")

#----------------------------------------------------------#

mp_full <- read.csv2(here::here("Data/output/results/B_01_ranger_all_data.csv"))
mp_split <- read.csv2(here::here("Data/output/results/B_02_ranger_atlas_separate.csv"))

vip_full <- read.csv2(here::here("Data/output/results/B_01_ranger_vimp_all_data.csv"))
vip_split <- read.csv2(here::here("Data/output/results/B_02_ranger_vimp_split_data.csv"))

int_full <- read.csv2(here::here("Data/output/results/B_01_Hstats_all_data.csv"))
int_split <- read.csv2(here::here("Data/output/results/B_02_Hstats_split_data.csv"))
#----------------------------------------------------------#



#----------------------------------------------------------#
# a) model performances (rsq, rmse)
#----------------------------------------------------------#

p_1 <-
  mp_full %>%
  pivot_longer(cols = c(rmse, rsq)) %>%
  ggplot() +
  geom_point(aes(x = response, y = value, shape = name), cex = 2, col = "black", bg = "black") +
  ggthemes::theme_few() +
  scale_shape_manual(values = c(rmse = 21, rsq = 1))


p_1

p_2 <-
  mp_split %>%
  pivot_longer(cols = c(rmse, rsq)) %>%
  ggplot() +
  geom_point(aes(x = factor(datasetID), y = value, shape = name), cex = 2, col = "black", bg = "black") +
  facet_wrap(~response) +
  ggthemes::theme_few() +
  scale_shape_manual(values = c(rmse = 21, rsq = 1)) +
  xlab("DatasetID")


p_2

ggsave(filename = "Model_performance_all_data.svg", plot = p_1, path = here::here("Figures/B_models/performance"), device = "svg", width = 5.5, height = 5)
ggsave(filename = "Model_performance_split_data.svg", plot = p_2, path = here::here("Figures/B_models/performance"), device = "svg", width = 5.5, height = 5)


#----------------------------------------------------------#
# b) Variable importances:
#----------------------------------------------------------#

vip_full <- vip_full %>%
  group_by(variable) %>%
  mutate(mean_importance_jaccard = mean(importance[response == "Jaccard_dissim"], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(variable = fct_reorder(variable, mean_importance_jaccard)) %>%
  mutate(hypo = case_when(
    variable %in% H1 ~ "H1",
    variable %in% H2 ~ "H2",
    variable %in% H3 ~ "H3"
  ))

p_3 <-
  vip_full %>%
  ggplot() +
  geom_col(aes(x = importance, y = variable, fill = hypo)) +
  facet_wrap(~response, scales = "free_x") + # free_x so y-axis stays fixed
  labs(y = "Variable", x = "Importance") +
  ggthemes::theme_few() +
  scale_fill_manual(values = c("#f1a340", "#7fbf7b", "#998ec3"))


p_3


vip_split <- vip_split %>%
  group_by(variable, datasetID, response) %>%
  mutate(mean_importance_jaccard = mean(importance[response == "Jaccard_dissim"], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(variable = fct_reorder(variable, mean_importance_jaccard)) %>%
  mutate(hypo = case_when(
    variable %in% H1 ~ "H1",
    variable %in% H2 ~ "H2",
    variable %in% H3 ~ "H3"
  ))

p_4 <-
  vip_split %>%
  ggplot() +
  geom_col(aes(x = importance, y = variable, fill = hypo)) +
  facet_wrap(~ datasetID + response, scales = "free_x", ncol = 3) + # free_x so y-axis stays fixed
  labs(y = "Variable", x = "Importance") +
  ggthemes::theme_few() +
  scale_fill_manual(values = c("#f1a340", "#7fbf7b", "#998ec3"))



p_4 # saved as 1200x800
ggsave(filename = "Variable_importance_all_data.svg", plot = p_3, path = here::here("Figures/B_models/performance"), device = "svg", width = 11, height = 5)
ggsave(filename = "Variable_importance_split_data.svg", plot = p_4, path = here::here("Figures/B_models/performance"), device = "svg", width = 11, height = 15)



#----------------------------------------------------------#
# split by response into 3 plots
#----------------------------------------------------------#

P_4a <- vip_split %>%
  filter(response == "Jaccard_dissim") %>%
  group_by(variable, datasetID, response) %>%
  mutate(mean_importance_jaccard = mean(importance[response == "Jaccard_dissim" && datasetID == 5], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(variable = fct_reorder(variable, mean_importance_jaccard)) %>%
  mutate(hypo = case_when(
    variable %in% H1 ~ "H1",
    variable %in% H2 ~ "H2",
    variable %in% H3 ~ "H3"
  )) %>%
  ggplot() +
  geom_col(aes(x = importance, y = variable, fill = hypo)) +
  facet_wrap(~datasetID, scales = "free_x", ncol = 2) + # free_x so y-axis stays fixed
  labs(y = "Variable", x = "Importance") +
  ggthemes::theme_few() +
  scale_fill_manual(values = c("#f1a340", "#7fbf7b", "#998ec3")) +
  ggtitle("Jaccard")
P_4a


P_4b <- vip_split %>%
  filter(response == "log_R2_1") %>%
  group_by(variable, datasetID, response) %>%
  mutate(mean_importance_lnRR = mean(importance[response == "log_R2_1" && datasetID == 5], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(variable = fct_reorder(variable, mean_importance_lnRR)) %>%
  mutate(hypo = case_when(
    variable %in% H1 ~ "H1",
    variable %in% H2 ~ "H2",
    variable %in% H3 ~ "H3"
  )) %>%
  ggplot() +
  geom_col(aes(x = importance, y = variable, fill = hypo)) +
  facet_wrap(~datasetID, scales = "free_x", ncol = 2) + # free_x so y-axis stays fixed
  labs(y = "Variable", x = "Importance") +
  ggthemes::theme_few() +
  scale_fill_manual(values = c("#f1a340", "#7fbf7b", "#998ec3")) +
  ggtitle("log ratio")


P_4b

P_4c <- vip_split %>%
  group_by(variable, datasetID, response) %>%
  mutate(mean_importance_lnRRy = mean(importance[response == "log_R2_1_per_year" && datasetID == 5], na.rm = TRUE)) %>%
  ungroup() %>%
  filter(response == "log_R2_1_per_year") %>%
  mutate(variable = fct_reorder(variable, mean_importance_lnRRy)) %>%
  mutate(hypo = case_when(
    variable %in% H1 ~ "H1",
    variable %in% H2 ~ "H2",
    variable %in% H3 ~ "H3"
  )) %>%
  ggplot() +
  geom_col(aes(x = importance, y = variable, fill = hypo)) +
  facet_wrap(~datasetID, scales = "free_x", ncol = 2) + # free_x so y-axis stays fixed
  labs(y = "Variable", x = "Importance") +
  ggthemes::theme_few() +
  scale_fill_manual(values = c("#f1a340", "#7fbf7b", "#998ec3")) +
  ggtitle("log ratio per year")

P_4c

ggsave(filename = "Figure_3_vip_J_split.svg", plot = P_4a, path = here::here("Figures/B_models/performance"), device = "svg", width = 10, height = 8)
ggsave(filename = "Figure_3_vip_lnR_split.svg", plot = P_4b, path = here::here("Figures/B_models/performance"), device = "svg", width = 10, height = 8)
ggsave(filename = "Figure_3_vip_lnRy_split.svg", plot = P_4c, path = here::here("Figures/B_models/performance"), device = "svg", width = 10, height = 8)


#----------------------------------------------------------#
# c) Partial dependencies:
#----------------------------------------------------------#

# Set variables:
H1 <- c("Mass", "GlobRangeSize_km2", "Migration", "Habitat_5", "Generalism", "Threatened", "pd")
H2 <- c("D_AOO_a", "mean_lnLac", "AOO", "joincount_delta", "circNorm", "minDist_toBorder_centr")
H3 <- c("datasetID")
predictors <- c(H1, H2, H3)
temp <- readRDS(here::here("Data/output/1_all_predictors_merged.rds"))

# Create lookup labels for clean plotting:
lookup <-
  c(
    "Mass" = "Body mass",
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
    "minDist_toBorder_centr" = "Smallest distance centroid-border"
  ) %>%
  as.data.frame()
lookup$curr_label <- row.names(lookup)
lookup$old_label <- lookup$.
row.names(lookup) <- NULL


# Categorical plot:
cat_vars <- c(
  "Migration",
  "Habitat_5",
  "Generalism",
  "Threatened",
  "datasetID"
)

# read back in
pdp_list <- readRDS(here::here("Data/output/temp/B_04_partial_dependencies.rds"))


# Merge to plotting df:
pdp_df_list <- replicate(15, list)
pdp_df_list_cat <- replicate(15, list)

for (model_i in seq_along(pdp_list)) {
  this_model <- pdp_list[[model_i]]
  pdp_df_list[[model_i]] <- list()
  pdp_df_list_cat[[model_i]] <- list()

  if (model_i %in% c(1, 4, 5, 6, 7)) response <- "Jaccard_dissim" else if (model_i %in% c(2, 8, 9, 10, 11)) response <- "log_R2_1" else if (model_i %in% c(3, 12, 13, 14, 15)) response <- "log_R2_1_per_year"

  if (model_i %in% c(1, 2, 3)) obj_ident <- "full" else if (model_i %in% c(4:7)) obj_ident <- "Jaccard_split" else if (model_i %in% c(8:11)) obj_ident <- "log_R2_1_split" else if (model_i %in% c(12:15)) obj_ident <- "log_R2_1_y_split"

  if (model_i %in% c(1, 2, 3)) datasetID <- "all" else if (model_i %in% c(4, 8, 12)) datasetID <- "5" else if (model_i %in% c(5, 9, 13)) datasetID <- "6" else if (model_i %in% c(6, 10, 14)) datasetID <- "13" else if (model_i %in% c(7, 11, 15)) datasetID <- "26"

  for (var_i in seq_along(this_model)) {
    pdp_df_i <- this_model[[var_i]]
    pdp_df_i$variable <- names(pdp_df_i)[1]

    pdp_df_i$value <- pdp_df_i[[1]]
    pdp_df_i$response <- response
    pdp_df_i$datasetID <- datasetID
    pdp_df_i$model <- obj_ident
    pdp_df_i[[1]] <- NULL

    if (unique(pdp_df_i$variable) %in% cat_vars) {
      pdp_df_list_cat[[model_i]][[var_i]] <- pdp_df_i
    } else {
      pdp_df_list[[model_i]][[var_i]] <- pdp_df_i
    }
  }
}

pdp_num <- pdp_df_list %>% bind_rows()
pdp_cat <- pdp_df_list_cat %>% bind_rows()

pdp_res <- list(
  numeric = pdp_num,
  categorical = pdp_cat
)



#--------------------------------------------------#
# Plots (ggplot2) ----
#--------------------------------------------------#

num_cat_i <- 1
model_i <- 1
resp_i <- 1

lookup <-
  c(
    "Mass" = "Body mass",
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
    "minDist_toBorder_centr" = "Smallest distance centroid-border"
  )




for (num_cat_i in seq_along(pdp_res)) { # there are 2 sub-lists (num, cat)

  if (num_cat_i == 1) { # numeric

    # chose sublist from results
    pdp_subset <- pdp_res[[num_cat_i]]

    # extract results from models dynamically
    for (model_i in seq_along(unique(pdp_subset$model))) { # there are 4 models

      this_model <- unique(pdp_subset$model)[model_i]
      model_dd <- pdp_subset %>% filter(model == this_model)

      for (resp_i in seq_along(unique(model_dd$response))) {
        this_resp <- unique(model_dd$response)[resp_i]

        print(paste0("Num/Cat = ", num_cat_i, "; model = ", this_model, "; response = ", this_resp))

        response_dd <- model_dd %>%
          filter(response == this_resp)

        p0 <- response_dd %>%
          ggplot() +
          geom_line(
            aes(
              y = yhat,
              x = value,
              col = datasetID,
              group = interaction(variable, datasetID)
            ),
            linewidth = 0.8,
            alpha = 0.8
          ) +
          labs(
            x = NULL,
            y = bquote(italic(.(this_resp)) ~ "(Partial)")
          ) +
          facet_wrap(~variable, scales = "free_x", labeller = labeller(variable = lookup)) +
          ggthemes::theme_par() +
          ggtitle(paste(this_resp), paste0(this_model, "-model")) +
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
            ylim(min(temp$log_R2_1), max(temp$log_R2_1)) +
            geom_hline(
              yintercept = 0,
              linewidth = 0.3,
              col = "darkgrey"
            )
        } else if (this_resp == "log_R2_1_per_year") {
          p1 <-
            p0 +
            ylim(min(temp$log_R2_1_per_year), max(temp$log_R2_1_per_year)) +
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
    for (model_i in seq_along(unique(pdp_subset$model))) { # there are 4 models

      this_model <- unique(pdp_subset$model)[model_i]
      model_dd <- pdp_subset %>% filter(model == this_model)

      for (resp_i in seq_along(unique(model_dd$response))) {
        this_resp <- unique(model_dd$response)[resp_i]

        print(paste0("Num/Cat = ", num_cat_i, "; model = ", this_model, "; response = ", this_resp))

        response_dd <- model_dd %>%
          filter(response == this_resp)

        p0 <- response_dd %>%
          ggplot() +
          geom_point(
            aes(
              y = yhat,
              x = value,
              col = datasetID
            ),
            alpha = 0.8
          ) +
          ggthemes::theme_par() +
          labs(y = bquote(italic(.(this_resp)) ~ "(Partial)")) +
          facet_wrap(~variable, scales = "free_x", labeller = labeller(variable = lookup)) +
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
            ylim(min(temp$log_R2_1), max(temp$log_R2_1)) +
            geom_hline(
              yintercept = 0,
              linewidth = 0.3,
              col = "darkgrey"
            )
        } else if (this_resp == "log_R2_1_per_year") {
          p1 <-
            p0 +
            ylim(min(temp$log_R2_1_per_year), max(temp$log_R2_1_per_year)) +
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
# Log Ratio per model:

plots_numeric <- list()
for (var_i in seq_along(unique(pdp_res[[1]]$variable))) {
  this_var <- unique(pdp_res[[1]]$variable)[var_i]
  this_resp <- "Jaccard_dissim"
  this_model <- "Jaccard_split"
  # Numeric Jaccard plots

  p0 <- pdp_res[[1]] %>%
    filter(model == "Jaccard_split" & variable == this_var) %>%
    ggplot() +
    geom_hline(yintercept = 0.5, linewidth = 0.3, col = "darkgrey") +
    geom_line(
      aes(
        x = value,
        y = yhat,
        group = datasetID,
        col = datasetID
      ),
      linewidth = 0.8
    ) +
    ggthemes::theme_few() +
    scale_color_manual(values = c(
      "all" = "black",
      "5" = "#E66101",
      "6" = "#FDB863",
      "13" = "#B2ABD2",
      "26" = "#5E3C99"
    )) +
    labs(
      x = lookup[this_var],
      y = bquote(italic(.(this_resp)) ~ "(Partial)"),
      col = "Dataset"
    ) +
    ggtitle(NULL, paste0(this_model, "-model")) +
    ylim(0, 1)

  p0
  plots_numeric[[var_i]] <- p0
  # ggsave(plot = last_plot(), path = here::here("Figures/B_models/performance"), filename = paste0("B_04_PartialPlot_Jaccard_split_",this_var,".svg" ))
}

plots_numeric

plots_categorical <- list()
for (var_i in seq_along(unique(pdp_res[[2]]$variable))) {
  this_var <- unique(pdp_res[[2]]$variable)[var_i]
  this_resp <- "Jaccard_dissim"
  this_model <- "Jaccard_split"
  # Categorical Jaccard plots

  p1 <- pdp_res[[2]] %>%
    filter(model == "Jaccard_split" & variable == this_var) %>%
    ggplot() +
    geom_hline(yintercept = 0.5, linewidth = 0.3, col = "darkgrey") +
    geom_point(
      aes(
        x = value,
        y = yhat,
        col = datasetID
      ),
      size = 3,
      alpha = 0.75
    ) +
    ggthemes::theme_few() +
    scale_color_brewer(palette = "PuOr") +
    labs(
      x = lookup[this_var],
      y = bquote(italic(.(this_resp)) ~ "(Partial)"),
      col = "Dataset"
    ) +
    ggtitle(NULL, paste0(this_model, "-model")) +
    ylim(0, 1)

  plots_categorical[[var_i]] <- p1
  # ggsave(plot = last_plot(), device = "svg", path = here::here("Figures/B_models/performance"), filename = paste0("B_04_PartialPlot_Jaccard_split_",this_var,".svg" ))
}

plots_categorical <- plots_categorical[c(1:4)]






patchwork::wrap_plots(plots_numeric, ncol = 3) +
  patchwork::plot_layout(tag_level = "keep", guides = "collect", axes = "collect_y") +
  patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")")

export::graph2ppt(
  file = here::here("Figures/B_models/performance/B_04_PartialSep_Jaccard.ppt"), width = 8, height = 10
)

patchwork::wrap_plots(plots_categorical, ncol = 3) +
  patchwork::plot_layout(tag_level = "keep", guides = "collect", axes = "collect_y") +
  patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")")

export::graph2ppt(
  file = here::here("Figures/B_models/performance/B_04_PartialSep_Categorical_Jaccard.ppt"), width = 8, height = (10 / 3) * 2
)





# Log Ratio per model:
plots_numeric <- list()
for (var_i in seq_along(unique(pdp_res[[1]]$variable))) {
  this_var <- unique(pdp_res[[1]]$variable)[var_i]
  this_resp <- "log_R2_1"
  this_model <- "log_R2_1_split"
  # Numeric Jaccard plots

  p0 <- pdp_res[[1]] %>%
    filter(model == "log_R2_1_split" & variable == this_var) %>%
    ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.3, col = "darkgrey") +
    geom_line(
      aes(
        x = value,
        y = yhat,
        group = datasetID,
        col = datasetID
      ),
      linewidth = 0.8
    ) +
    ggthemes::theme_few() +
    scale_color_manual(values = c(
      "all" = "black",
      "5" = "#E66101",
      "6" = "#FDB863",
      "13" = "#B2ABD2",
      "26" = "#5E3C99"
    )) +
    labs(
      x = lookup[this_var],
      y = bquote(italic(.(this_resp)) ~ "(Partial)"),
      col = "Dataset"
    ) +
    ggtitle(NULL, paste0(this_model, "-model")) +
    ylim(-1, 4)

  p0
  plots_numeric[[var_i]] <- p0
  # ggsave(plot = last_plot(), path = here::here("Figures/B_models/performance"), filename = paste0("B_04_PartialPlot_Jaccard_split_",this_var,".svg" ))
}

plots_numeric

plots_categorical <- list()
for (var_i in seq_along(unique(pdp_res[[2]]$variable))) {
  this_var <- unique(pdp_res[[2]]$variable)[var_i]
  this_resp <- "log_R2_1"
  this_model <- "log_R2_1_split"
  # Categorical Jaccard plots

  p1 <- pdp_res[[2]] %>%
    filter(model == "log_R2_1_split" & variable == this_var) %>%
    ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.3, col = "darkgrey") +
    geom_point(
      aes(
        x = value,
        y = yhat,
        col = datasetID
      ),
      size = 3,
      alpha = 0.75
    ) +
    ggthemes::theme_few() +
    scale_color_brewer(palette = "PuOr") +
    labs(
      x = lookup[this_var],
      y = bquote(italic(.(this_resp)) ~ "(Partial)"),
      col = "Dataset"
    ) +
    ggtitle(NULL, paste0(this_model, "-model")) +
    ylim(-1, 4)

  plots_categorical[[var_i]] <- p1
  # ggsave(plot = last_plot(), device = "svg", path = here::here("Figures/B_models/performance"), filename = paste0("B_04_PartialPlot_Jaccard_split_",this_var,".svg" ))
}

plots_categorical <- plots_categorical[c(1:4)]






patchwork::wrap_plots(plots_numeric, ncol = 3) +
  patchwork::plot_layout(tag_level = "keep", guides = "collect", axes = "collect_y") +
  patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")")

export::graph2ppt(
  file = here::here("Figures/B_models/performance/B_04_PartialSep_logR.ppt"), width = 8, height = 10
)

patchwork::wrap_plots(plots_categorical, ncol = 3) +
  patchwork::plot_layout(tag_level = "keep", guides = "collect", axes = "collect_y") +
  patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")")

export::graph2ppt(
  file = here::here("Figures/B_models/performance/B_04_PartialSep_Categorical_logR.ppt"), width = 8, height = (10 / 3) * 2
)





#----------------------------------------------------------#
# d) Interactions:
# #----------------------------------------------------------#


p_5 <-
  int_full %>%
  filter(!c(test == "total_interaction_strength")) %>%
  ggplot() +
  geom_col(aes(x = H2, y = reorder_within(variable, H2, response))) +
  facet_wrap(~ test + response, scales = "free") +
  scale_y_reordered() +
  ggthemes::theme_few()

p_5 # saved as 1200x800

p_6 <-
  int_split %>%
  filter(!c(test == "total_interaction_strength")) %>%
  ggplot() +
  geom_col(aes(x = H2, y = reorder_within(variable, H2, response), fill = factor(datasetID))) +
  facet_wrap(~ test + response, scales = "free") +
  scale_y_reordered() +
  ggthemes::theme_few()
p_6 # saved as 1500x1500


# save as ppt for postediting
esquisse::ggplot_to_ppt() # saved to folder here::here("Figures/B_models/performance/")
