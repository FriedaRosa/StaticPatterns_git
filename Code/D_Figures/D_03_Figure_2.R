#----------------------------------------------------#
# Load libraries
#----------------------------------------------------#

library(dplyr)
library(ggplot2)

# Source 00_Configuration.R
source(here::here("Code/00_Configuration.R"))
package_list <- c(package_list, "ranger", "tidymodels", "caret", "skimr", "DALEXtra")
lapply(package_list, require, character = TRUE)
tidymodels_prefer()

#----------------------------------------------------#
# Start clean
#----------------------------------------------------#

rm(list = ls())
gc()
set.seed(123)

#---------------------------------------------------#
# Plot settings
#---------------------------------------------------#

# Color vector for Log Ratio
my_cols_LR <- c(
  "strong decrease (> halfing)" = "#d7191c",
  "weak decrease (< halfing)" = "#fdae61",
  "stable" = "#e0e0e0",
  "weak increase (< doubling)" = "#a6d96a",
  "strong increase (> doubling)" = "#1a9641"
)

# Color vector for Jaccard
my_cols_J <- c(
  "stable (< 0.1)" = "#B35806",
  "weak turnover (0.1 - 0.25)" = "#F1A340",
  "weak intermediate turnover (0.25 - 0.5)" = "#FEE0B6",
  "strong intermediate turnover (0.5 - 0.75)" = "#D8DAEB",
  "strong turnover (> 0.75)" = "#998EC3",
  "complete turnover (> 0.9)" = "#542788"
)

# Standard Theme
common_theme <- theme_classic() +
  theme(
    axis.line = element_line(colour = "azure4", linetype = "solid"),
    axis.ticks = element_line(colour = "gray13"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(vjust = 0.9, hjust = 0.9),
    # plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 2),
    legend.position = "none"
  )



#---------------------------------------------------#
#
# Part 1: Get model results
#
#---------------------------------------------------#


# Set variables for model
sp_id <- c("verbatimIdentification", "scientificName")
H1 <- c("Mass", "GlobRangeSize_km2", "Migration", "Habitat_5", "Generalism", "Threatened", "pd")
H2 <- c("D_AOO_a", "mean_lnLac", "AOO", "joincount_delta", "circNorm", "minDist_toBorder_centr")
H3 <- c("datasetID")
predictors <- c(H1, H2, H3)
responses <- c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year")

# filter data to necessary columns
dat <-
  readRDS(here::here("Data/output/1_all_predictors_merged.rds")) %>%
  select(all_of(c(sp_id, responses, H3, H1, H2, "samplingPeriodID"))) %>%
  ungroup()

# Load models
all_models <-
  readRDS(here::here("Data/output/B_models/B_01_list_all_results_rf.rds"))

fit_J <- all_models$res[[1]] %>%
  extract_fit_parsnip()

fit_lnRR <- all_models$res[[2]] %>%
  extract_fit_parsnip()

fit_lnRR_y <- all_models$res[[3]] %>%
  extract_fit_parsnip()

# Predictions from the models:
comp_df <- dat %>%
  filter(samplingPeriodID == 2) %>%
  select(
    -samplingPeriodID,
    -Jaccard_dissim,
    -log_R2_1,
    -log_R2_1_per_year,
    #-verbatimIdentification,
    #-scientificName
  )

comp_df$Jacc_pred <-
  round(
    stats::predict(
      all_models$res[[1]] %>%
        tune::extract_workflow(),
      new_data = comp_df
    ),
    4
  )$.pred

comp_df$LogRatio_pred <-
  round(
    predict(
      all_models$res[[2]] %>%
        tune::extract_workflow(),
      new_data = comp_df
    ),
    4
  )$.pred

comp_df$LogRatio_y_pred <-
  round(
    predict(
      all_models$res[[3]] %>%
        tune::extract_workflow(),
      new_data = comp_df
    ),
    4
  )$.pred

# checks:
comp_df %>%
  select(verbatimIdentification, datasetID, Jacc_pred, LogRatio_y_pred, LogRatio_pred) %>%
  head()


# Merge predictions to the data
dat_mod <- full_join(comp_df, dat) %>%
  mutate(
    Jaccard = case_when(
      samplingPeriodID == 2 ~ Jacc_pred,
      .default = Jaccard_dissim
    ),
    log_R2_1 = case_when(
      samplingPeriodID == 2 ~ LogRatio_pred,
      .default = log_R2_1
    ),
    log_R2_1_y = case_when(
      samplingPeriodID == 2 ~ LogRatio_y_pred,
      .default = log_R2_1_per_year
    )
  ) %>%
  select(datasetID, samplingPeriodID, verbatimIdentification, log_R2_1, log_R2_1_y, Jaccard)

# Checks:
dat_mod %>%
  filter(samplingPeriodID == 1 & Jaccard != 0) %>%
  glimpse()



#
# Part 2: Make trends categorical
#

# set factors and cutoff values for categories:

# Define factor levels
factor_levels_LR <- c(
  "strong decrease (> halfing)",
  "weak decrease (< halfing)",
  "stable",
  "weak increase (< doubling)",
  "strong increase (> doubling)"
)
# Define the cutoff points for the ranges in log scale
cutoff_points_LR <- c(-Inf, log(0.5), log(0.9), log(1.1), log(2), Inf)

# Define factor levels
factor_levels_J <- c(
  "stable (< 0.1)",
  "weak turnover (0.1 - 0.25)",
  "weak intermediate turnover (0.25 - 0.5)",
  "strong intermediate turnover (0.5 - 0.75)",
  "strong turnover (> 0.75)",
  "complete turnover (> 0.9)"
)
# Define the cutoff points for the ranges in log scale
cutoff_points_J <- c(-Inf, 0.1, 0.25, 0.5, 0.75, 0.9, Inf)


# make copy from data:
dat_cat <- dat

# prepare data for plotting:


# Convert the log ratios and Jaccard to factor variables based on the ranges
dat_cat$trend_LR <-
  cut(dat_cat$log_R2_1,
    breaks = cutoff_points_LR,
    labels = factor_levels_LR
  )

dat_cat$trend_J <-
  cut(dat_cat$Jaccard,
    breaks = cutoff_points_J,
    labels = factor_levels_J
  )

dat_cat <- dat_cat %>%
  rstatix::reorder_levels(trend_LR, factor_levels_LR) %>%
  rstatix::reorder_levels(trend_J, factor_levels_J) %>%
  na.omit()

dat_cat_wide <- dat_cat %>%
  tidyr::pivot_wider(
    id_cols = c(datasetID, verbatimIdentification),
    names_from = samplingPeriodID,
    values_from = c(log_R2_1, trend_LR, Jaccard_dissim, trend_J)
  )

print(colnames(dat_cat_wide))
lookup <- c(
  "LR_Observed" = "log_R2_1_1",
  "trend_LR_Observed" = "trend_LR_1",
  "LR_Predicted" = "log_R2_1_2",
  "trend_LR_Predicted" = "trend_LR_2",
  "J_Observed" = "Jaccard_dissim_1",
  "trend_J_Observed" = "trend_J_1",
  "J_Predicted" = "Jaccard_dissim_2",
  "trend_J_Predicted" = "trend_J_2"
)
dat_cat_wide <- dat_cat_wide %>%
  dplyr::rename(all_of(lookup)) %>%
  arrange(datasetID)



names_data <- unique(dat_cat_wide$datasetID)

# Reorder columns
dat_cat_wide <- dat_cat_wide[, c(
  "datasetID",
  "verbatimIdentification",
  "J_Predicted",
  "trend_J_Predicted",
  "J_Observed",
  "trend_J_Observed",
  "LR_Predicted",
  "trend_LR_Predicted",
  "LR_Observed",
  "trend_LR_Observed"
)]

# Plot: ----

library(grid)
J_LR_plot <- ggplot(
  dat_cat %>% filter(samplingPeriodID == 1),
  aes(y = Jaccard_dissim, x = log_R2_1, fill = trend_LR)
) +
  ## Prepare background colors ============================================== ##
  geom_hline(yintercept = 0.25, lwd = 0.2, col = "darkgrey") +
  geom_hline(yintercept = 0.5, lwd = 0.1, col = "red") +
  geom_hline(yintercept = 0.75, lwd = 0.2, col = "darkgrey") +
  annotation_custom(
    grob = rectGrob(gp = gpar(fill = "#fdae61", col = "#fdae61", alpha = 0.3)), xmin = log(0.5), xmax = log(0.9), ymin = 1, ymax = 0.1
  ) + # top-left
  annotation_custom(
    grob = rectGrob(gp = gpar(fill = "#a6d96a", col = "#a6d96a", alpha = 0.3)), xmin = log(1.1), xmax = log(2), ymin = 1, ymax = 0.1
  ) + # top-right
  annotation_custom(
    grob = rectGrob(gp = gpar(fill = "#d7191c", col = "#d7191c", alpha = 0.3)), xmin = -2, xmax = log(0.5), ymin = 1, ymax = 0.5
  ) + # middle-left
  annotation_custom(
    grob = rectGrob(gp = gpar(fill = "#1a9641", col = "#1a9641", alpha = 0.3)),
    xmin = log(2), xmax = 2, ymin = 1, ymax = 0.5
  ) + # bottom-right
  annotation_custom(
    grob = rectGrob(gp = gpar(fill = "#d7191c", col = "#d7191c", alpha = 0.5)),
    xmin = -Inf, xmax = -2, ymin = 1, ymax = 0.75
  ) + # bottom-left
  annotation_custom(
    grob = rectGrob(gp = gpar(fill = "#1a9641", col = "#1a9641", alpha = 0.5)),
    xmin = 2, xmax = Inf, ymin = 1, ymax = 0.75
  ) + # middle-right
  annotation_custom(
    grob = rectGrob(gp = gpar(fill = "gray", col = "gray", alpha = 0.3)),
    xmin = log(0.9), xmax = log(1.1), ymin = 0, ymax = 1
  ) + # bottom-right

  ## Add points ============================================================= ##
  geom_jitter(shape = 21, col = "darkgrey", show.legend = F) +

  ## Theme ================================================================== ##
  ggthemes::theme_few() +



  ## Scales ================================================================= ##
  # scale_x_continuous(n.breaks = 9, limits = c(-4, 4)) +
  scale_fill_manual(values = my_cols_LR) +
  labs(
    y = "Jaccard dissimilarity (Observed)",
    x = "Log Ratio (Observed)",
    fill = "Direction"
  )


J_LR_plot2 <- ggExtra::ggMarginal(
  p = J_LR_plot,
  type = "densigram",
  margins = "both",
  size = 5,
  colour = "gray34",
  fill = "#EBEBEBC5"
)

J_LR_plot2
export::graph2ppt(J_LR_plot2, "Figures/A_data/D_03_LogRatio_Jaccard_Marginal.ppt")
ggsave("Figures/A_data/D_03_LogRatio_Jaccard_Marginal.pdf", J_LR_plot2, width = 5, height = 5)
ggsave("Figures/A_data/D_03_LogRatio_Jaccard_Marginal.png", J_LR_plot2, width = 5, height = 5)

# esquisse::ggplot_to_ppt("J_LR_plot2")
#
## Legend:

legend <- ggplot(
  dat_cat %>% filter(samplingPeriodID == 1),
  aes(y = Jaccard_dissim, x = log_R2_1, fill = trend_LR)
) +
  geom_jitter(shape = 21, col = "darkgrey", show.legend = T) +

  ## Theme ================================================================== ##
  ggthemes::theme_few() +
  scale_fill_manual(values = my_cols_LR) +
  labs(
    y = "Jaccard dissimilarity (Observed)",
    x = "Log Ratio (Observed)",
    fill = "Direction"
  ) +
  theme(legend.position = "right")
ggsave("Figures/A_data/D_03_LogRatio_Jaccard_Marginal_legend.svg", legend, width = 10, height = 5)
