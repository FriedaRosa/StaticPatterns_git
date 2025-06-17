## Zero change species ##
library(dplyr)
library(here)
here::here()

zero1_2 <- readRDS(here::here("Data/output/1_all_predictors_merged.rds")) %>%
  filter(log_R2_1 >= log(0.9) & log_R2_1 <= log(1.1))
zero3 <- readRDS(here::here("Data/output/2_big_table_3.rds")) %>%
  filter(log_R3_2 >= log(0.9) & log_R3_2 <= log(1.1)) %>%
  mutate(
    datasetID = as.factor(as.character(datasetID)),
    samplingPeriodID = as.factor(as.character(samplingPeriodID))
  )

library(ggplot2)
library(ggthemes)

log_ratios <- zero1_2 %>%
  mutate(log_ratio = log_R2_1) %>%
  full_join(zero3 %>% mutate(log_ratio = log_R3_2)) %>%
  select(-log_R2_1, -log_R3_2)


p <- ggplot(log_ratios) +
  aes(x = log_ratio, y = Jaccard_dissim, fill = datasetID, color = datasetID) +
  geom_smooth(method = "gam", se = TRUE, alpha = 0.2) +
  scale_fill_manual(values = c(
    "5" = "#E66101",
    "6" = "#FDB863",
    "13" = "#B2ABD2",
    "26" = "#5E3C99"
  )) +
  scale_color_manual(values = c(
    "5" = "#E66101",
    "6" = "#FDB863",
    "13" = "#B2ABD2",
    "26" = "#5E3C99"
  )) +
  labs(
    x = "log ratio AOO",
    y = "J",
    title = "Zero-change species",
    subtitle = "Jaccard~log ratio",
    fill = "Dataset",
    color = "Dataset"
  ) +
  ggthemes::theme_par() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(face = "italic"),
    axis.title.x = element_text(face = "italic")
  ) +
  facet_wrap(vars(samplingPeriodID))



ggplot(log_ratios) +
  aes(x = log_ratio, y = Jaccard_dissim, fill = samplingPeriodID, color = samplingPeriodID) +
  geom_smooth(method = "gam", se = TRUE, alpha = 0.2) +
  scale_fill_manual(values = c(
    "1" = "#E66101",
    "2" = "#FDB863",
    "3" = "#B2ABD2"
    )) +
  scale_color_manual(values = c(
    "1" = "#E66101",
    "2" = "#FDB863",
    "3" = "#B2ABD2"
  )) +
  labs(
    x = "log ratio AOO",
    y = "J",
    title = "Zero-change species",
    #subtitle = "Jaccard~log ratio",
    fill = "Sampling Period",
    color = "Sampling Period"
  ) +
  ggthemes::theme_par() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(face = "italic"),
    axis.title.x = element_text(face = "italic")
  ) +
  facet_wrap(vars(datasetID),
             labeller = labeller(datasetID = c("5" = "Czechia",
                                               "6" = "New York",
                                               "13" = "Japan",
                                               "26" = "Europe")))


log_ratios %>% filter(log_ratio == 0) %>% group_by(datasetID, samplingPeriodID) %>%
  summarize(n_sp = n_distinct(verbatimIdentification),
            n_sp_BL = n_distinct(scientificName))


library(rstatix)

zero1_2 %>% get_summary_stats(type = "common")



get_mode_from_hist <- function(x, bins = 100) {
  h <- hist(x, breaks = bins, plot = TRUE)
  mode_bin_center <- h$mids[which.max(h$counts)]
  return(mode_bin_center)
}

# Example usage
mode_estimate <- get_mode_from_hist(log_ratios$log_ratio)
