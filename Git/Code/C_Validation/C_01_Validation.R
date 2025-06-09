#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                   C_01_Validation.R
#                
#
#                    Friederike Wölke 
#                           2025
#
#----------------------------------------------------------#

library(here)
library(tidymodels)


#----------------------------------------------------------#

res_full <-
  readRDS(here::here("Data/output/B_models/B_01_list_all_results_rf.rds"))

res_split <-
  readRDS(here::here("Data/output/B_models/B_02_rf_res_atlas_split.rds"))

tp3 <- readRDS(here::here("Data/output/2_big_table_3.rds")) %>%
  select(datasetID, verbatimIdentification, Jaccard_dissim, log_R3_2) %>%
  mutate(datasetID = as.factor(datasetID))

dta2 <- readRDS(here::here("Data/output/1_all_predictors_merged.rds")) %>%
  filter(samplingPeriodID == 2) %>%
  filter(datasetID %in% c(5,13))

#----------------------------------------------------------#

Jaccard_validation <- res_full$res[[1]] %>%
  extract_fit_parsnip() %>%
  predict(new_data = dta2) %>%
  bind_cols(dta2 %>% select(verbatimIdentification, datasetID)) %>%
  left_join(tp3) %>%
  mutate(diff = Jaccard_dissim - .pred)

logRatio_validation <- res_full$res[[2]] %>%
  extract_fit_parsnip() %>%
  predict(new_data = dta2) %>%
  bind_cols(dta2 %>% select(verbatimIdentification, datasetID)) %>%
  left_join(tp3) %>%
  mutate(diff = log_R3_2 - .pred)


#----------------------------------------------------------#
# Plots:
#----------------------------------------------------------#

hist(logRatio_validation$diff)
hist(Jaccard_validation$diff)

Jaccard_validation %>%
ggplot()+
  geom_histogram(aes(x = diff, fill = datasetID), alpha = 0.5, position = "identity")+
  ggthemes::theme_par()

logRatio_validation %>%
  ggplot()+
  geom_histogram(aes(x = diff, fill = datasetID), alpha = 0.5, position = "identity")+
  ggthemes::theme_par()


p_1 <- Jaccard_validation %>%
  ggplot(aes(x = .pred, y = Jaccard_dissim, col = datasetID))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(alpha = 0.5)+
  geom_smooth(aes(col = datasetID), method= "lm")+
  ggthemes::theme_few()+
  scale_color_manual(values = c("#f1a340", "#998ec3"))

p_2 <- logRatio_validation %>%
  ggplot(aes(x = .pred, y = log_R3_2, col = datasetID))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(alpha = 0.5)+
  xlim(-4,2)+
  #geom_smooth()+
  ggthemes::theme_few()+
  scale_color_manual(values = c("#f1a340", "#998ec3"))

p_1
p_2
ggsave("Figures/C_validation/C_01_Validation_Jaccard.svg", p_1, width = 5, height = 4)
ggsave("Figures/C_validation/C_01_Validation_logRatio.svg", p_2, width = 5, height = 4)

esquisse::ggplot_to_ppt() # saved to here("Figures/C_validation/Prediction_against_replication_3_plots.pptx")


#----------------------------------------------------------#
# Stats #
#----------------------------------------------------------#

library(yardstick)

# For Jaccard dissimilarity model
metrics_jaccard <- Jaccard_validation %>% group_by(datasetID) %>%
  metrics(truth = Jaccard_dissim, estimate = .pred)

# For log ratio model
metrics_logRatio <- logRatio_validation %>% group_by(datasetID) %>%
  metrics(truth = log_R3_2, estimate = .pred)

print(metrics_jaccard)
print(metrics_logRatio)


# Pearson correlation for Jaccard dissimilarity
cor_jaccard <- cor.test(Jaccard_validation$.pred, Jaccard_validation$Jaccard_dissim, method = "pearson")

# Pearson correlation for log ratio
cor_logRatio <- cor.test(logRatio_validation$.pred, logRatio_validation$log_R3_2, method = "spearman")
cor.test(logRatio_validation$.pred, logRatio_validation$log_R3_2, method = "pearson")


print(cor_jaccard)
print(cor_logRatio)




# correctly predicted species # 
logRatio_validation %>% filter(abs(diff) < 0.01) %>% distinct() %>% skimr::skim()


# jaccard histograms for each region #

Jaccard_validation %>%
  filter(Jaccard_dissim > 0.75) %>% distinct() %>% arrange(desc(Jaccard_dissim))
  group_by(datasetID) %>% 
  ggplot()+
  geom_histogram(aes(Jaccard_dissim))+
  facet_wrap(~ datasetID)

  
  
  ####
  ####
dta <- readRDS(here("Data/output/1_all_predictors_merged.rds"))

dta %>% 
  filter(Jaccard_dissim > 0.9) %>% 
  filter(samplingPeriodID == 1)%>%
  select(Jaccard_dissim, log_R2_1, a,b,c,d, rel_AOO, datasetID, verbatimIdentification) %>%
  distinct() %>% 
  arrange(desc(Jaccard_dissim)) %>% 
  group_by(datasetID) %>% 
  skimr::skim()



dta %>%
  filter(samplingPeriodID == 1) %>%
  group_by(datasetID) %>%
  filter(rel_AOO >= 0.8) %>%
  select(Jaccard_dissim, log_R2_1, a,b,c,d, rel_AOO, datasetID, verbatimIdentification) %>%
  filter(log_R2_1 != 0) %>%
  distinct() %>% 
  arrange(desc(Jaccard_dissim)) %>% 
  group_by(datasetID) %>% 
  skimr::skim()


dta %>%
  filter(samplingPeriodID == 1) %>%
  group_by(datasetID) %>%
  select(verbatimIdentification, datasetID, rel_AOO) %>%
  slice_max(rel_AOO)

dta %>% filter(datasetID == 5) %>%
  slice_max(abs(log_R2_1))

dta %>% 
  filter(samplingPeriodID == 1) %>%
  group_by(datasetID) %>%
  mutate(abs_logR = abs(log_R2_1)) %>%
  filter(abs_logR <= 0.1 ) %>%
  select(verbatimIdentification, Jaccard_dissim, abs_logR, log_R2_1, datasetID, rel_AOO) %>% # %>% skimr::skim() 
filter(Jaccard_dissim != 0)

dta %>% 
  filter(samplingPeriodID == 1) %>%
  group_by(datasetID) %>%
  mutate(abs_logR = abs(log_R2_1)) %>%
  filter(abs_logR <= 0.1 ) %>%
  select(verbatimIdentification, Jaccard_dissim, abs_logR, log_R2_1, datasetID, rel_AOO) %>%
ggplot(aes(y = Jaccard_dissim, x= datasetID))+
  geom_boxplot()


dta %>% filter(samplingPeriodID == 1) %>%
  mutate(abs_logR = abs(log_R2_1),
         abs_logR_y = abs(log_R2_1_per_year)) %>%
  ggplot(aes(color = datasetID, x = abs_logR_y, y = Jaccard_dissim))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "gam")+
  xlim(0,0.3)+
  geom_vline(xintercept =0)+
  ylim(0,1)+
 facet_wrap(~datasetID, scales = "free_x")+
  ggthemes::theme_par()

dta %>% 
  filter(samplingPeriodID == 1) %>%
  group_by(datasetID) %>%
  mutate(abs_logR = abs(log_R2_1),
         abs_logR_y = abs(log_R2_1_per_year)) %>%
  filter(abs_logR_y <= 0.001 ) %>%
  select(verbatimIdentification, abs_logR_y, Jaccard_dissim, abs_logR, log_R2_1, datasetID, rel_AOO) %>%
  filter(Jaccard_dissim <=0.1)
