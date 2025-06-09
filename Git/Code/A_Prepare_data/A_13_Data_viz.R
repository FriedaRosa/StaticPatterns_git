#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                 Data_vizualization.R
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
lapply(package_list, require, character = TRUE)

#----------------------------------------------------------#
# Read the data -----
#----------------------------------------------------------#

dta <- readRDS(here::here("Data/output/1_all_predictors_merged.rds")) %>%
  filter(samplingPeriodID == 1)


#----------------------------------------------------------#
# Set variables -----
#----------------------------------------------------------#


# Set variables for model
sp_id <- c("verbatimIdentification", "scientificName")
H1 <- c("Mass", "GlobRangeSize_km2", "Migration", "Habitat_5", "Generalism", "Threatened", "pd")
H2 <- c("D_AOO_a", "mean_lnLac", "AOO", "joincount_delta", "circNorm", "minDist_toBorder_centr")
H3 <- c("datasetID")
predictors <- c(H1, H2, H3)
responses <- c("Jaccard_dissim", "log_R2_1", "log_R2_1_per_year")

#----------------------------------------------------------#
# Reduce data to variables ----
#----------------------------------------------------------#


# modify data for modeling
dta_new <- dta %>%
  select(all_of(c(sp_id, responses, H3, H1, H2))) %>%
  ungroup()

#----------------------------------------------------------#
# Plots ----
#----------------------------------------------------------#

library(inspectdf)
library(summarytools)
library(caret)
library(ggplot2)
library(DataExplorer)
library(explore)

create_report(dta_new, output_dir =here::here("Figures/A_data/"), output_file = "A_13_report.html")
#explore(dta_new)
GGally::ggpairs(ggplot2::aes(colour = datasetID), data= dta_new %>% select(-verbatimIdentification, -scientificName)%>% as.data.frame())
inspect_types(dta_new) %>% show_plot()
inspect_mem(dta_new) %>% show_plot()
inspect_num(dta_new) %>% show_plot()
inspect_imb(dta_new) %>% show_plot()
inspect_cat(dta_new %>% select(-verbatimIdentification, -scientificName) %>% group_by(datasetID)) %>% show_plot()
inspect_cor(dta_new) %>% show_plot()


freq(dta_new)
ctable(x = dta_new$Habitat_5,
       y = dta_new$Threatened)

ctable(x = dta_new$datasetID,
       y = dta_new$Threatened)
descr(dta_new)
view(dfSummary(dta_new))


dfSummary(dta_new)
stby(data = dta_new,
     INDICES = dta_new$datasetID,
     FUN = descr,
     stats = "common",
     transpose = T)

# 1. Feature plot of relationships between all vars
#caret::featurePlot
trellis.par.set(theme = col.whitebg(), warn = FALSE)

# H1 ~ Jaccard
featurePlot(x = dta_new %>% select(datasetID, all_of(H1)),
            y =dta_new$Jaccard_dissim,
            group = dta_new$datasetID,
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "Scatterplot Matrix of traits (H1) - Jaccard 1",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))

# H2 ~ Jaccard
featurePlot(x = dta_new %>% select(datasetID, all_of(H2)),
            y =dta_new$Jaccard_dissim,
            group = dta_new$datasetID,
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "Scatterplot Matrix of range geometry (H2) - Jaccard 1",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))
#-------------------------------------------------------------#
# H1 ~ log Ratio
featurePlot(x = dta_new %>% select(datasetID, all_of(H1)),
            y =dta_new$log_R2_1,
            group = dta_new$datasetID,
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "Scatterplot Matrix of traits (H1) - Log Ratio",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))

# H2 ~ log Ratio
featurePlot(x = dta_new %>% select(datasetID, all_of(H2)),
            y =dta_new$log_R2_1,
            group = dta_new$datasetID,
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "Scatterplot Matrix of range geometry (H2) - log Ratio",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))
#--------------------------------------------------------------#

# Split range geometries (H2) for each atlas:
featurePlot(x = dta_new %>% filter(datasetID == 5) %>%
              select(all_of(c(H2))),
            y =  dta_new %>% filter(datasetID == 5) %>%
              select(Jaccard_dissim),
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "Czechia: Scatterplot Matrix of traits (H1) - Jaccard 1",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))

featurePlot(x = dta_new %>% filter(datasetID == 6) %>%
              select(all_of(c(H2))),
            y =  dta_new %>% filter(datasetID == 6) %>%
              select(Jaccard_dissim),
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "New York: Scatterplot Matrix of range geometry (H2) - Jaccard 1",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))

featurePlot(x = dta_new %>% filter(datasetID == 13) %>%
              select(all_of(c(H2))),
            y =  dta_new %>% filter(datasetID == 13) %>%
              select(Jaccard_dissim),
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "Japan: Scatterplot Matrix of range geometry (H2) - Jaccard 1",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))

featurePlot(x = dta_new %>% filter(datasetID == 26) %>%
              select(all_of(c(H2))),
            y =  dta_new %>% filter(datasetID == 26) %>%
              select(Jaccard_dissim),
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "Europe: Scatterplot Matrix of range geometry (H2) - Jaccard 1",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))


#-----------------------------------------------------#
# Split range geometries (H2) for each atlas:
featurePlot(x = dta_new %>% filter(datasetID == 5) %>%
              select(all_of(c(H2))),
            y =  dta_new %>% filter(datasetID == 5) %>%
              select(log_R2_1),
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "Czechia: Scatterplot Matrix of traits (H1) - Jaccard 1",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))

featurePlot(x = dta_new %>% filter(datasetID == 6) %>%
              select(all_of(c(H2))),
            y =  dta_new %>% filter(datasetID == 6) %>%
              select(log_R2_1),
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "New York: Scatterplot Matrix of range geometry (H2) - Jaccard 1",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))

featurePlot(x = dta_new %>% filter(datasetID == 13) %>%
              select(all_of(c(H2))),
            y =  dta_new %>% filter(datasetID == 13) %>%
              select(log_R2_1),
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "Japan: Scatterplot Matrix of range geometry (H2) - Jaccard 1",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))

featurePlot(x = dta_new %>% filter(datasetID == 26) %>%
              select(all_of(c(H2))),
            y =  dta_new %>% filter(datasetID == 26) %>%
              select(log_R2_1),
            plot = "pairs",
            pch = 16,
            alpha = 0.3,
            cex = 0.5,
            xlab = "Europe: Scatterplot Matrix of range geometry (H2) - Jaccard 1",
            auto.key = list(columns = 4),
            par.settings =
              list(fontsize = list(text = 6)))
#-----------------------------------------------------#

ggplot(dta_new, aes(x = Threatened, fill = datasetID))+
  geom_bar(posititon = "identity")


ggplot(dta_new, aes(x = Habitat_5, fill = datasetID))+
  geom_bar(posititon = "identity")

