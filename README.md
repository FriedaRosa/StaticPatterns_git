# Static Predictors - Read Me

# 0. Meta information:

- **Project title**: Static Predictors

- **Author**: Friederike Johanna Rosa Wölke, MSc

- **Date:** 2025-05-28

- **Location**: Prague, Czech Republic

- **License:** CC BY-NC-ND 4.0 (until publication)
  (https://creativecommons.org/licenses/by-nc-nd/4.0/)

- **R Package versions**: registered in file `renv.lock`

- **Computational demands:**

  - *Estimated total run time*: 50 h  
    locally on laptop without parallelization;  
    can be significantly enhanced by running predictor scripts in
    parallel or enable multiple cores to approx. 6 h

- **Manual to the folder:** `Folder_metadata.xlsx`

  - has a list of A) scripts, input and output files, figure locations,
    run times, etc. and B) Files and their sources

------------------------------------------------------------------------

### How to use this folder:

- install `renv` and `here` packages if not already installed
- open the Git.Rproj in RStudio or VScode and set `here::here()` as
  working directory to the root folder (“Git”)
- use `renv::restore()` to restore packages & dependencies from the
  lockfile (this will lead in a huge downloading session of packages)
- there are rendered .html files for all code (see/click links below to see)
- click [*here*](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/11e14ef4f0208185b1110de9ea50ad1081532c07/Project_Description.html) for advanced Project Description (incl. figs, results, predictor table)

------------------------------------------------------------------------
# Rendered code in .html
## Example: New York State
- [A_01_NY_Preprocess_data.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Demo_NewYork/Code/A_01_NY_Preprocess_data.html)
- [A_02_NY_Calculate_atlas_variables.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Demo_NewYork/Code/A_02_NY_Calculate_atlas_variables.html)
- [A_03_NY_Calculate_predictors.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Demo_NewYork/Code/A_03_NY_Calculate_predictors.html)
- [B_01_NY_Machine_learning.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Demo_NewYork/Code/B_01_NY_Machine_learning.html)

## All code (Japan, New York, Czechia, Europe)

### A: Prepare data
- [A_01_Get_data.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_01_Get_data.html)
- [A_02_Preprocess_data.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_02_Preprocess_data.html)
- [A_03_Taxon_matching_across_datasets.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_03_Taxon_matching_across_datasets.html)
- [A_04_Calculate_atlas_variables.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_04_Calculate_atlas_variables.html)
- [A_05_Predictors_avonet.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_05_Predictors_avonet.html)
- [A_06_Predictors_phylogenetic_metrics.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_06_Predictors_phylogenetic_metrics.html)
- [A_07_Predictors_IUCN.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_07_Predictors_IUCN.html)
- [A_08_Predictors_geometry.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_08_Predictors_geometry.html)
- [A_09_Predictors_spatial_autocorrelation.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_09_Predictors_spatial_autocorrelation.html)
- [A_10a_Predictors_lacunarity_rasterize_ranges.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_10a_Predictors_lacunarity_rasterize_ranges.html)
- [A_10b_Predictors_lacunarity_calculation.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_10b_Predictors_lacunarity_calculation.html)
- [A_11_Predictors_global_range_size.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_11_Predictors_global_range_size.html)
- [A_12_Merge_predictors.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_12_Merge_predictors.html)
- [A_13_Data_viz.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/A_Prepare_data/A_13_Data_viz.html)

### B: Machine learning
- [B_01_RandomForest_all_data.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/11aa958f64686347a34d8ae0110bcbd409c1eae0/Code/B_Machine_learning/B_01_RandomForest_all_data.html)
- [B_02_RandomForest_split_data.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/11aa958f64686347a34d8ae0110bcbd409c1eae0/Code/B_Machine_learning/B_02_RandomForest_split_data.html)
- [B_03_Phylogenetic_autocorrelation.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/B_Machine_learning/B_03_Phylogenetic_autocorrelation.html)
- [B_04_Model_performance_figures.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/B_Machine_learning/B_04_Model_performance_figures.html)

### C: Validation
- [C_01_Validation.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/C_Validation/C_01_Validation.html)
- [C_02_Temporal_change_simulations.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/C_Validation/C_02_Temporal_change_simulations.html)
- [C_03_Zero_trend_species.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/C_Validation/C_03_Zero_trend_species.html)

### D: Figures
- [D_01_Make_change_maps.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/D_Figures/D_01_Make_change_maps.html)
- [D_02_Figure_1.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/D_Figures/D_02_Figure_1.html)
- [D_03_Figure_2.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/D_Figures/D_03_Figure_2.html)
- [D_04_Figure_3.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/D_Figures/D_04_Figure_3.html)
- [D_05_Figure_S12a.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/D_Figures/D_05_Figure_S12a.html)
- [D_05_Figure_S12b.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/D_Figures/D_05_Figure_S12b.html)
- [D_05_Figure_S12c.html](https://rawcdn.githack.com/FriedaRosa/StaticPatterns_git/5e629d71a4a8bd07e130f08a0d9360f80f748137/Code/D_Figures/D_05_Figure_S12c.html)

------------------------------------------------------------------------
# 1. Project/File structure:

``` r
fs::dir_tree(here::here(), recurse = F)
```

    C:/Users/wolke/OneDrive - CZU v Praze/Frieda_PhD_files/02_StaticPatterns/Git
    ├── Code
    ├── Data
    ├── Demo_NewYork
    ├── Figures
    ├── Folder_metadata.xlsx
    ├── Git.Rproj
    ├── Project_Description.html
    ├── Project_Description.qmd
    ├── README.md
    ├── README.qmd
    ├── README.rmarkdown
    ├── README_files
    ├── renv
    ├── renv.lock
    └── StaticPatterns_Results_all.xlsx

There are three sections in this project: The first part (A) produces
the predictors and data needed for modelling. It starts by grabbing data
from the database, cleaning it, filtering it, ad then producing the
predictors. The second part (B) uses the predictors to train a
randomForest model and evaluate it using xAI (explainable AI),
interaction effects and variable importance. In the last part (C), the
model predictions are checked against the latest replication of the
empirical data.

Additionally, the project contains several sensitivity analyses and
robustness checks, which are not part of the main analysis but were used
to aid interpretation of the results and determine patterns of
stochasticity in the data.

### Script nomenclature:

Each R script is labelled by part (A,B,C) and script sequence (1-14).  
The 00_Configuration.R script is needed for almost all other scripts. It
ensures that packages are installed and has file paths and global
variables and lookup tables needed for many steps.

# 2. Methods summary

## A) Description of steps

1.  Get data from MOBI database for first two replications (Cz, Ny, Jp,
    Eu)

2.  Remove cells and species that were not sampled twice; filter species
    based on expert knowledge and introduced status

3.  Prepare predictors for H1 and H2, use datasetID as H3 to determine
    effect of atlas in the ‘full model’  
    <u>H1:</u> Body mass, Habitat_5, Threatened_01, Generalism_01,
    Phylodistinct, Migration_123, Global range size  
    <u>H2:</u> Fractal dimension, Lacunarity, Spatial autocorrelation,
    circularity, AOO, minimum distance to the border from the centroid  
    <u>H3:</u> datasetID

4.  Calculate responses:  

    1)  Jaccard_dissimilarity,  
    2)  log Ratio AOO,  
    3)  log Ratio AOO per year

5.  Make simulations of Jaccard_dissimilarity based on different
    combinations of parameters and evaluate the effect of these on the
    Jaccard values. Certain combinations of parameters restrict
    Jaccard_dissimilarity to a range of values. This can be used to
    determine the effect of mathematical constraints on the
    Jaccard_dissimilarity values.  
    <u>Parameters:</u>  

    1)  initially occupied cells,  
    2)  total number of cells possible,  
    3)  number of changes

6.  Train model with  

    1)  ‘all data’ and  
    2)  subsets for each datasetID (‘split data’) using random forest

7.  Extract for all three responses:  

    1)  rsq, rmse,  
    2)  hyper-parameters,  
    3)  predictions,  
    4)  variable importance,  
    5)  interactions  
    6)  partial dependence plots

8.  Test for phylogenetic autocorrelation for each datasetID in the
    model residuals (and the raw data)

9.  Calculate responses from third atlas replication (Cz, Jp), use
    predictors calculated from second period to predict responses for
    the third period and get residuals

## B) Modeling settings:

- 80/20 split (80 training, 20 testing)

- 3x repeated 10-fold cross validation

- permutation importance (not impurity)

- always split variables : datasetID

- respect unordered factors = T

- Bayesian hyperparameter tuning:

  - mtry = 2-10

  - min_n = 5-15

  - trees = 1000-5000

  - initial values = 5

  - iterations = 50

  - no improve = 10

  - set a seed.

## C) Data overview

### Responses:

- Jaccard_dissimilarity

- log_R2_1 (log ratio between sampling period 2 and 1)

- log_R2_1_per_year (log ratio between sampling period 2 and 1 divided
  by the number of years between sampling)

#### Notes:

- The higher J_dissim, the more variable log_R2_1 and log_R2_1_per_year

- The smaller AOO, the more variable log_R2_1 and log_R2_1_per_year

- The lower D, the more variable (and more positive?) log_R2_1

- The higher mean_lnLac, the more variable (and more positive?) log_R2_1

- The higher mean_lnLac, the lower Jaccard_dissim

- Species in New York and Japan are more dissimilar than species in
  Czechia and Europe


  
