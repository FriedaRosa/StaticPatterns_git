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

------------------------------------------------------------------------

# 1. Project/File structure:

``` r
fs::dir_tree(here::here(), recurse = F)
```

    C:/Users/wolke/OneDrive - CZU v Praze/Frieda_PhD_files/02_StaticPatterns/Git
    ├── Code
    ├── Data
    ├── Demo_NewYork
    ├── Documentation
    ├── Figures
    ├── Folder_metadata.xlsx
    ├── Git.Rproj
    ├── images
    ├── Project_Description.qmd
    ├── README.html
    ├── README.md
    ├── README.qmd
    ├── README.rmarkdown
    ├── README_files
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

## A) Workflow diagram

## B) Description of steps

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

#### Modeling settings:

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
