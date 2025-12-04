# Bayesian Analysis for Predicting Radiation-Induced Liver Disease (RILD)

This project implements a Bayesian statistical analysis to predict Radiation-Induced Liver Disease (RILD) using clinical, laboratory, and dosimetric data.

## Overview

This project implements **advanced Bayesian methods** for predicting RILD, including:

### Core Methods (Advanced):
- **Fully Bayesian missing data handling** (within-model imputation using Stan)
- **Multilevel/hierarchical models** with random intercepts and slopes
- **Bayesian model averaging** and ensemble predictions
- **Bayesian variable selection** using horseshoe priors
- **Posterior predictive simulations** for clinical scenarios

### Standard Bayesian Methods:
- **Bayesian logistic regression** using Stan (via brms)
- **Model diagnostics** including MCMC convergence checks
- **Posterior inference** with credible intervals
- **Model performance evaluation** (ROC curves, confusion matrices)
- **Sensitivity analysis** with alternative priors
- **Leave-one-out cross-validation** for model comparison

## Data Files

- `data/baseline_labs_anonymized.csv`: Laboratory values (platelet count, ALP, bilirubin, albumin, AST, ALT, INR)
- `data/Clinical_Features_anonymized.csv`: Patient demographics and clinical features
- `data/dose_volume_histogram_5gy_bins_anonymized.csv`: Dose-volume histogram data
- `data/patients_to_process.csv`: RILD outcome variable (0 = no RILD, 1 = RILD)

## Required R Packages

**First, install all required packages:**

### Option 1: Run the setup script (Recommended)
In R or RStudio, run:
```r
source("setup_packages.R")
```

### Option 2: Install manually
```r
# Install from CRAN
install.packages(c("tidyverse", "brms", "bayesplot", "mice", "VIM", 
                   "corrplot", "pROC", "loo"), 
                 repos = "https://cran.rstudio.com/")
```

### Option 3: Use the install script
```bash
Rscript install_packages.R
```

**Important Notes:**
- `brms` requires RStan, which may take 10-15 minutes to install
- If you get "package not found" errors, try running `setup_packages.R` in your R session
- Some packages may require compilation (Xcode Command Line Tools on Mac)
- See [brms documentation](https://github.com/paul-buerkner/brms) for troubleshooting

## Running the Analysis

1. Ensure all data files are in the `data/` folder
2. Open `bayesian_rild_analysis.R` in R or RStudio
3. Run the script (this may take 10-30 minutes depending on your system)

The script will:
- Load and clean the data
- Perform Bayesian multiple imputation
- Fit Bayesian logistic regression models
- Generate diagnostic plots and performance metrics
- Save model objects for later use

## Output Files

The analysis generates several output files:

- **Diagnostic plots:**
  - `missing_data_pattern.png`: Visualization of missing data
  - `mcmc_trace_plots.png`: MCMC convergence diagnostics
  - `posterior_density_plots.png`: Posterior distributions
  - `energy_plot.png`: NUTS energy diagnostics
  - `roc_curve.png`: ROC curve for model performance
  - `coefficient_posteriors.png`: Effect sizes with credible intervals
  - `posterior_predictive_check.png`: Model fit assessment

- **Model objects:**
  - `bayesian_rild_model.rds`: Main Bayesian model
  - `bayesian_rild_model_informative.rds`: Model with informative priors

## Key Features

### Fully Bayesian Missing Data Handling
- **Within-model imputation**: Missing data handled directly in Stan (more sophisticated than pre-imputation)
- Uses multivariate normal model for missing covariates
- Uncertainty in missing values properly propagated through the model

### Multilevel/Hierarchical Models
- **Random intercepts by Diagnosis**: Allows diagnosis-specific baseline RILD risk
- **Random slopes**: Allows key predictors (e.g., EUD) to have diagnosis-specific effects
- Properly accounts for clustering and heterogeneity across diagnosis groups
- Uses LKJ priors for correlation structure in random effects

### Bayesian Model Averaging
- Combines multiple models using LOO-CV weights
- Ensemble predictions improve over individual models
- Properly accounts for model uncertainty

### Bayesian Variable Selection
- Uses horseshoe priors for automatic variable selection
- Shrinks unimportant variables toward zero
- Provides sparsity-inducing regularization

### Standard Bayesian Methods
- Uses weakly informative priors (normal(0, 1) for coefficients)
- 4 MCMC chains with 2000 iterations (1000 warmup)
- Includes convergence diagnostics (R-hat, trace plots, energy plots)

### Statistical Analysis
- Posterior summaries with credible intervals
- Effect sizes (odds ratios) with uncertainty quantification
- Model performance metrics (AUC, sensitivity, specificity)
- Leave-one-out cross-validation for model comparison
- Sensitivity analysis with alternative priors

## Models Implemented

1. **Fixed Effects Logistic Regression**: Standard Bayesian logistic regression
2. **Fully Bayesian Missing Data Model**: Handles missing data within Stan
3. **Multilevel Model (Random Intercepts)**: Hierarchical structure by Diagnosis
4. **Multilevel Model (Random Slopes)**: Random slopes for EUD by Diagnosis
5. **Variable Selection Model**: Horseshoe priors for automatic selection
6. **Ensemble Model**: Bayesian model averaging of all models

## Model Variables

**Predictors included:**
- Age
- Gender
- Diagnosis type (used for random effects)
- Laboratory values: platelet count, bilirubin, albumin, INR
- Dosimetric parameters: PTV volume, EUD, mean dose, V30
- Derived features: liver volume ratio, AST/ALT ratio

**Outcome:**
- RILD (binary: 0 = no RILD, 1 = RILD)

## Advanced Statistical Methods

This project demonstrates **advanced Bayesian methods** suitable for graduate-level coursework:

1. **Missing Data (Level 8)**: Fully Bayesian imputation within the model
2. **Multilevel Models (Level 6)**: Hierarchical structure with random effects
3. **Model Averaging**: Bayesian ensemble methods
4. **Variable Selection**: Sparsity-inducing priors
5. **Posterior Predictive Simulations**: Clinical scenario analysis

## Notes

- The analysis handles special characters in lab data (e.g., "<12.", "See Note")
- Missing data handled both via pre-imputation (mice) and fully Bayesian (within Stan)
- Model convergence checked using R-hat statistics (target: < 1.01)
- All plots saved as PNG files for easy inclusion in reports
- Multiple models compared using LOO-CV and information criteria
- Random effects properly account for diagnosis-specific heterogeneity

## Troubleshooting

If you encounter issues:

1. **RStan installation problems:** Follow the [RStan installation guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)

2. **Memory issues:** Reduce the number of iterations or chains in the `brm()` calls

3. **Convergence warnings:** Increase `adapt_delta` or `max_treedepth` in the `control` argument

4. **Missing packages:** Install any missing packages using `install.packages()`

## Citation

This analysis was developed for an Applied Bayesian Analysis course project.

