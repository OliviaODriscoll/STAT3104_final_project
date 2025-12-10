# Test Script: Baseline Model Without Radiomics
# This script tests just the baseline model with diagnostics

# ============================================================================
# 1. SETUP AND DATA LOADING
# ============================================================================

library(tidyverse)
library(brms)          # Bayesian regression models with Stan (supports Laplace priors)
library(bayesplot)     # MCMC diagnostics and posterior plots
library(mice)          # Multiple imputation
library(pROC)          # ROC curves
library(loo)           # Leave-One-Out cross-validation
library(dagitty)       # Directed Acyclic Graphs for causal inference
library(posterior)     # For extracting posterior samples

set.seed(42)

cat("=== BASELINE MODEL (NO RADIOMICS) ===\n\n")
baseline_labs <- read.csv("data/baseline_labs_anonymized.csv", stringsAsFactors = FALSE)
clinical_features <- read.csv("data/Clinical_Features_anonymized.csv", stringsAsFactors = FALSE)
dvh_data <- read.csv("data/dose_volume_histogram_5gy_bins_anonymized.csv", stringsAsFactors = FALSE)
rild_outcome <- read.csv("data/patients_to_process.csv", stringsAsFactors = FALSE)


clean_numeric <- function(x) {
  x <- gsub("<", "", x)
  x <- gsub("See Note", NA, x)
  x <- gsub(",", "", x)
  as.numeric(x)
}

baseline_labs_clean <- baseline_labs %>%
  mutate_at(vars(-MRN), clean_numeric) %>%
  rename(
    platelet = `PLATELET.COUNT`,
    alp = ALP,
    bilirubin = `BILIRUBIN.TOTAL`,
    albumin = `ALBUMIN.LEVEL`,
    ast = AST,
    alt = ALT,
    inr = INR
  )

clinical_features_clean <- clinical_features %>%
  mutate(
    Age = as.numeric(Age),
    Gender = factor(Gender),
    Diagnosis = factor(Diagnosis),
    PTV_volume_cc = as.numeric(PTV_volume_cc),
    Liver_minus_PTV_volume_cc = as.numeric(Liver_minus_PTV_volume_cc),
    EUD_alpha_0.06 = as.numeric(EUD_alpha_0.06)
  )

dvh_clean <- dvh_data %>%
  mutate_at(vars(-MRN), as.numeric)

rild_clean <- rild_outcome %>%
  mutate(RILD = as.numeric(RILD))

merged_base <- baseline_labs_clean %>%
  full_join(clinical_features_clean, by = "MRN") %>%
  full_join(dvh_clean, by = "MRN") %>%
  full_join(rild_clean, by = "MRN") %>%
  filter(!is.na(RILD))

# Calculate derived features
dvh_vars <- grep("^[0-9]", names(merged_base), value = TRUE)
if (length(dvh_vars) == 0) {
  dvh_vars <- grep("Gy$", names(merged_base), value = TRUE)
}

if (length(dvh_vars) > 0 && nrow(merged_base) > 0) {
  dvh_mat <- as.matrix(merged_base[, dvh_vars, drop = FALSE])
  bin_centers <- seq(2.5, 2.5 + (ncol(dvh_mat) - 1) * 5, by = 5)
  weighted_sum <- rowSums(sweep(dvh_mat, 2, bin_centers, "*"), na.rm = TRUE)
  total_volume <- rowSums(dvh_mat, na.rm = TRUE)
  mean_dose_calc <- ifelse(total_volume > 0, weighted_sum / total_volume, NA)
} else {
  mean_dose_calc <- rep(NA, nrow(merged_base))
}

merged_base <- merged_base %>%
  mutate(
    mean_dose = mean_dose_calc,
    liver_volume_ratio = Liver_minus_PTV_volume_cc / 
                        (PTV_volume_cc + Liver_minus_PTV_volume_cc),
    ast_alt_ratio = ast / alt
  ) %>%
  mutate_all(~ifelse(is.infinite(.), NA, .))

cat("Data:", nrow(merged_base), "patients,", 
    sum(merged_base$RILD == 1, na.rm = TRUE), "RILD cases\n")


# One-hot encode function
# Binary variables (2 levels) are converted to numeric 0/1
# Multi-category variables (>2 levels) are one-hot encoded
one_hot_encode <- function(data, categorical_vars = c("Gender", "Diagnosis")) {
  data_encoded <- data
  for (var in categorical_vars) {
    if (var %in% names(data_encoded)) {
      if (!is.factor(data_encoded[[var]])) {
        data_encoded[[var]] <- factor(data_encoded[[var]])
      }
      levels <- levels(data_encoded[[var]])
      n_levels <- length(levels)
      if (n_levels > 1) {
        if (n_levels == 2) {
          # Binary variable: convert to numeric 0/1 (keep original variable name)
          data_encoded[[var]] <- as.numeric(data_encoded[[var]]) - 1
        } else {
          # Multi-category variable: one-hot encode (create new columns)
          for (i in 2:n_levels) {
            level_name <- levels[i]
            var_name <- paste0(var, "_", gsub("[^A-Za-z0-9_]", "_", level_name))
            data_encoded[[var_name]] <- as.numeric(data_encoded[[var]] == level_name)
          }
          data_encoded[[var]] <- NULL
        }
      }
    }
  }
  return(data_encoded)
}

key_vars <- c("Age", "Gender", "Diagnosis", "platelet", "alp", "bilirubin", 
              "albumin", "ast", "alt", "inr", "PTV_volume_cc", 
              "Liver_minus_PTV_volume_cc", "EUD_alpha_0.06", "mean_dose", 
              "liver_volume_ratio", "ast_alt_ratio")

key_vars <- key_vars[key_vars %in% names(merged_base)]

# Prepare data - KEEP NAs for Bayesian imputation
model_data <- merged_base %>%
  select(MRN, RILD, all_of(key_vars)) %>%
  filter(!is.na(RILD))

# Report missing data pattern
cat("Missing data pattern (before Bayesian imputation):\n")
missing_counts <- colSums(is.na(model_data[, key_vars, drop = FALSE]))
if (sum(missing_counts) > 0) {
  cat("Missing values:", sum(missing_counts), "across", sum(missing_counts > 0), "variables\n")
}
numeric_vars <- model_data %>%
  select(-MRN, -RILD) %>%  # Explicitly exclude outcome
  select_if(is.numeric) %>%
  names()


# Store standardization parameters for use in model
standardization_params <- list()

for (var in numeric_vars) {
  if (is.numeric(model_data[[var]])) {
    # Compute mean/sd from PREDICTORS ONLY (no outcome variable)
    var_mean <- mean(model_data[[var]], na.rm = TRUE)
    var_sd <- sd(model_data[[var]], na.rm = TRUE)
    standardization_params[[var]] <- list(mean = var_mean, sd = var_sd)
    
    if (!is.na(var_sd) && var_sd > 1e-10) {
      # Standardize observed values, keep NAs as NA
      model_data[[var]] <- (model_data[[var]] - var_mean) / var_sd
    } else {
      # Zero variance - set to 0, but keep NAs
      model_data[[var]][!is.na(model_data[[var]])] <- 0
    }
  }
}


# One-hot encode
model_data <- one_hot_encode(model_data, categorical_vars = c("Gender", "Diagnosis"))

# Remove MRN and prepare final data
model_data_final <- model_data %>%
  select(-MRN) %>%
  mutate(RILD = factor(RILD, levels = c(0, 1)))

# Verify no radiomics variables are present
all_vars <- names(model_data_final)
radiomics_vars <- grep("radiomic|texture|shape|firstorder|glcm|glrlm|glszm|original_|wavelet_|log_|CTORMid|Primary", 
                       all_vars, ignore.case = TRUE, value = TRUE)
if (length(radiomics_vars) > 0) {
  model_data_final <- model_data_final %>% select(-all_of(radiomics_vars))
}

# Prepare data for ulam with Bayesian imputation
# Keep NAs - they will be estimated as parameters
n_obs <- nrow(model_data_final)
data_list <- list()
data_list$RILD <- as.numeric(model_data_final$RILD) - 1
data_list$N <- n_obs

# Get predictor names (exclude RILD)
predictor_names <- setdiff(names(model_data_final), "RILD")

# Track which variables have missing values
vars_with_missing <- c()
missing_indicators <- list()

for (col in predictor_names) {
  if (is.factor(model_data_final[[col]])) {
    data_list[[col]] <- as.numeric(model_data_final[[col]]) - 1
  } else {
    data_list[[col]] <- as.numeric(model_data_final[[col]])
  }
  
    if (any(is.na(data_list[[col]]))) {
      vars_with_missing <- c(vars_with_missing, col)
      missing_indicators[[col]] <- as.numeric(is.na(data_list[[col]]))
    }
}

cat("\nVariables with missing values:", length(vars_with_missing), "\n")
if (length(vars_with_missing) > 0) {
  cat("Missing values will be estimated as parameters in the Stan model.\n\n")
}

# Sanitize variable names
name_mapping <- data.frame(
  original = predictor_names,
  sanitized = paste0("X", 1:length(predictor_names)),
  stringsAsFactors = FALSE
)

data_list_sanitized <- list()
data_list_sanitized$RILD <- data_list$RILD

for (i in 1:nrow(name_mapping)) {
  orig_name <- name_mapping$original[i]
  sanitized_name <- name_mapping$sanitized[i]
  if (orig_name %in% names(data_list)) {
    data_list_sanitized[[sanitized_name]] <- data_list[[orig_name]]
    cat("  ", orig_name, "->", sanitized_name, "\n")
  }
}

predictor_names_sanitized <- name_mapping$sanitized

# Calculate prevalence for intercept
rild_prevalence <- mean(data_list_sanitized$RILD)
intercept_mean <- if (rild_prevalence > 0 && rild_prevalence < 1) {
  qlogis(rild_prevalence)
} else {
  0
}

cat("\nRILD prevalence:", round(rild_prevalence, 3), 
    "-> intercept prior mean:", round(intercept_mean, 3), "\n\n")


cat("\n=== BAYESIAN IMPUTATION SETUP ===\n\n")

cat("Imputing missing values with MICE...\n")


# Prepare data for MICE imputation
# Convert data_list to data frame for MICE
# Only include predictor variables that exist in data_list and have the correct length
predictor_names_in_list <- intersect(predictor_names, names(data_list))
# Verify all have the same length
valid_predictors <- sapply(predictor_names_in_list, function(x) {
  length(data_list[[x]]) == n_obs
})
predictor_names_valid <- predictor_names_in_list[valid_predictors]

if (length(predictor_names_valid) == 0) {
  stop("No valid predictor variables found in data_list for MICE imputation")
}

impute_df <- as.data.frame(data_list[predictor_names_valid])
impute_df$RILD <- data_list$RILD  # Include outcome for better imputation

# Check missing data pattern
missing_counts <- colSums(is.na(impute_df))
if (sum(missing_counts) > 0) {
  cat("Missing data pattern before MICE:\n")
  print(missing_counts[missing_counts > 0])
  cat("\nTotal missing values:", sum(missing_counts), "\n\n")
  
  cat("Running MICE imputation (this may take a minute)...\n")
  
  mice_result <- tryCatch({
    mice(
      impute_df,
      m = 1,                    # Single imputation (we'll use first imputation)
      method = "pmm",           # Predictive Mean Matching for numeric
      maxit = 20,               # 20 iterations (Gibbs sampling iterations)
      seed = 42,
      printFlag = FALSE
    )
  }, error = function(e) {
    cat("MICE failed with error:", e$message, "\n")
    cat("Falling back to mean imputation...\n")
    NULL
  })
  
  if (!is.null(mice_result)) {
    # Extract first imputation
    impute_df_complete <- complete(mice_result, 1)
    
    # Update data_list with imputed values
    for (col in predictor_names_valid) {
      if (col %in% names(impute_df_complete)) {
        if (col %in% names(missing_counts) && missing_counts[col] > 0) {
          n_imputed <- missing_counts[col]
          data_list[[col]] <- impute_df_complete[[col]]
          cat("  ", col, ": imputed", n_imputed, "values\n")
        } else {
          data_list[[col]] <- impute_df_complete[[col]]
        }
      }
    }
    
    cat("\n✓ MICE imputation completed successfully\n")
    
    # Verify no NAs remain
    remaining_nas <- sum(is.na(impute_df_complete[predictor_names_valid]))
    if (remaining_nas > 0) {
      cat("WARNING:", remaining_nas, "NAs still remain after MICE\n")
      cat("These will be imputed with mean values\n")
      
      # Fallback: mean imputation for any remaining NAs
      for (col in predictor_names) {
        if (any(is.na(data_list[[col]]))) {
          obs_vals <- data_list[[col]][!is.na(data_list[[col]])]
          if (length(obs_vals) > 0) {
            data_list[[col]][is.na(data_list[[col]])] <- mean(obs_vals, na.rm = TRUE)
          } else {
            data_list[[col]][is.na(data_list[[col]])] <- 0
          }
        }
      }
    }
  } else {
    # Fallback to mean imputation if MICE fails
    cat("Using mean imputation as fallback...\n")
    for (col in predictor_names) {
      if (col %in% names(data_list)) {
        if (any(is.na(data_list[[col]]))) {
          obs_vals <- data_list[[col]][!is.na(data_list[[col]])]
          if (length(obs_vals) > 0) {
            data_list[[col]][is.na(data_list[[col]])] <- mean(obs_vals, na.rm = TRUE)
          } else {
            data_list[[col]][is.na(data_list[[col]])] <- 0
          }
        }
      }
    }
  }
} else {
  cat("✓ No missing values found - skipping imputation\n")
}

# Update sanitized data list with imputed values
for (i in 1:nrow(name_mapping)) {
  orig_name <- name_mapping$original[i]
  sanitized_name <- name_mapping$sanitized[i]
  if (orig_name %in% names(data_list)) {
    data_list_sanitized[[sanitized_name]] <- data_list[[orig_name]]
  }
}


lambda_lasso <- 0.5
coef_terms <- paste0("b_", predictor_names_sanitized, " * ", predictor_names_sanitized)
coef_terms_str <- paste(coef_terms, collapse = " + ")

laplace_scale <- 1 / lambda_lasso
model_data_brms <- as.data.frame(data_list_sanitized)
model_data_brms$RILD <- as.factor(model_data_brms$RILD)

# Build formula string
formula_str <- paste("RILD ~", paste(predictor_names_sanitized, collapse = " + "))

laplace_scale_num <- as.numeric(round(laplace_scale, 3))
intercept_mean_num <- as.numeric(round(intercept_mean, 3))


laplace_prior <- eval(parse(text = paste0("prior(double_exponential(0, ", laplace_scale_num, "), class = 'b')")))
intercept_prior <- eval(parse(text = paste0("prior(normal(", intercept_mean_num, ", 2.5), class = 'Intercept')")))

# Combine priors
model_priors <- c(laplace_prior, intercept_prior)

X_df <- as.data.frame(data_list_sanitized[predictor_names_sanitized])
X_matrix <- as.matrix(sapply(X_df, as.numeric))
y <- as.numeric(data_list_sanitized$RILD)

cor_matrix <- cor(X_matrix, use = "pairwise.complete.obs")
high_cor <- which(abs(cor_matrix) > 0.9 & abs(cor_matrix) < 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  for (i in 1:min(10, nrow(high_cor))) {
    row_idx <- high_cor[i, 1]
    col_idx <- high_cor[i, 2]
    if (row_idx < col_idx) {  # Avoid duplicates
      var1 <- name_mapping$original[row_idx]
      var2 <- name_mapping$original[col_idx]
      cor_val <- cor_matrix[row_idx, col_idx]
      cat(sprintf("    %s <-> %s: %.3f\n", var1, var2, cor_val))
    }
  }
} 


  
model <- brm(
  formula = as.formula(formula_str),
  data = model_data_brms,
  family = bernoulli(link = "logit"),
  prior = model_priors,
  chains = 4,
  iter = 5000,
  warmup = 2000,
  cores = 4,
  seed = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  file = "model_full_all_predictors_brms"  # Save model
)


# Save the full model
saveRDS(model, "model_full_all_predictors.rds")
cat("Full model saved to: model_full_all_predictors.rds\n")
cat("  - Contains all", length(predictor_names_sanitized), "predictors\n")
cat("  - Model name: 'model' (full model with all predictors)\n\n")



# Extract samples for correlation analysis (brms)
post <- as_draws_df(model)
post_df <- as.data.frame(post)

# brms uses b_Intercept for intercept, ensure we have 'a' for compatibility
if ("b_Intercept" %in% names(post_df) && !"a" %in% names(post_df)) {
  post_df$a <- post_df$b_Intercept
}

# Check parameter correlations in posterior
coef_params <- grep("^b_", names(post_df), value = TRUE)
if (length(coef_params) > 1) {
  coef_matrix <- as.matrix(post_df[coef_params])
  post_cor <- cor(coef_matrix)
  
  # Find high correlations
  high_cor_post <- which(abs(post_cor) > 0.7 & abs(post_cor) < 1, arr.ind = TRUE)
  if (nrow(high_cor_post) > 0) {
    cat("  WARNING: High posterior correlations detected (>0.7):\n")
    shown <- 0
    for (i in 1:nrow(high_cor_post)) {
      row_idx <- high_cor_post[i, 1]
      col_idx <- high_cor_post[i, 2]
      if (row_idx < col_idx && shown < 10) {  # Avoid duplicates, limit output
        param1 <- coef_params[row_idx]
        param2 <- coef_params[col_idx]
        cor_val <- post_cor[row_idx, col_idx]
        orig1 <- name_mapping$original[name_mapping$sanitized == gsub("^b_", "", param1)]
        orig2 <- name_mapping$original[name_mapping$sanitized == gsub("^b_", "", param2)]
        cat(sprintf("    %s <-> %s: %.3f\n", 
                    ifelse(length(orig1) > 0 && !is.na(orig1), orig1, param1),
                    ifelse(length(orig2) > 0 && !is.na(orig2), orig2, param2),
                    cor_val))
        shown <- shown + 1
      }
    }
} 

# Extract samples for diagnostics (brms uses as_draws_df from posterior package)
post <- as_draws_df(model)

# Convert to data frame format compatible with existing code
post_df <- as.data.frame(post)
# Rename columns to match expected format (brms uses b_Intercept, b_X1, etc.)
# Remove chain/iteration columns for coefficient extraction
coef_cols <- grep("^b_", names(post_df), value = TRUE)
post_df <- post_df[, c(coef_cols, "lp__"), drop = FALSE]
if ("b_Intercept" %in% names(post_df) && !"a" %in% names(post_df)) {
  post_df$a <- post_df$b_Intercept
}

# Get parameter names
param_names <- names(post_df)
all_params <- c("a", grep("^b_", param_names, value = TRUE))

# Get total number of samples (4 chains * 3000 post-warmup iterations)
n_chains <- 4
n_warmup <- 2000
n_post_warmup <- 3000
total_samples <- n_chains * n_post_warmup

model_summary <- summary(model)
print(model_summary)
precis_output <- model_summary$fixed  # Extract fixed effects for compatibility


# Convert to format compatible with existing code
if (!is.null(precis_output) && is.data.frame(precis_output)) {
  # Create a data frame with R-hat and ESS in expected format
  # brms uses "Rhat" and "Bulk_ESS" column names
  precis_df <- data.frame(
    row.names = rownames(precis_output),
    mean = precis_output$Estimate,
    sd = precis_output$Est.Error,
    Rhat = if("Rhat" %in% colnames(precis_output)) precis_output$Rhat else NA,
    Rhat4 = if("Rhat" %in% colnames(precis_output)) precis_output$Rhat else NA,
    ess_bulk = if("Bulk_ESS" %in% colnames(precis_output)) precis_output$Bulk_ESS else NA,
    n_eff = if("Bulk_ESS" %in% colnames(precis_output)) precis_output$Bulk_ESS else NA,
    stringsAsFactors = FALSE
  )
    
  # Print R-hat and ESS summary
  if ("Rhat" %in% colnames(precis_output)) {
    rhat_vals <- precis_output$Rhat
    cat("\nR-hat summary:\n")
    cat("  Max R-hat:", round(max(rhat_vals, na.rm = TRUE), 4), "\n")
    cat("  Min R-hat:", round(min(rhat_vals, na.rm = TRUE), 4), "\n")
    high_rhat <- rhat_vals[rhat_vals > 1.01 & !is.na(rhat_vals)]
    if (length(high_rhat) > 0) {
      cat("  WARNING:", length(high_rhat), "parameters with R-hat > 1.01\n")
    } else {
      cat("  ✓ All R-hat values < 1.01 (good convergence)\n")
    }
  }
  
  if ("Bulk_ESS" %in% colnames(precis_output)) {
    ess_vals <- precis_output$Bulk_ESS
    cat("\nEffective Sample Size (ESS) summary:\n")
    cat("  Min ESS:", round(min(ess_vals, na.rm = TRUE), 0), "\n")
    cat("  Max ESS:", round(max(ess_vals, na.rm = TRUE), 0), "\n")
    low_ess <- ess_vals[ess_vals < 1200 & !is.na(ess_vals)]
    if (length(low_ess) > 0) {
      cat("  WARNING:", length(low_ess), "parameters with low ESS (< 1200)\n")
    } else {
      cat("  ✓ All ESS values >= 1200 (good sampling)\n")
    }
  }
} else {
  precis_df <- NULL
  cat("\nCould not extract R-hat and ESS from model summary.\n")
}

# Effective sample size
if (!is.null(precis_output) && "n_eff" %in% names(precis_output)) {
  neff_vals <- precis_output$n_eff
  cat("\nEffective sample size:\n")
  cat("  Min Neff:", round(min(neff_vals, na.rm = TRUE), 0), "\n")
  cat("  Max Neff:", round(max(neff_vals, na.rm = TRUE), 0), "\n")
  low_neff <- neff_vals[neff_vals < 1200]  # 1200 = 0.1 * 12000 (4 chains * 3000 samples)
  if (length(low_neff) > 0) {
    cat("  WARNING:", length(low_neff), "parameters with low effective sample size (< 400)\n")
  } else {
    cat("  ✓ All Neff values > 400 (good sampling)\n")
  }
}

cat("\n=== COEFFICIENT SUMMARIES ===\n\n")

# Create comprehensive coefficient table
# brms uses b_Intercept for intercept, and b_X1, b_X2, etc. for coefficients
coef_params <- grep("^b_", names(post_df), value = TRUE)
coef_params <- coef_params[coef_params != "b_Intercept"]  # Remove intercept from coef list
all_params <- c("a", coef_params)  # "a" will map to b_Intercept

# Extract precis output for R-hat and ESS
precis_df <- as.data.frame(precis_output)

# Create coefficient summary table
coef_table <- data.frame(
  Parameter = all_params,
  Original_Name = NA,
  Mean = NA,
  SD = NA,
  CI_Lower_2.5 = NA,
  CI_Upper_97.5 = NA,
  CI_Lower_5.5 = NA,
  CI_Upper_94.5 = NA,
  Rhat = NA,
  ESS = NA,
  stringsAsFactors = FALSE
)

# Fill in the table
for (i in 1:nrow(coef_table)) {
  param <- coef_table$Parameter[i]
  
  # Map parameter names: "a" -> "b_Intercept" for brms
  param_name <- if(param == "a") "b_Intercept" else param
  
  if (param_name %in% names(post_df)) {
    param_samples <- post_df[[param_name]]
    
    # Basic statistics
    coef_table$Mean[i] <- mean(param_samples)
    coef_table$SD[i] <- sd(param_samples)
    coef_table$CI_Lower_2.5[i] <- quantile(param_samples, 0.025)
    coef_table$CI_Upper_97.5[i] <- quantile(param_samples, 0.975)
    coef_table$CI_Lower_5.5[i] <- quantile(param_samples, 0.055)
    coef_table$CI_Upper_94.5[i] <- quantile(param_samples, 0.945)
    
    # Original name
    if (param == "a") {
      coef_table$Original_Name[i] <- "Intercept"
    } else {
      sanitized <- gsub("^b_", "", param)
      orig <- name_mapping$original[name_mapping$sanitized == sanitized]
      if (length(orig) > 0 && !is.na(orig) && orig != "") {
        coef_table$Original_Name[i] <- orig
      } else {
        coef_table$Original_Name[i] <- sanitized
      }
    }
    
    if (param == "a") {
      precis_row_name <- "Intercept"
    } else {
      # Remove "b_" prefix: "b_X6" -> "X6"
      precis_row_name <- gsub("^b_", "", param)
    }
    
    # Try to extract from precis_output first (original brms summary)
    if (!is.null(precis_output) && is.data.frame(precis_output) && 
        precis_row_name %in% rownames(precis_output)) {
      if ("Rhat" %in% colnames(precis_output)) {
        coef_table$Rhat[i] <- as.numeric(precis_output[precis_row_name, "Rhat"])
      }
      if ("Bulk_ESS" %in% colnames(precis_output)) {
        coef_table$ESS[i] <- as.numeric(precis_output[precis_row_name, "Bulk_ESS"])
      }
    } else if (!is.null(precis_df) && precis_row_name %in% rownames(precis_df)) {
      # Fallback to precis_df if precis_output doesn't have it
      if ("Rhat" %in% colnames(precis_df)) {
        coef_table$Rhat[i] <- as.numeric(precis_df[precis_row_name, "Rhat"])
      } else if ("Rhat4" %in% colnames(precis_df)) {
        coef_table$Rhat[i] <- as.numeric(precis_df[precis_row_name, "Rhat4"])
      }
      if ("ess_bulk" %in% colnames(precis_df)) {
        coef_table$ESS[i] <- as.numeric(precis_df[precis_row_name, "ess_bulk"])
      } else if ("n_eff" %in% colnames(precis_df)) {
        coef_table$ESS[i] <- as.numeric(precis_df[precis_row_name, "n_eff"])
      }
    }
  }
}

# Sort by absolute mean (excluding intercept)
coef_table_sorted <- coef_table %>%
  mutate(Abs_Mean = abs(Mean)) %>%
  arrange(desc(Abs_Mean))


for (i in 1:nrow(coef_table_sorted)) {
  row <- coef_table_sorted[i, ]
  cat(sprintf("%-15s %-30s %10.4f %10.4f %12.4f %12.4f %8.4f %10.0f\n",
              row$Parameter,
              substr(row$Original_Name, 1, 29),
              row$Mean,
              row$SD,
              row$CI_Lower_2.5,
              row$CI_Upper_97.5,
              ifelse(is.na(row$Rhat), NA, row$Rhat),
              ifelse(is.na(row$ESS), NA, row$ESS)))
}


}

# Save coefficient table to CSV
write.csv(coef_table_sorted, "full_model_coefficients.csv", row.names = FALSE)

# Save full model
saveRDS(model, "model_full_all_predictors.rds")

# Intercept summary
intercept_row <- coef_table_sorted[coef_table_sorted$Parameter == "a", ]

# Plot 1: Coefficient means with credible intervals
png("coefficient_plot_with_intervals.png", width = 2000, height = 1200, res = 150)

# Sort coefficients by mean value (excluding intercept)
coef_plot_data <- coef_table_sorted %>%
  filter(Parameter != "a") %>%
  arrange(Mean)

n_coefs <- nrow(coef_plot_data)
y_pos <- 1:n_coefs

par(mar = c(5, 12, 4, 2))
plot(coef_plot_data$Mean, y_pos, 
     xlim = c(min(coef_plot_data$CI_Lower_2.5) * 1.1, 
              max(coef_plot_data$CI_Upper_97.5) * 1.1),
     xlab = "Coefficient Value", ylab = "",
     yaxt = "n", type = "n",
     main = "Coefficient Estimates with 95% Credible Intervals")

# Add credible intervals
segments(coef_plot_data$CI_Lower_2.5, y_pos,
         coef_plot_data$CI_Upper_97.5, y_pos,
         col = "gray70", lwd = 2)

# Add points for means
points(coef_plot_data$Mean, y_pos, pch = 19, col = "steelblue", cex = 1.2)

# Add vertical line at 0
abline(v = 0, col = "red", lty = 2, lwd = 2)

# Add labels
axis(2, at = y_pos, labels = coef_plot_data$Original_Name, las = 2, cex.axis = 0.7)

# Color code by significance
sig_indices <- which((coef_plot_data$CI_Lower_2.5 > 0 & coef_plot_data$CI_Upper_97.5 > 0) |
                     (coef_plot_data$CI_Lower_2.5 < 0 & coef_plot_data$CI_Upper_97.5 < 0))
if (length(sig_indices) > 0) {
  points(coef_plot_data$Mean[sig_indices], y_pos[sig_indices], 
         pch = 19, col = "darkgreen", cex = 1.5)
}

legend("bottomright", 
       legend = c("Mean", "95% CI", "Significant (CI excludes 0)", "Zero line"),
       col = c("steelblue", "gray70", "darkgreen", "red"),
       pch = c(19, NA, 19, NA),
       lty = c(NA, 1, NA, 2),
       lwd = c(NA, 2, NA, 2))

dev.off()
cat("Coefficient plot saved to: coefficient_plot_with_intervals.png\n\n")

# Create individual trace plots for all parameters
cat("Generating individual trace plots for all parameters...\n")
for (param in all_params) {
  if (param %in% names(post_df)) {
    param_samples <- post_df[[param]]
    orig_name <- coef_table$Original_Name[coef_table$Parameter == param]
    
    # Create filename
    filename <- paste0("traceplot_", param, ".png")
    png(filename, width = 1200, height = 800, res = 150)
    par(mar = c(5, 5, 4, 2))
    
    # Create trace plot
    plot(1:length(param_samples), param_samples, type = "l", 
         main = paste("Trace Plot:", param, "\n", orig_name),
         xlab = "Sample", ylab = "Value", 
         col = "steelblue", lwd = 1)
    
    # Add mean line
    abline(h = mean(param_samples), col = "red", lty = 2, lwd = 2)
    
    # Add 95% credible interval
    ci_lower <- quantile(param_samples, 0.025)
    ci_upper <- quantile(param_samples, 0.975)
    abline(h = c(ci_lower, ci_upper), col = "orange", lty = 3, lwd = 1.5)
    
    # Add summary stats
    legend_text <- sprintf("Mean=%.3f\nSD=%.3f\n95%% CI: [%.3f, %.3f]", 
                          mean(param_samples), sd(param_samples),
                          ci_lower, ci_upper)
    legend("topright", legend = legend_text, bty = "n", cex = 1)
    
    dev.off()
  }
}
cat("All parameter trace plots saved as individual files\n\n")

# Save coefficient table to CSV
# Save full model coefficient table
write.csv(coef_table_sorted, "full_model_coefficients.csv", row.names = FALSE)

# Get predictions (brms uses posterior_linpred)
pred_linpred <- posterior_linpred(model, transform = FALSE)
pred_probs <- plogis(pred_linpred)  # plogis is base R inverse logit
mean_probs <- colMeans(pred_probs)
actual <- as.numeric(model_data_final$RILD) - 1


waic_val <- NA
waic_se <- NA
loo_val <- NA
loo_se <- NA

tryCatch({
  # Extract posterior samples (brms)
  post_samples <- as_draws_df(model)
  post_samples_df <- as.data.frame(post_samples)
  # Get number of samples (rows in the draws data frame)
  n_samples <- nrow(post_samples_df)
  n_obs <- length(data_list_sanitized$RILD)
  
  # Get predictions for all observations across all samples (brms)
  # This gives us a matrix: rows = samples, cols = observations
  pred_linpred_all <- posterior_linpred(model, transform = FALSE)
  pred_probs_all <- plogis(pred_linpred_all)  # plogis is base R inverse logit
  
  # Calculate log-likelihood for each observation across all samples
  actual_vec <- as.numeric(data_list_sanitized$RILD)
  
  # Log-likelihood matrix: samples x observations (required by loo package)
  log_lik_matrix <- matrix(NA, nrow = nrow(pred_probs_all), ncol = ncol(pred_probs_all))
  
  for (i in 1:ncol(pred_probs_all)) {
    # For each observation, calculate log-likelihood across all samples
    probs_i <- pred_probs_all[, i]
    y_i <- actual_vec[i]
    
    # Log-likelihood: y*log(p) + (1-y)*log(1-p)
    # Handle edge cases where p is 0 or 1
    probs_i <- pmax(probs_i, 1e-10)
    probs_i <- pmin(probs_i, 1 - 1e-10)
    
    log_lik_matrix[, i] <- y_i * log(probs_i) + (1 - y_i) * log(1 - probs_i)
  }
  
  # Use loo package to calculate WAIC and LOO
  cat("  Computing WAIC and LOO from log-likelihood matrix...\n")
  
  # Calculate WAIC using loo package
  waic_result <- waic(log_lik_matrix)
  waic_val <- waic_result$estimates["waic", "Estimate"]
  waic_se <- waic_result$estimates["waic", "SE"]
  
  cat("  WAIC:", round(waic_val, 2), "±", round(waic_se, 2), "\n")
  cat("  Effective number of parameters (p_waic):", 
      round(waic_result$estimates["p_waic", "Estimate"], 2), "\n")
  
  # Calculate LOO using loo package
  loo_result <- loo(log_lik_matrix)
  loo_val <- loo_result$estimates["looic", "Estimate"]
  loo_se <- loo_result$estimates["looic", "SE"]
  
  cat("  LOO:", round(loo_val, 2), "±", round(loo_se, 2), "\n")
  
  # Check for problematic observations
  if (any(loo_result$pareto_k > 0.7)) {
    n_problematic <- sum(loo_result$pareto_k > 0.7)
    cat("  Warning:", n_problematic, "observations with high Pareto k (>0.7)\n")
    cat("  These may have high influence on LOO estimates\n")
  }
  
}, error = function(e) {
  cat("  Error calculating WAIC/LOO:", e$message, "\n")
  cat("  Attempting fallback calculation...\n")
  
  # Fallback: manual WAIC calculation
  tryCatch({
    pred_linpred_all <- posterior_linpred(model, transform = FALSE)
    pred_probs_all <- plogis(pred_linpred_all)
    actual_vec <- as.numeric(data_list_sanitized$RILD)
    
    log_lik_matrix <- matrix(NA, nrow = nrow(pred_probs_all), ncol = ncol(pred_probs_all))
    for (i in 1:ncol(pred_probs_all)) {
      probs_i <- pmax(pmin(pred_probs_all[, i], 1 - 1e-10), 1e-10)
      log_lik_matrix[, i] <- actual_vec[i] * log(probs_i) + (1 - actual_vec[i]) * log(1 - probs_i)
    }
    
    pointwise_log_lik <- colMeans(log_lik_matrix)
    pointwise_var_log_lik <- apply(log_lik_matrix, 2, var)
    lppd <- sum(pointwise_log_lik)
    p_waic <- sum(pointwise_var_log_lik)
    waic_val <<- -2 * (lppd - p_waic)
    
    cat("  WAIC (fallback):", round(waic_val, 2), "\n")
  }, error = function(e2) {
    cat("  Could not calculate WAIC:", e2$message, "\n")
  })
})


# Create initial comparison table with full model metrics
full_model_summary <- data.frame(
  Model = "Full Model (All Predictors)",
  Number_Predictors = as.integer(length(predictor_names_sanitized)),
  WAIC = if (exists("waic_val") && !is.na(waic_val)) as.numeric(waic_val)[1] else NA_real_,
  WAIC_SE = if (exists("waic_se") && !is.na(waic_se)) as.numeric(waic_se)[1] else NA_real_,
  LOO = if (exists("loo_val") && !is.na(loo_val)) as.numeric(loo_val)[1] else NA_real_,
  LOO_SE = if (exists("loo_se") && !is.na(loo_se)) as.numeric(loo_se)[1] else NA_real_,
  AUC = if (exists("auc_val") && !is.na(auc_val)) as.numeric(auc_val)[1] else NA_real_,
  stringsAsFactors = FALSE
)

# Save full model summary
write.csv(full_model_summary, "full_model_summary.csv", row.names = FALSE)
cat("Full model summary saved to: full_model_summary.csv\n\n")

# ROC curve and AUC
cat("\nCalculating AUC (Area Under ROC Curve)...\n")
roc_obj <- roc(actual, mean_probs, quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))
cat("  AUC:", round(auc_val, 4), "\n")

# Optimal threshold
optimal_threshold <- tryCatch({
  coords(roc_obj, "best", ret = "threshold", best.method = "youden")$threshold
}, error = function(e) {
  mean(actual)
})

if (is.na(optimal_threshold) || optimal_threshold <= 0 || optimal_threshold >= 1) {
  optimal_threshold <- 0.5
}

cat("Optimal threshold:", round(optimal_threshold, 4), "\n\n")

# Confusion matrix
pred_class <- ifelse(mean_probs > optimal_threshold, 1, 0)
cm <- table(Actual = actual, Predicted = pred_class)
print(cm)

if (nrow(cm) == 2 && ncol(cm) == 2) {
  sens <- cm[2,2] / sum(cm[2,])
  spec <- cm[1,1] / sum(cm[1,])
  ppv <- cm[2,2] / sum(cm[,2])
  npv <- cm[1,1] / sum(cm[,1])
  }

# Update and save full model summary metrics with all performance metrics
full_model_summary_metrics <- data.frame(
  Model = "Full Model (All Predictors)",
  Number_Predictors = as.integer(length(predictor_names_sanitized)),
  WAIC = if (exists("waic_val") && !is.na(waic_val)) as.numeric(waic_val)[1] else NA_real_,
  WAIC_SE = if (exists("waic_se") && !is.na(waic_se)) as.numeric(waic_se)[1] else NA_real_,
  LOO = if (exists("loo_val") && !is.na(loo_val)) as.numeric(loo_val)[1] else NA_real_,
  LOO_SE = if (exists("loo_se") && !is.na(loo_se)) as.numeric(loo_se)[1] else NA_real_,
  AUC = if (exists("auc_val") && !is.na(auc_val)) as.numeric(auc_val)[1] else NA_real_,
  Sensitivity = if (exists("sens") && !is.na(sens)) as.numeric(sens)[1] else NA_real_,
  Specificity = if (exists("spec") && !is.na(spec)) as.numeric(spec)[1] else NA_real_,
  PPV = if (exists("ppv") && !is.na(ppv)) as.numeric(ppv)[1] else NA_real_,
  NPV = if (exists("npv") && !is.na(npv)) as.numeric(npv)[1] else NA_real_,
  Accuracy = if (exists("actual") && exists("pred_class")) as.numeric(round(mean(actual == pred_class), 4)) else NA_real_,
  stringsAsFactors = FALSE
)
write.csv(full_model_summary_metrics, "full_model_summary_metrics.csv", row.names = FALSE)


# Save summary to file
summary_metrics <- c(
  AUC = round(auc_val, 4),
  WAIC = ifelse(is.na(waic_val), NA, round(waic_val, 2)),
  WAIC_SE = ifelse(is.na(waic_se), NA, round(waic_se, 2)),
  Sensitivity = ifelse(exists("sens"), round(sens, 4), NA),
  Specificity = ifelse(exists("spec"), round(spec, 4), NA),
  PPV = ifelse(exists("ppv"), round(ppv, 4), NA),
  NPV = ifelse(exists("npv"), round(npv, 4), NA),
  Accuracy = round(mean(actual == pred_class), 4),
  N_Predictors = length(predictor_names),
  N_Samples = nrow(model_data_final)
)

summary_results <- data.frame(
  Metric = names(summary_metrics),
  Value = as.numeric(summary_metrics)
)

write.csv(summary_results, "model_summary_metrics.csv", row.names = FALSE)
