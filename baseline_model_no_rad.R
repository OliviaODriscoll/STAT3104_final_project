# Test Script: Baseline Model Without Radiomics
# This script tests just the baseline model with diagnostics

# ============================================================================
# 1. SETUP AND DATA LOADING
# ============================================================================

library(tidyverse)
library(rethinking)    # Bayesian regression models using Stan (ulam)
library(bayesplot)     # MCMC diagnostics and posterior plots
library(mice)          # Multiple imputation
library(pROC)          # ROC curves
library(glmnet)        # Lasso/Elastic Net for feature selection

set.seed(42)

cat("=== TESTING BASELINE MODEL (NO RADIOMICS) ===\n\n")

# Load data files
cat("=== LOADING DATA ===\n\n")
baseline_labs <- read.csv("data/baseline_labs_anonymized.csv", stringsAsFactors = FALSE)
clinical_features <- read.csv("data/Clinical_Features_anonymized.csv", stringsAsFactors = FALSE)
dvh_data <- read.csv("data/dose_volume_histogram_5gy_bins_anonymized.csv", stringsAsFactors = FALSE)
rild_outcome <- read.csv("data/patients_to_process.csv", stringsAsFactors = FALSE)

# ============================================================================
# 2. DATA CLEANING
# ============================================================================

cat("=== DATA CLEANING ===\n\n")

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

# ============================================================================
# 3. FEATURE ENGINEERING
# ============================================================================

cat("=== FEATURE ENGINEERING ===\n\n")

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

v30_cols <- c("30-35Gy", "35-40Gy", "40-45Gy", "45-50Gy", 
              "50-55Gy", "55-60Gy", "60-65Gy", "65-70Gy", 
              "70-75Gy", "75-80Gy", "80-85Gy", "85-90Gy", 
              "90-95Gy", "95-100Gy")
v30_cols_exist <- v30_cols[v30_cols %in% names(merged_base)]

v20_cols <- c("20-25Gy", "25-30Gy", "30-35Gy", "35-40Gy", 
              "40-45Gy", "45-50Gy", "50-55Gy", "55-60Gy", 
              "60-65Gy", "65-70Gy", "70-75Gy", "75-80Gy", 
              "80-85Gy", "85-90Gy", "90-95Gy", "95-100Gy")
v20_cols_exist <- v20_cols[v20_cols %in% names(merged_base)]

V30_calc <- if (length(v30_cols_exist) > 0) {
  rowSums(merged_base[, v30_cols_exist, drop = FALSE], na.rm = TRUE)
} else {
  rep(NA, nrow(merged_base))
}

V20_calc <- if (length(v20_cols_exist) > 0) {
  rowSums(merged_base[, v20_cols_exist, drop = FALSE], na.rm = TRUE)
} else {
  rep(NA, nrow(merged_base))
}

merged_base <- merged_base %>%
  mutate(
    mean_dose = mean_dose_calc,
    V30 = V30_calc,
    V20 = V20_calc,
    liver_volume_ratio = Liver_minus_PTV_volume_cc / 
                        (PTV_volume_cc + Liver_minus_PTV_volume_cc),
    ast_alt_ratio = ast / alt
  ) %>%
  mutate_all(~ifelse(is.infinite(.), NA, .))

cat("Base dataset:", nrow(merged_base), "patients\n")
cat("RILD cases:", sum(merged_base$RILD == 1, na.rm = TRUE), "\n")
cat("RILD controls:", sum(merged_base$RILD == 0, na.rm = TRUE), "\n\n")

# ============================================================================
# 4. HELPER FUNCTIONS
# ============================================================================

# One-hot encode function
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
        for (i in 2:n_levels) {
          level_name <- levels[i]
          var_name <- paste0(var, "_", gsub("[^A-Za-z0-9_]", "_", level_name))
          data_encoded[[var_name]] <- as.numeric(data_encoded[[var]] == level_name)
        }
        data_encoded[[var]] <- NULL
      }
    }
  }
  return(data_encoded)
}

# ============================================================================
# 5. PREPARE DATA (SIMPLIFIED VERSION)
# ============================================================================

cat("=== PREPARING DATA ===\n\n")

# Use a subset of the prepare_and_standardize logic, simplified
# Select key variables
key_vars <- c("Age", "Gender", "Diagnosis", "platelet", "alp", "bilirubin", 
              "albumin", "ast", "alt", "inr", "PTV_volume_cc", 
              "Liver_minus_PTV_volume_cc", "EUD_alpha_0.06", "mean_dose", 
              "V30", "V20", "liver_volume_ratio", "ast_alt_ratio")

key_vars <- key_vars[key_vars %in% names(merged_base)]

# Prepare data
model_data <- merged_base %>%
  select(MRN, RILD, all_of(key_vars)) %>%
  filter(!is.na(RILD))

# Simple imputation (median for numeric, mode for factors)
for (var in key_vars) {
  if (var %in% names(model_data)) {
    if (is.numeric(model_data[[var]])) {
      model_data[[var]][is.na(model_data[[var]])] <- median(model_data[[var]], na.rm = TRUE)
    } else if (is.factor(model_data[[var]])) {
      most_common <- names(sort(table(model_data[[var]]), decreasing = TRUE))[1]
      if (!is.null(most_common)) {
        model_data[[var]][is.na(model_data[[var]])] <- most_common
      }
    }
  }
}

# Standardize numeric features
numeric_vars <- model_data %>%
  select(-MRN, -RILD) %>%
  select_if(is.numeric) %>%
  names()

for (var in numeric_vars) {
  if (is.numeric(model_data[[var]])) {
    var_mean <- mean(model_data[[var]], na.rm = TRUE)
    var_sd <- sd(model_data[[var]], na.rm = TRUE)
    if (!is.na(var_sd) && var_sd > 1e-10) {
      model_data[[var]] <- (model_data[[var]] - var_mean) / var_sd
    } else {
      model_data[[var]] <- 0
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
  cat("WARNING: Found radiomics variables in data:\n")
  print(radiomics_vars)
  cat("Removing radiomics variables...\n")
  model_data_final <- model_data_final %>%
    select(-all_of(radiomics_vars))
} else {
  cat("✓ Confirmed: No radiomics variables in baseline model\n")
}

cat("\nFinal data:", nrow(model_data_final), "patients,", 
    ncol(model_data_final) - 1, "predictors\n")
cat("RILD: 0 =", sum(model_data_final$RILD == 0), 
    ", 1 =", sum(model_data_final$RILD == 1), "\n")
cat("Variables:", paste(setdiff(names(model_data_final), "RILD"), collapse = ", "), "\n\n")

# ============================================================================
# 6. TRAIN MODEL
# ============================================================================

cat("=== TRAINING BASELINE MODEL ===\n\n")

# Get predictor names
predictor_names <- setdiff(names(model_data_final), "RILD")
cat("Predictors:", length(predictor_names), "\n")

# Prepare data for ulam
n_obs <- nrow(model_data_final)
data_list <- list()
data_list$RILD <- as.numeric(model_data_final$RILD) - 1

for (col in predictor_names) {
  if (is.factor(model_data_final[[col]])) {
    data_list[[col]] <- as.numeric(model_data_final[[col]]) - 1
  } else {
    data_list[[col]] <- as.numeric(model_data_final[[col]])
  }
  if (any(is.na(data_list[[col]]))) {
    data_list[[col]][is.na(data_list[[col]])] <- 0
  }
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

# Build model formula
coef_terms <- paste0("b_", predictor_names_sanitized, " * ", predictor_names_sanitized)
coef_terms_str <- paste(coef_terms, collapse = " + ")

model_code <- paste0(
  "alist(\n",
  "  RILD ~ bernoulli_logit(mu),\n",
  "  mu <- a + ", coef_terms_str, ",\n",
  "  a ~ normal(", intercept_mean, ", 2.5)"
)

for (pred in predictor_names_sanitized) {
  model_code <- paste0(model_code, ",\n  b_", pred, " ~ normal(0, 1)")
}
model_code <- paste0(model_code, "\n)")

model_formula <- eval(parse(text = model_code))

cat("Fitting model with ulam...\n")
cat("This may take 5-15 minutes...\n\n")

# Fit model
model <- ulam(
  model_formula,
  data = data_list_sanitized,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  seed = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

cat("\n=== MODEL FITTING COMPLETE ===\n\n")

# ============================================================================
# 7. MCMC DIAGNOSTICS
# ============================================================================

cat("=== MCMC DIAGNOSTICS ===\n\n")

# Extract samples for diagnostics
post <- extract.samples(model)
post_df <- as.data.frame(post)

# Trace plots for key parameters
cat("Generating trace plots...\n")

# Get parameter names
param_names <- names(post_df)
key_params <- c("a", grep("^b_", param_names, value = TRUE)[1:min(10, length(grep("^b_", param_names)))])
key_params <- key_params[key_params %in% param_names]

cat("Plotting trace plots for:", length(key_params), "parameters\n")

# Create trace plots from extracted samples
# extract.samples() combines all chains, so we'll create trace plots showing the combined samples
cat("  Creating trace plots from extracted samples...\n")

# Get total number of samples (4 chains * 1000 post-warmup iterations)
n_chains <- 4
n_warmup <- 1000
n_post_warmup <- 1000
total_samples <- n_chains * n_post_warmup

# Create trace plots
png("baseline_model_traceplots.png", width = 2000, height = 1200, res = 150)
n_params <- length(key_params)
n_cols <- 2
n_rows <- ceiling(n_params / n_cols)
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 2))

for (param in key_params) {
  if (param %in% names(post_df)) {
    param_samples <- post_df[[param]]
    
    # Create trace plot
    plot(1:length(param_samples), param_samples, type = "l", 
         main = paste("Trace:", param),
         xlab = "Sample", ylab = "Value", 
         col = "steelblue", lwd = 0.5)
    
    # Add mean line
    abline(h = mean(param_samples), col = "red", lty = 2, lwd = 2)
    
    # Add credible interval
    ci_lower <- quantile(param_samples, 0.025)
    ci_upper <- quantile(param_samples, 0.975)
    abline(h = c(ci_lower, ci_upper), col = "orange", lty = 3, lwd = 1.5)
    
    # Add text with summary stats
    legend_text <- sprintf("Mean=%.3f\nSD=%.3f", 
                          mean(param_samples), sd(param_samples))
    legend("topright", legend = legend_text, bty = "n", cex = 0.8)
  }
}
dev.off()

cat("Trace plots saved to: baseline_model_traceplots.png\n\n")

# Also create density plots
png("baseline_model_density_plots.png", width = 2000, height = 1200, res = 150)
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 2))

for (param in key_params) {
  if (param %in% names(post_df)) {
    param_samples <- post_df[[param]]
    
    # Create density plot
    dens <- density(param_samples)
    plot(dens, main = paste("Posterior Density:", param),
         xlab = "Value", ylab = "Density", 
         col = "steelblue", lwd = 2)
    
    # Add credible interval
    ci_lower <- quantile(param_samples, 0.025)
    ci_upper <- quantile(param_samples, 0.975)
    abline(v = c(ci_lower, ci_upper), col = "red", lty = 2, lwd = 2)
    abline(v = mean(param_samples), col = "darkgreen", lty = 1, lwd = 2)
    
    # Add text
    legend_text <- sprintf("Mean=%.3f\n95%% CI: [%.3f, %.3f]", 
                          mean(param_samples), ci_lower, ci_upper)
    legend("topright", legend = legend_text, bty = "n", cex = 0.8)
  }
}
dev.off()

cat("Density plots saved to: baseline_model_density_plots.png\n\n")

# R-hat and effective sample size using precis from rethinking
cat("Convergence diagnostics (using precis):\n")
precis_output <- precis(model, depth = 2)
print(precis_output)

# Try to extract R-hat from precis if available
if (!is.null(precis_output) && "Rhat4" %in% names(precis_output)) {
  rhat_vals <- precis_output$Rhat4
  cat("\nR-hat summary:\n")
  cat("  Max R-hat:", round(max(rhat_vals, na.rm = TRUE), 4), "\n")
  cat("  Min R-hat:", round(min(rhat_vals, na.rm = TRUE), 4), "\n")
  high_rhat <- rhat_vals[rhat_vals > 1.01]
  if (length(high_rhat) > 0) {
    cat("  WARNING:", length(high_rhat), "parameters with R-hat > 1.01\n")
  } else {
    cat("  ✓ All R-hat values < 1.01 (good convergence)\n")
  }
} else {
  cat("\nR-hat values are in the precis output above.\n")
}

# Effective sample size
if (!is.null(precis_output) && "n_eff" %in% names(precis_output)) {
  neff_vals <- precis_output$n_eff
  cat("\nEffective sample size:\n")
  cat("  Min Neff:", round(min(neff_vals, na.rm = TRUE), 0), "\n")
  cat("  Max Neff:", round(max(neff_vals, na.rm = TRUE), 0), "\n")
  low_neff <- neff_vals[neff_vals < 400]  # 400 = 0.1 * 4000 (4 chains * 1000 samples)
  if (length(low_neff) > 0) {
    cat("  WARNING:", length(low_neff), "parameters with low effective sample size (< 400)\n")
  } else {
    cat("  ✓ All Neff values > 400 (good sampling)\n")
  }
}

# ============================================================================
# 8. POSTERIOR SUMMARIES
# ============================================================================

cat("\n=== POSTERIOR SUMMARIES ===\n\n")

# Summary of intercept
if ("a" %in% names(post_df)) {
  cat("Intercept (a):\n")
  cat("  Mean:", round(mean(post_df$a), 4), "\n")
  cat("  SD:", round(sd(post_df$a), 4), "\n")
  cat("  2.5%:", round(quantile(post_df$a, 0.025), 4), "\n")
  cat("  97.5%:", round(quantile(post_df$a, 0.975), 4), "\n\n")
}

# Top coefficients by absolute mean
coef_means <- sapply(grep("^b_", names(post_df), value = TRUE), function(x) {
  mean(abs(post_df[[x]]))
})
coef_means <- sort(coef_means, decreasing = TRUE)

cat("Top 5 coefficients by absolute mean:\n")
for (i in 1:min(5, length(coef_means))) {
  coef_name <- names(coef_means)[i]
  orig_name <- name_mapping$original[name_mapping$sanitized == gsub("^b_", "", coef_name)]
  cat(sprintf("  %s (%s): %.4f\n", coef_name, orig_name, coef_means[i]))
}

# ============================================================================
# 9. PREDICTIONS AND PERFORMANCE
# ============================================================================

cat("\n=== PREDICTIONS AND PERFORMANCE ===\n\n")

# Get predictions
pred_linpred <- link(model, data = data_list_sanitized)
pred_probs <- inv_logit(pred_linpred)
mean_probs <- colMeans(pred_probs)
actual <- as.numeric(model_data_final$RILD) - 1

cat("Prediction summary:\n")
cat("  Min:", round(min(mean_probs), 4), "\n")
cat("  Max:", round(max(mean_probs), 4), "\n")
cat("  Mean:", round(mean(mean_probs), 4), "\n")
cat("  Median:", round(median(mean_probs), 4), "\n\n")

cat("Actual outcomes:\n")
cat("  0 =", sum(actual == 0), "\n")
cat("  1 =", sum(actual == 1), "\n\n")

# ROC curve
roc_obj <- roc(actual, mean_probs, quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))
cat("AUC:", round(auc_val, 4), "\n")

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
  
  cat("\nPerformance metrics:\n")
  cat("  Sensitivity:", round(sens, 4), "\n")
  cat("  Specificity:", round(spec, 4), "\n")
  cat("  PPV:", round(ppv, 4), "\n")
  cat("  NPV:", round(npv, 4), "\n")
  cat("  Accuracy:", round(mean(actual == pred_class), 4), "\n")
}

cat("\n=== TEST COMPLETE ===\n")

