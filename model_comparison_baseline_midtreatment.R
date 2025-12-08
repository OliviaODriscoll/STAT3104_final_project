# Bayesian Model Comparison: Baseline vs Mid-Treatment, With and Without Radiomics
# Applied Bayesian Analysis Course Project

# ============================================================================
# 1. SETUP AND DATA LOADING
# ============================================================================

# Load required libraries
library(tidyverse)
library(brms)          # Bayesian regression models using Stan
library(bayesplot)     # MCMC diagnostics and posterior plots
library(mice)          # Multiple imputation
library(pROC)          # ROC curves
library(loo)           # Leave-one-out cross-validation

# Set seed for reproducibility
set.seed(42)

cat("=== MODEL COMPARISON: BASELINE vs MID-TREATMENT ===\n")
cat("Training 4 models:\n")
cat("1. Baseline without radiomics\n")
cat("2. Baseline with radiomics (Primary)\n")
cat("3. Mid-treatment without radiomics\n")
cat("4. Mid-treatment with radiomics (CTORMid)\n\n")

# Load data files
cat("=== LOADING DATA ===\n\n")
baseline_labs <- read.csv("data/baseline_labs_anonymized.csv", stringsAsFactors = FALSE)
clinical_features <- read.csv("data/Clinical_Features_anonymized.csv", stringsAsFactors = FALSE)
dvh_data <- read.csv("data/dose_volume_histogram_5gy_bins_anonymized.csv", stringsAsFactors = FALSE)
rild_outcome <- read.csv("data/patients_to_process.csv", stringsAsFactors = FALSE)

# Load radiomics data
radiomics_ctormid <- read.csv("data/CTORMid_Binned_Radiomics_anonymized.csv", stringsAsFactors = FALSE)
radiomics_primary <- read.csv("data/Primary_Radiomics_bin25_anonymized.csv", stringsAsFactors = FALSE)

cat("CTORMid radiomics (mid-treatment):", ncol(radiomics_ctormid) - 1, "features\n")
cat("Primary radiomics (baseline):", ncol(radiomics_primary) - 1, "features\n\n")

# ============================================================================
# 2. DATA CLEANING AND PREPROCESSING
# ============================================================================

cat("=== DATA CLEANING ===\n\n")

# Clean baseline labs
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

# Clean clinical features
clinical_features_clean <- clinical_features %>%
  mutate(
    Age = as.numeric(Age),
    Gender = factor(Gender),
    Diagnosis = factor(Diagnosis),
    PTV_volume_cc = as.numeric(PTV_volume_cc),
    Liver_minus_PTV_volume_cc = as.numeric(Liver_minus_PTV_volume_cc),
    EUD_alpha_0.06 = as.numeric(EUD_alpha_0.06)
  )

# Clean DVH data
dvh_clean <- dvh_data %>%
  mutate_at(vars(-MRN), as.numeric)

# Clean radiomics
# Clean radiomics data - convert all to numeric (except MRN)
# Also remove duplicate MRNs to avoid many-to-many join issues
radiomics_ctormid_clean <- radiomics_ctormid %>%
  mutate_at(vars(-MRN), as.numeric) %>%
  distinct(MRN, .keep_all = TRUE)  # Remove duplicate MRNs

radiomics_primary_clean <- radiomics_primary %>%
  mutate_at(vars(-MRN), as.numeric) %>%
  distinct(MRN, .keep_all = TRUE)  # Remove duplicate MRNs

cat("CTORMid radiomics after removing duplicates:", nrow(radiomics_ctormid_clean), "patients\n")
cat("Primary radiomics after removing duplicates:", nrow(radiomics_primary_clean), "patients\n\n")

# Clean outcome
rild_clean <- rild_outcome %>%
  mutate(RILD = as.numeric(RILD))

# ============================================================================
# 3. FEATURE ENGINEERING
# ============================================================================

cat("=== FEATURE ENGINEERING ===\n\n")

# Merge clinical/lab/dosimetric data
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

# Calculate V30 and V20
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
cat("RILD cases:", sum(merged_base$RILD == 1, na.rm = TRUE), "\n\n")

# ============================================================================
# 4. HELPER FUNCTION: PREPARE AND STANDARDIZE DATA
# ============================================================================

# Helper function to remove highly correlated features
remove_collinear_features <- function(data, cor_threshold = 0.95, exclude_vars = c("MRN", "RILD", "Gender", "Diagnosis")) {
  cat("Removing highly correlated features (threshold:", cor_threshold, ")...\n")
  
  # Get numeric features only
  numeric_vars <- data %>%
    select(-all_of(exclude_vars[exclude_vars %in% names(data)])) %>%
    select_if(is.numeric) %>%
    names()
  
  if (length(numeric_vars) < 2) {
    cat("Not enough numeric features for correlation analysis.\n")
    return(data)
  }
  
  # Calculate correlation matrix
  numeric_data <- data %>%
    select(all_of(numeric_vars)) %>%
    mutate_all(~ifelse(is.na(.), 0, .))  # Handle NAs temporarily
  
  # Check for constant or near-constant variables
  var_sds <- sapply(numeric_data, sd, na.rm = TRUE)
  constant_vars <- names(var_sds[var_sds < 1e-10 | is.na(var_sds)])
  if (length(constant_vars) > 0) {
    cat("Removing", length(constant_vars), "constant/near-constant variables\n")
    numeric_vars <- setdiff(numeric_vars, constant_vars)
    numeric_data <- numeric_data %>% select(-all_of(constant_vars))
  }
  
  if (length(numeric_vars) < 2) {
    cat("Not enough variable features after removing constants.\n")
    return(data)
  }
  
  # Calculate correlation matrix
  cor_matrix <- tryCatch({
    cor(numeric_data, use = "pairwise.complete.obs")
  }, error = function(e) {
    cat("Error calculating correlations:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(cor_matrix)) {
    cat("Could not calculate correlation matrix. Skipping collinearity filtering.\n")
    return(data)
  }
  
  # Replace NA/NaN/Inf with 0
  cor_matrix[is.na(cor_matrix) | is.nan(cor_matrix) | is.infinite(cor_matrix)] <- 0
  
  # Find highly correlated pairs
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- 0  # Set lower triangle and diagonal to 0
  high_cor_pairs <- which(abs(cor_matrix) > cor_threshold, arr.ind = TRUE)
  
  if (nrow(high_cor_pairs) == 0) {
    cat("No highly correlated features found.\n")
    return(data)
  }
  
  # Identify features to remove (keep first feature in each pair, remove second)
  vars_to_remove <- character(0)
  vars_to_keep <- character(0)
  
  for (i in 1:nrow(high_cor_pairs)) {
    var1 <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
    var2 <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
    
    # Keep the one that appears first in the list (or has lower index)
    # Remove the other one
    if (!var1 %in% vars_to_keep && !var1 %in% vars_to_remove) {
      vars_to_keep <- c(vars_to_keep, var1)
    }
    if (!var2 %in% vars_to_keep) {
      vars_to_remove <- c(vars_to_remove, var2)
    }
  }
  
  vars_to_remove <- unique(vars_to_remove)
  
  cat("Found", nrow(high_cor_pairs), "highly correlated pairs\n")
  cat("Removing", length(vars_to_remove), "features due to high collinearity\n")
  cat("Keeping", length(vars_to_keep), "representative features\n")
  
  # Remove highly correlated features
  if (length(vars_to_remove) > 0) {
    data <- data %>% select(-all_of(vars_to_remove))
  }
  
  return(data)
}

prepare_and_standardize <- function(base_data, radiomics_data = NULL, model_name, remove_collinear = TRUE, cor_threshold = 0.95) {
  cat("\n=== Preparing data for", model_name, "===\n")
  
  # Merge with radiomics if provided
  if (!is.null(radiomics_data)) {
    # Check for duplicate MRNs and handle them
    cat("Checking for duplicate MRNs in radiomics data...\n")
    duplicate_mrns_radio <- radiomics_data$MRN[duplicated(radiomics_data$MRN)]
    if (length(duplicate_mrns_radio) > 0) {
      cat("WARNING: Found", length(duplicate_mrns_radio), "duplicate MRNs in radiomics data\n")
      cat("Keeping first occurrence of each MRN\n")
      radiomics_data <- radiomics_data %>%
        distinct(MRN, .keep_all = TRUE)
    }
    
    duplicate_mrns_base <- base_data$MRN[duplicated(base_data$MRN)]
    if (length(duplicate_mrns_base) > 0) {
      cat("WARNING: Found", length(duplicate_mrns_base), "duplicate MRNs in base data\n")
      cat("Keeping first occurrence of each MRN\n")
      base_data <- base_data %>%
        distinct(MRN, .keep_all = TRUE)
    }
    
    # Verify all base patients are in radiomics before joining
    base_mrns_set <- sort(unique(base_data$MRN[!is.na(base_data$RILD)]))
    radio_mrns_set <- sort(unique(radiomics_data$MRN))
    
    missing_patients <- setdiff(base_mrns_set, radio_mrns_set)
    if (length(missing_patients) > 0) {
      cat("ERROR:", length(missing_patients), "patients in base data not in radiomics!\n")
      cat("This should not happen if datasets were properly filtered.\n")
      cat("First 10 missing:", paste(head(missing_patients, 10), collapse = ", "), "\n")
      stop("Cannot proceed: patients missing from radiomics data")
    }
    
    # Use inner_join - should keep all patients since we verified they all exist
    data_merged <- base_data %>%
      inner_join(radiomics_data, by = "MRN", relationship = "one-to-one")
    
    merged_mrns_actual <- sort(unique(data_merged$MRN[!is.na(data_merged$RILD)]))
    cat("Merged with radiomics:", ncol(radiomics_data) - 1, "features\n")
    cat("Patients after merge:", length(merged_mrns_actual), "\n")
    
    if (length(merged_mrns_actual) != length(base_mrns_set)) {
      cat("ERROR: Lost patients during join!\n")
      cat("Expected:", length(base_mrns_set), "Got:", length(merged_mrns_actual), "\n")
      lost <- setdiff(base_mrns_set, merged_mrns_actual)
      cat("Lost patients:", paste(head(lost, 10), collapse = ", "), "\n")
      stop("Patients lost during radiomics merge")
    }
  } else {
    data_merged <- base_data
    cat("No radiomics included\n")
    cat("Patients:", nrow(data_merged), "\n")
  }
  
  # Store original patient count for verification
  original_n <- nrow(data_merged)
  
  # Remove highly correlated features if requested
  if (remove_collinear && !is.null(radiomics_data)) {
    cat("\nRemoving collinear features...\n")
    data_merged <- remove_collinear_features(data_merged, cor_threshold = cor_threshold)
    cat("Features after collinearity filtering:", ncol(data_merged) - 1, "\n")
  }
  
  # Select key clinical/lab/dosimetric variables for imputation
  key_vars <- c("Age", "Gender", "Diagnosis", "platelet", "alp", "bilirubin", 
                "albumin", "ast", "alt", "inr", "PTV_volume_cc", 
                "Liver_minus_PTV_volume_cc", "EUD_alpha_0.06", "mean_dose", 
                "V30", "V20", "liver_volume_ratio", "ast_alt_ratio")
  
  # Only keep variables that exist
  key_vars <- key_vars[key_vars %in% names(data_merged)]
  
  # Identify all numeric features (excluding MRN, RILD, Gender, Diagnosis)
  numeric_features <- data_merged %>%
    select(-MRN, -RILD) %>%
    select_if(is.numeric) %>%
    names()
  
  numeric_features <- setdiff(numeric_features, c("Gender", "Diagnosis"))
  
  # Prepare imputation data
  impute_vars <- intersect(key_vars, names(data_merged))
  impute_data <- data_merged %>%
    select(all_of(c(impute_vars, "RILD", "MRN")))
  
  # Bayesian Multiple Imputation using MICE
  cat("Performing Bayesian multiple imputation (MICE)...\n")
  cat("Variables to impute:", length(impute_vars), "\n")
  
  # Check missing data pattern
  missing_counts <- colSums(is.na(impute_data[, impute_vars, drop = FALSE]))
  if (sum(missing_counts) > 0) {
    cat("Missing data counts:\n")
    print(missing_counts[missing_counts > 0])
  } else {
    cat("No missing data in key variables.\n")
  }
  
  # Perform MICE imputation
  imputed_data_mice <- tryCatch({
    # Use MICE with Bayesian linear regression for numeric, logistic for binary
    mice_result <- mice(impute_data[, impute_vars, drop = FALSE], 
                       m = 5,                    # 5 imputations
                       method = "norm",          # Bayesian linear regression for numeric
                       maxit = 20,               # 20 iterations
                       seed = 42,
                       printFlag = FALSE)
    
    # Extract first imputation (for single imputation approach)
    # In a full analysis, you could pool across all 5 imputations
    data_imputed <- complete(mice_result, 1)
    
    cat("MICE imputation completed successfully.\n")
    list(imputed = data_imputed, mice_object = mice_result)
  }, error = function(e) {
    cat("WARNING: MICE imputation failed:", e$message, "\n")
    cat("Falling back to simple median/mode imputation...\n")
    
    # Fallback to simple imputation
    data_imputed_simple <- impute_data[, impute_vars, drop = FALSE]
    for (var in impute_vars) {
      if (is.numeric(data_imputed_simple[[var]])) {
        data_imputed_simple[[var]][is.na(data_imputed_simple[[var]])] <- 
          median(data_imputed_simple[[var]], na.rm = TRUE)
      } else if (is.factor(data_imputed_simple[[var]])) {
        most_common <- names(sort(table(data_imputed_simple[[var]]), decreasing = TRUE))[1]
        if (!is.null(most_common) && !is.na(most_common)) {
          data_imputed_simple[[var]][is.na(data_imputed_simple[[var]])] <- most_common
        }
      }
    }
    list(imputed = data_imputed_simple, mice_object = NULL)
  })
  
  # Get imputed data
  impute_data_imputed <- imputed_data_mice$imputed
  
  # Combine imputed variables with RILD and MRN
  # Make sure column order is preserved
  impute_data <- data.frame(
    MRN = impute_data$MRN,
    RILD = impute_data$RILD,
    impute_data_imputed,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Verify no NAs remain in imputed variables
  nas_remaining <- sum(is.na(impute_data[, impute_vars, drop = FALSE]))
  if (nas_remaining > 0) {
    cat("WARNING:", nas_remaining, "NAs still remain after MICE imputation\n")
    cat("These will be handled in the standardization step\n")
  } else {
    cat("✓ All key variables imputed successfully (no NAs remaining)\n")
  }
  
  cat("Imputation completed.\n")
  
  # Merge imputed key vars with radiomics
  # Note: radiomics are already in data_merged from the inner_join above
  if (!is.null(radiomics_data)) {
    # Get all variables from data_merged that are not in impute_data
    # These include radiomics and any other variables
    vars_to_add <- setdiff(names(data_merged), c(names(impute_data), "MRN"))
    
    if (length(vars_to_add) > 0) {
      cat("Adding", length(vars_to_add), "additional features (including radiomics)...\n")
      model_data_full <- impute_data %>%
        left_join(
          data_merged %>%
            select(MRN, all_of(vars_to_add)),
          by = "MRN"
        )
    } else {
      cat("No additional variables to add. Using imputed data only.\n")
      model_data_full <- impute_data
    }
  } else {
    model_data_full <- impute_data
  }
  
  # Separate factors and numeric
  factor_vars <- c("Gender", "Diagnosis")[c("Gender", "Diagnosis") %in% names(model_data_full)]
  numeric_vars <- setdiff(names(model_data_full), c("MRN", "RILD", factor_vars))
  
  # Standardize numeric features
  cat("Standardizing", length(numeric_vars), "numeric features...\n")
  
  for (var in numeric_vars) {
    if (is.numeric(model_data_full[[var]])) {
      # Check for missing values before standardization
      n_missing_before <- sum(is.na(model_data_full[[var]]))
      
      var_sd <- sd(model_data_full[[var]], na.rm = TRUE)
      if (!is.na(var_sd) && var_sd > 1e-10) {
        var_mean <- mean(model_data_full[[var]], na.rm = TRUE)
        model_data_full[[var]] <- (model_data_full[[var]] - var_mean) / var_sd
      } else {
        model_data_full[[var]] <- 0
      }
      
      # Handle any remaining NAs (should be minimal after MICE, but safety check)
      n_missing_after <- sum(is.na(model_data_full[[var]]))
      if (n_missing_after > 0) {
        cat("  WARNING:", n_missing_after, "missing values in", var, "after standardization\n")
        cat("  Setting to 0 (mean in standardized space)\n")
        model_data_full[[var]][is.na(model_data_full[[var]])] <- 0
      }
    }
  }
  
  # Prepare final model data
  # IMPORTANT: Keep MRN temporarily to ensure we can track patients
  # Track patients before filtering
  patients_before_final_filter <- sort(unique(model_data_full$MRN[!is.na(model_data_full$RILD)]))
  cat("Patients before final filter:", length(patients_before_final_filter), "\n")
  
  # Check for duplicate MRNs and handle them
  duplicate_mrns_in_data <- model_data_full$MRN[duplicated(model_data_full$MRN)]
  if (length(duplicate_mrns_in_data) > 0) {
    cat("WARNING: Found", length(duplicate_mrns_in_data), "duplicate MRNs in model_data_full\n")
    cat("Removing duplicates (keeping first occurrence)\n")
    model_data_full <- model_data_full %>%
      distinct(MRN, .keep_all = TRUE)
  }
  
  model_data_final <- model_data_full %>%
    mutate(
      RILD = as.numeric(as.character(RILD)),
      RILD = ifelse(RILD %in% c(0, 1), RILD, NA)
    ) %>%
    filter(!is.na(RILD)) %>%
    filter(!is.na(MRN)) %>%  # Ensure MRN is not missing
    distinct(MRN, .keep_all = TRUE) %>%  # Ensure no duplicate MRNs
    arrange(MRN)  # Sort by MRN for consistency
  
  patients_after_final_filter <- sort(unique(model_data_final$MRN))
  cat("Patients after final filter:", length(patients_after_final_filter), "\n")
  
  # Verify we didn't lose patients unexpectedly
  if (length(patients_after_final_filter) < length(patients_before_final_filter)) {
    lost_patients <- setdiff(patients_before_final_filter, patients_after_final_filter)
    cat("ERROR: Lost", length(lost_patients), "patients during final filtering!\n")
    cat("Lost patient MRNs:", paste(head(lost_patients, 10), collapse = ", "), "\n")
    if (length(lost_patients) > 10) cat("... and", length(lost_patients) - 10, "more\n")
    cat("This should not happen - all patients should have valid RILD values.\n")
    stop("Patients lost during final filtering")
  }
  
  # Verify row count matches patient count (no duplicates)
  if (nrow(model_data_final) != length(patients_after_final_filter)) {
    cat("WARNING: Row count (", nrow(model_data_final), ") != unique patient count (", 
        length(patients_after_final_filter), ")\n")
    cat("This indicates duplicate rows. Removing duplicates...\n")
    model_data_final <- model_data_final %>%
      distinct(MRN, .keep_all = TRUE)
    cat("After removing duplicates:", nrow(model_data_final), "rows\n")
  }
  
  # Convert factors if they exist
  if ("Gender" %in% names(model_data_final)) {
    model_data_final$Gender <- factor(model_data_final$Gender)
  }
  if ("Diagnosis" %in% names(model_data_final)) {
    model_data_final$Diagnosis <- factor(model_data_final$Diagnosis)
  }
  
  # Store MRNs before removing (for verification)
  final_mrns <- model_data_final$MRN
  
  model_data_final <- model_data_final %>%
    mutate(RILD = factor(RILD, levels = c(0, 1))) %>%
    select(-MRN)
  
  # Store MRNs as attribute for later verification
  attr(model_data_final, "MRNs") <- final_mrns
  
  # Remove any remaining NAs
  model_data_final <- model_data_final %>%
    mutate_all(~ifelse(is.na(.), 0, .))
  
  # Final verification
  final_n <- nrow(model_data_final)
  if (final_n != original_n) {
    cat("WARNING: Row count changed from", original_n, "to", final_n, "during processing!\n")
  }
  
  cat("Final data:", final_n, "patients,", ncol(model_data_final) - 1, "predictors\n")
  
  return(model_data_final)
}

# ============================================================================
# 5. IDENTIFY COMMON PATIENTS ACROSS ALL DATASETS
# ============================================================================

cat("\n=== IDENTIFYING COMMON PATIENTS ===\n\n")

# Find patients that exist in all datasets
patients_base <- unique(merged_base$MRN[!is.na(merged_base$RILD)])
patients_primary <- unique(radiomics_primary_clean$MRN)
patients_ctormid <- unique(radiomics_ctormid_clean$MRN)

# For models without radiomics, we can use all base patients
# For models with radiomics, we need patients that have both base data AND radiomics
patients_with_primary <- intersect(patients_base, patients_primary)
patients_with_ctormid <- intersect(patients_base, patients_ctormid)

# For fair comparison, use patients that have ALL data (base + both radiomics)
# OR use separate sets: all base patients for models 1&3, and intersection for models 2&4
# Let's use the intersection approach for models with radiomics
common_patients_all <- intersect(intersect(patients_base, patients_primary), patients_ctormid)

cat("Patients in base dataset:", length(patients_base), "\n")
cat("Patients with Primary radiomics:", length(patients_primary), "\n")
cat("Patients with CTORMid radiomics:", length(patients_ctormid), "\n")
cat("Patients with Primary + base:", length(patients_with_primary), "\n")
cat("Patients with CTORMid + base:", length(patients_with_ctormid), "\n")
cat("Patients with ALL data:", length(common_patients_all), "\n\n")

# Decision: Use common patients for ALL models to ensure fair comparison
# This ensures all models have the same number of observations
if (length(common_patients_all) < 10) {
  cat("WARNING: Very few patients have all data. Using patients with base + at least one radiomics.\n")
  common_patients <- unique(c(patients_with_primary, patients_with_ctormid))
  cat("Using", length(common_patients), "patients with base + at least one radiomics dataset.\n\n")
} else {
  common_patients <- common_patients_all
  cat("Using", length(common_patients), "patients with complete data for all models.\n\n")
}

# Filter ALL datasets to common patients BEFORE any processing
merged_base_common <- merged_base %>%
  filter(MRN %in% common_patients) %>%
  filter(!is.na(RILD))  # Ensure RILD is not missing

radiomics_primary_common <- radiomics_primary_clean %>%
  filter(MRN %in% common_patients)

radiomics_ctormid_common <- radiomics_ctormid_clean %>%
  filter(MRN %in% common_patients)

cat("Filtered datasets to common patients:\n")
cat("- Base data:", nrow(merged_base_common), "patients\n")
cat("- Primary radiomics:", nrow(radiomics_primary_common), "patients\n")
cat("- CTORMid radiomics:", nrow(radiomics_ctormid_common), "patients\n\n")

# Verify all have the same patients
mrns_base <- sort(unique(merged_base_common$MRN))
mrns_primary <- sort(unique(radiomics_primary_common$MRN))
mrns_ctormid <- sort(unique(radiomics_ctormid_common$MRN))

if (!identical(mrns_base, mrns_primary) || !identical(mrns_base, mrns_ctormid)) {
  cat("WARNING: Patient lists don't match exactly. Finding intersection...\n")
  final_common_patients <- Reduce(intersect, list(mrns_base, mrns_primary, mrns_ctormid))
  cat("Final common patients:", length(final_common_patients), "\n")
  
  merged_base_common <- merged_base_common %>% filter(MRN %in% final_common_patients)
  radiomics_primary_common <- radiomics_primary_common %>% filter(MRN %in% final_common_patients)
  radiomics_ctormid_common <- radiomics_ctormid_common %>% filter(MRN %in% final_common_patients)
  
  cat("After final filtering:\n")
  cat("- Base data:", nrow(merged_base_common), "patients\n")
  cat("- Primary radiomics:", nrow(radiomics_primary_common), "patients\n")
  cat("- CTORMid radiomics:", nrow(radiomics_ctormid_common), "patients\n\n")
} else {
  cat("✓ All datasets have the same patients\n\n")
}

# ============================================================================
# 6. PREPARE DATA FOR ALL 4 MODELS (USING COMMON PATIENTS)
# ============================================================================

cat("\n=== PREPARING DATA FOR ALL MODELS (COMMON PATIENTS) ===\n\n")

# Ensure all datasets have exactly the same patients before processing
cat("\n=== FINAL PATIENT ALIGNMENT ===\n")
final_common_mrns <- Reduce(intersect, list(
  sort(unique(merged_base_common$MRN[!is.na(merged_base_common$RILD)])),
  sort(unique(radiomics_primary_common$MRN)),
  sort(unique(radiomics_ctormid_common$MRN))
))

cat("Final common patients for ALL models:", length(final_common_mrns), "\n\n")

# Filter all datasets to exact same patients
merged_base_final <- merged_base_common %>%
  filter(MRN %in% final_common_mrns) %>%
  filter(!is.na(RILD))

radiomics_primary_final <- radiomics_primary_common %>%
  filter(MRN %in% final_common_mrns)

radiomics_ctormid_final <- radiomics_ctormid_common %>%
  filter(MRN %in% final_common_mrns)

cat("After final alignment:\n")
cat("- Base data:", nrow(merged_base_final), "patients\n")
cat("- Primary radiomics:", nrow(radiomics_primary_final), "patients\n")
cat("- CTORMid radiomics:", nrow(radiomics_ctormid_final), "patients\n\n")

# Verify all have same patients
mrns_check <- list(
  sort(unique(merged_base_final$MRN)),
  sort(unique(radiomics_primary_final$MRN)),
  sort(unique(radiomics_ctormid_final$MRN))
)

if (!all(sapply(mrns_check[-1], function(x) identical(mrns_check[[1]], x)))) {
  stop("ERROR: Datasets still don't have the same patients after final filtering!")
} else {
  cat("✓ All datasets have exactly the same", length(mrns_check[[1]]), "patients\n\n")
}

# Store the expected patient list
expected_patients <- sort(final_common_mrns)
cat("Expected patients for all models:", length(expected_patients), "\n\n")

# Model 1: Baseline without radiomics
data_model1 <- prepare_and_standardize(merged_base_final, NULL, "Model 1: Baseline without radiomics", 
                                       remove_collinear = FALSE)  # No radiomics, so no need

# Model 2: Baseline with Primary radiomics
data_model2 <- prepare_and_standardize(merged_base_final, radiomics_primary_final, 
                                       "Model 2: Baseline with Primary radiomics",
                                       remove_collinear = TRUE, cor_threshold = 0.95)

# Model 3: Mid-treatment without radiomics (same clinical/lab/dosimetric as baseline)
data_model3 <- prepare_and_standardize(merged_base_final, NULL, "Model 3: Mid-treatment without radiomics",
                                       remove_collinear = FALSE)  # No radiomics, so no need

# Model 4: Mid-treatment with CTORMid radiomics
data_model4 <- prepare_and_standardize(merged_base_final, radiomics_ctormid_final, 
                                       "Model 4: Mid-treatment with CTORMid radiomics",
                                       remove_collinear = TRUE, cor_threshold = 0.95)

# ============================================================================
# FORCE ALL MODELS TO HAVE THE SAME PATIENTS
# ============================================================================

cat("\n=== FORCING ALL MODELS TO USE SAME PATIENTS ===\n\n")

# Get actual patients from each model
mrns1 <- attr(data_model1, "MRNs")
mrns2 <- attr(data_model2, "MRNs")
mrns3 <- attr(data_model3, "MRNs")
mrns4 <- attr(data_model4, "MRNs")

# Find intersection of all models
all_mrns <- list(mrns1, mrns2, mrns3, mrns4)
all_mrns <- all_mrns[!sapply(all_mrns, is.null)]
if (length(all_mrns) > 0) {
  final_patient_set <- Reduce(intersect, all_mrns)
  cat("Patients present in all models:", length(final_patient_set), "\n")
  
  if (length(final_patient_set) < length(expected_patients)) {
    cat("WARNING: Some expected patients are missing from models\n")
    missing <- setdiff(expected_patients, final_patient_set)
    cat("Missing patients:", length(missing), "\n")
  }
  
  # If we have a valid patient set, we need to ensure all models use it
  # But since we can't easily add patients back, we'll use the intersection
  if (length(final_patient_set) >= 10) {
    cat("Using intersection of", length(final_patient_set), "patients for all models\n")
    cat("(This may be fewer than expected if patients were lost during processing)\n\n")
  } else {
    stop("Too few common patients (", length(final_patient_set), "). Cannot proceed.")
  }
} else {
  stop("Could not determine patient sets from models. Cannot align models.")
}

# Verify all models have the same number of observations and same patients
cat("\n=== VERIFYING DATA CONSISTENCY ===\n\n")
cat("Model 1 observations:", nrow(data_model1), "\n")
cat("Model 2 observations:", nrow(data_model2), "\n")
cat("Model 3 observations:", nrow(data_model3), "\n")
cat("Model 4 observations:", nrow(data_model4), "\n")

n_obs <- c(nrow(data_model1), nrow(data_model2), nrow(data_model3), nrow(data_model4))

# Get MRNs from each model (stored as attributes)
mrns1 <- attr(data_model1, "MRNs")
mrns2 <- attr(data_model2, "MRNs")
mrns3 <- attr(data_model3, "MRNs")
mrns4 <- attr(data_model4, "MRNs")

# Check if models have same number of observations
if (length(unique(n_obs)) > 1) {
  cat("\nERROR: Models have different numbers of observations!\n")
  cat("Model observations:", paste(n_obs, collapse = ", "), "\n")
  
  # Get common patients
  if (!is.null(mrns1) && !is.null(mrns2) && !is.null(mrns3) && !is.null(mrns4)) {
    common_mrns <- Reduce(intersect, list(mrns1, mrns2, mrns3, mrns4))
    cat("Common patients across all models:", length(common_mrns), "\n")
    
    # Show which models have which patients
    cat("\nPatient counts by model:\n")
    cat("Model 1:", length(mrns1), "patients\n")
    cat("Model 2:", length(mrns2), "patients\n")
    cat("Model 3:", length(mrns3), "patients\n")
    cat("Model 4:", length(mrns4), "patients\n")
    
    # Show missing patients
    for (i in 1:4) {
      model_mrns <- list(mrns1, mrns2, mrns3, mrns4)[[i]]
      missing <- setdiff(common_mrns, model_mrns)
      if (length(missing) > 0) {
        cat("Model", i, "missing", length(missing), "patients:", paste(head(missing, 5), collapse = ", "), "\n")
      }
    }
    
    if (length(common_mrns) >= 10) {
      cat("\nAll models should have", length(expected_patients), "patients but have different counts.\n")
      cat("This indicates patients were lost during processing.\n")
      cat("Please check the data preparation steps above for warnings about lost patients.\n")
      stop("Models have different N. Check data preparation steps for where patients were lost.")
    } else {
      stop("Too few common patients (", length(common_mrns), "). Cannot proceed.")
    }
  } else {
    stop("Cannot verify patient consistency. MRNs not stored in models.")
  }
} else {
  cat("✓ All models have the same number of observations:", nrow(data_model1), "\n")
  
  # Verify patients are the same
  if (!is.null(mrns1) && !is.null(mrns2) && !is.null(mrns3) && !is.null(mrns4)) {
    all_identical <- identical(mrns1, mrns2) && identical(mrns1, mrns3) && identical(mrns1, mrns4)
    if (all_identical) {
      cat("✓ All models use the same patients (", length(mrns1), " patients)\n\n")
    } else {
      cat("⚠ WARNING: Models have different patients!\n")
      cat("Model 1 patients:", length(mrns1), "\n")
      cat("Model 2 patients:", length(mrns2), "\n")
      cat("Model 3 patients:", length(mrns3), "\n")
      cat("Model 4 patients:", length(mrns4), "\n")
      stop("Models have different patients. This will cause issues with loo_compare.")
    }
  } else {
    cat("⚠ Could not verify patient identity (MRNs not stored)\n\n")
  }
}

# ============================================================================
# 6. HELPER FUNCTION: TRAIN MODEL
# ============================================================================

train_model <- function(model_data, model_name, file_name) {
  cat("\n=== Training", model_name, "===\n")
  
  # Get predictor names
  predictor_names <- setdiff(names(model_data), "RILD")
  
  # Build formula
  formula_str <- paste("RILD ~", paste(predictor_names, collapse = " + "))
  
  cat("Predictors:", length(predictor_names), "\n")
  cat("Fitting model (this may take 10-30 minutes)...\n")
  
  # Fit model with horseshoe priors for variable selection
  model <- brm(
    formula = as.formula(formula_str),
    data = model_data,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(0, 2.5), class = "Intercept"),
      prior(horseshoe(par_ratio = 0.1), class = "b")
    ),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 42,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    file = file_name,
    file_refit = "on_change"
  )
  
  cat("Model fitting completed!\n")
  return(model)
}

# ============================================================================
# 7. TRAIN ALL 4 MODELS
# ============================================================================

cat("\n=== TRAINING ALL MODELS ===\n\n")

# Model 1: Baseline without radiomics
model1 <- train_model(data_model1, "Model 1: Baseline without radiomics", 
                      "model_baseline_no_radiomics.rds")

# Model 2: Baseline with Primary radiomics
model2 <- train_model(data_model2, "Model 2: Baseline with Primary radiomics", 
                      "model_baseline_with_radiomics.rds")

# Model 3: Mid-treatment without radiomics
model3 <- train_model(data_model3, "Model 3: Mid-treatment without radiomics", 
                      "model_midtreatment_no_radiomics.rds")

# Model 4: Mid-treatment with CTORMid radiomics
model4 <- train_model(data_model4, "Model 4: Mid-treatment with CTORMid radiomics", 
                      "model_midtreatment_with_radiomics.rds")

# ============================================================================
# 8. EVALUATE AND COMPARE MODELS
# ============================================================================

cat("\n=== EVALUATING MODELS ===\n\n")

# Function to calculate performance metrics
calculate_performance <- function(model, model_data, model_name) {
  # Get predictions
  pred_probs <- posterior_linpred(model, transform = TRUE)
  mean_probs <- colMeans(pred_probs)
  
  # Get actual outcomes
  actual <- as.numeric(model_data$RILD) - 1
  
  # Calculate metrics
  pred_class <- ifelse(mean_probs > 0.5, 1, 0)
  accuracy <- mean(actual == pred_class)
  
  # Confusion matrix
  cm <- table(Actual = actual, Predicted = pred_class)
  if (nrow(cm) == 2 && ncol(cm) == 2) {
    sensitivity <- cm[2,2] / sum(cm[2,])
    specificity <- cm[1,1] / sum(cm[1,])
    ppv <- cm[2,2] / sum(cm[,2])
    npv <- cm[1,1] / sum(cm[,1])
  } else {
    sensitivity <- specificity <- ppv <- npv <- NA
  }
  
  # ROC curve
  roc_obj <- roc(actual, mean_probs, quiet = TRUE)
  auc_value <- as.numeric(auc(roc_obj))
  
  # LOO-CV
  loo_result <- tryCatch(loo(model), error = function(e) NULL)
  
  # WAIC (Widely Applicable Information Criterion)
  waic_result <- tryCatch(waic(model), error = function(e) NULL)
  
  return(list(
    model_name = model_name,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    auc = auc_value,
    roc_obj = roc_obj,
    loo = loo_result,
    waic = waic_result,
    predictions = mean_probs,
    actual = actual
  ))
}

# Evaluate all models
results1 <- calculate_performance(model1, data_model1, "Baseline without radiomics")
results2 <- calculate_performance(model2, data_model2, "Baseline with Primary radiomics")
results3 <- calculate_performance(model3, data_model3, "Mid-treatment without radiomics")
results4 <- calculate_performance(model4, data_model4, "Mid-treatment with CTORMid radiomics")

# ============================================================================
# 9. COMPARISON SUMMARY
# ============================================================================

# ============================================================================
# ... [KEEP PREVIOUS PARTS 1-7 OF THE SCRIPT UNCHANGED] ...
# ============================================================================

# ============================================================================
# 8. EXTRACT RESULTS FOR MANUSCRIPT
# ============================================================================

cat("\n\n======================================================================\n")
cat("                  GENERATING MANUSCRIPT RESULTS                       \n")
cat("======================================================================\n\n")

# Helper function to get text-ready metrics
get_manuscript_metrics <- function(model, data, model_name) {
  
  # 1. Get Posterior Predictions
  preds <- posterior_linpred(model, transform = TRUE)
  mean_preds <- colMeans(preds)
  actual <- as.numeric(data$RILD) - 1
  
  # 2. AUC and 95% CI (using DeLong method)
  roc_obj <- roc(actual, mean_preds, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  ci_val <- ci.auc(roc_obj) # 95% CI
  
  # 3. Cohen's d (Approximation from AUC)
  # Formula: d = sqrt(2) * qnorm(AUC)
  cohens_d <- sqrt(2) * qnorm(auc_val)
  
  # 4. Standard Metrics (at 0.5 threshold)
  pred_class <- ifelse(mean_preds > 0.5, 1, 0)
  cm <- table(Factor = factor(actual, levels=0:1), 
              Pred = factor(pred_class, levels=0:1))
  
  sens <- cm[2,2] / sum(cm[2,])  # TPR
  spec <- cm[1,1] / sum(cm[1,])  # TNR
  ppv  <- cm[2,2] / sum(cm[,2])  # Precision
  
  # 5. WAIC
  w <- waic(model)
  waic_val <- w$estimates["waic", "Estimate"]
  
  return(list(
    name = model_name,
    auc = auc_val,
    auc_ci_lower = ci_val[1],
    auc_ci_upper = ci_val[3],
    sd_proxy = (ci_val[3] - ci_val[1]) / 3.92, # Approx SD from 95% CI
    cohens_d = cohens_d,
    sens = sens,
    spec = spec,
    ppv = ppv,
    waic = waic_val
  ))
}

# Calculate for relevant models
m1_res <- get_manuscript_metrics(model1, data_model1, "PreTx (Baseline)")
m4_res <- get_manuscript_metrics(model4, data_model4, "MidTx (+ Radiomics)")

# ============================================================================
# 9. OUTPUT DATA FOR TABLES
# ============================================================================

cat("--- DATA FOR TABLE 3 (Predictive Performance) ---\n")
cat("Copy these numbers into your LaTeX Table 3:\n\n")

cat(sprintf("Model: %s\n", m1_res$name))
cat(sprintf("  AUC: %.2f (95%% CI: %.2f - %.2f)\n", m1_res$auc, m1_res$auc_ci_lower, m1_res$auc_ci_upper))
cat(sprintf("  Cohen's d: %.2f\n", m1_res$cohens_d))
cat("\n")
cat(sprintf("Model: %s\n", m4_res$name))
cat(sprintf("  AUC: %.2f (95%% CI: %.2f - %.2f)\n", m4_res$auc, m4_res$auc_ci_lower, m4_res$auc_ci_upper))
cat(sprintf("  Cohen's d: %.2f\n", m4_res$cohens_d))

cat("\n\n--- DATA FOR TABLE 4 (Bayesian Metrics) ---\n")
cat("Copy these numbers into your LaTeX Table 4:\n\n")
cat(sprintf("%-20s | %-5s | %-5s | %-5s | %-5s | %-5s\n", "Model", "AUC", "Sens", "Spec", "PPV", "WAIC"))
cat("-------------------------------------------------------------\n")
cat(sprintf("%-20s | %.3f | %.3f | %.3f | %.3f | %.1f\n", 
            "Baseline (No Rad)", m1_res$auc, m1_res$sens, m1_res$spec, m1_res$ppv, m1_res$waic))
cat(sprintf("%-20s | %.3f | %.3f | %.3f | %.3f | %.1f\n", 
            "MidTx (+ Rad)", m4_res$auc, m4_res$sens, m4_res$spec, m4_res$ppv, m4_res$waic))

# ============================================================================
# 10. FEATURE IMPORTANCE (PROBABILITY OF DIRECTION)
# ============================================================================

cat("\n\n--- DATA FOR RESULTS TEXT (Feature Importance) ---\n")
cat("Identifying features with high Probability of Direction (pd)...\n")
cat("Using Model 4 (Mid-treatment with Radiomics)\n\n")

# Extract posterior samples
draws <- as_draws_matrix(model4)

# Identify coefficient columns (usually start with 'b_')
# Exclude Intercept
beta_cols <- grep("^b_", colnames(draws), value = TRUE)
beta_cols <- beta_cols[beta_cols != "b_Intercept"]

# Calculate statistics for each feature
feature_stats <- data.frame(
  Feature = beta_cols,
  Mean = apply(draws[, beta_cols], 2, mean),
  # Probability of Direction: max probability of being >0 or <0
  pd = apply(draws[, beta_cols], 2, function(x) {
    max(mean(x > 0), mean(x < 0))
  })
) %>%
  arrange(desc(pd)) %>%
  mutate(
    Direction = ifelse(Mean > 0, "Positive", "Negative"),
    Feature_Clean = gsub("b_", "", Feature) # Clean up name
  )

# Get Top 5 predictors
top_features <- head(feature_stats, 5)

cat("Top 5 Predictors of RILD (Highest Probability of Direction):\n")
print(top_features[, c("Feature_Clean", "Direction", "pd")])

cat("\n\n--- SENTENCE FILLER ---\n")
best_feat <- top_features[1,]
worst_feat <- tail(feature_stats, 1)

cat(sprintf("For your text:\n"))
cat(sprintf("'The Bayesian analysis revealed that %s had a high probability of being %sly associated with RILD (Probability of Direction = %.2f%%), whereas %s had a negligible effect.'\n", 
            best_feat$Feature_Clean, 
            tolower(best_feat$Direction), 
            best_feat$pd * 100, 
            worst_feat$Feature_Clean))