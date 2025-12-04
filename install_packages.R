# Script to install all required R packages for the Bayesian RILD analysis

# List of required packages
required_packages <- c(
  "tidyverse",    # Data manipulation and visualization
  "brms",         # Bayesian regression models using Stan
  "bayesplot",    # MCMC diagnostics and posterior plots
  "mice",         # Multiple imputation
  "VIM",          # Visualization of missing data
  "corrplot",     # Correlation plots
  "pROC",         # ROC curves
  "loo"           # Leave-one-out cross-validation
)

# Function to check and install packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages, dependencies = TRUE)
  } else {
    cat("All required packages are already installed.\n")
  }
}

# Install missing packages
install_if_missing(required_packages)

# Special handling for brms (requires RStan)
if (!requireNamespace("brms", quietly = TRUE)) {
  cat("\nInstalling brms (this may take a while as it includes Stan)...\n")
  install.packages("brms", dependencies = TRUE)
}

# Check if RStan is available
if (!requireNamespace("rstan", quietly = TRUE)) {
  cat("\nRStan is required for brms. Installing RStan...\n")
  cat("Note: RStan installation may require additional system dependencies.\n")
  cat("See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details.\n")
  install.packages("rstan", dependencies = TRUE)
}

# Verify installations
cat("\n=== Package Installation Status ===\n")
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "is installed\n")
  } else {
    cat("✗", pkg, "is NOT installed\n")
  }
}

cat("\nInstallation check complete!\n")
cat("If any packages failed to install, please install them manually.\n")

