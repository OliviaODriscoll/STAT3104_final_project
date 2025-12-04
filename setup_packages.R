# Setup script to install all required packages
# Run this script first: source("setup_packages.R")

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

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

cat("=== Installing Required R Packages ===\n\n")

# Function to install packages if missing
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat("✓", pkg, "is already installed\n")
  }
}

# Special handling for brms and rstan
cat("\n=== Checking Stan-related packages ===\n")
if (!require("rstan", quietly = TRUE)) {
  cat("Installing rstan (this may take 10-15 minutes)...\n")
  install.packages("rstan", dependencies = TRUE)
} else {
  cat("✓ rstan is installed\n")
}

if (!require("brms", quietly = TRUE)) {
  cat("Installing brms (this may take a few minutes)...\n")
  install.packages("brms", dependencies = TRUE)
} else {
  cat("✓ brms is installed\n")
}

# Final verification
cat("\n=== Final Package Check ===\n")
all_installed <- TRUE
for (pkg in required_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "\n")
  } else {
    cat("✗", pkg, "FAILED TO INSTALL\n")
    all_installed <- FALSE
  }
}

if (all_installed) {
  cat("\n✓ All packages installed successfully!\n")
  cat("You can now run bayesian_rild_analysis.R\n")
} else {
  cat("\n✗ Some packages failed to install.\n")
  cat("Please install them manually:\n")
  cat("install.packages(c('", paste(required_packages, collapse = "', '"), "'))\n", sep = "")
}

