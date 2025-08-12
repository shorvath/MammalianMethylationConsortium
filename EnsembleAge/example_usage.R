#!/usr/bin/env Rscript
#' EnsembleAge: Simple Usage Example
#' 
#' This script demonstrates how to use EnsembleAge with simulated data.
#' Run this script to test the package installation and basic functionality.

# Load required libraries
library(EnsembleAge)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)

# Source the functions directly (needed for development)
source('R/ultimate_predict.R')
source('R/easy_workflow.R')
source('R/data_preprocessing.R')
source('R/predict_age.R') 
source('R/user_prediction_functions.R')
source('R/core_functions.R')
source('R/efficient_data_loading.R')
source('R/platform_utils.R')
source('R/data_orientation.R')

cat("🧬 EnsembleAge Simple Example\n")
cat("============================\n\n")

# Step 1: Create simulated methylation data
cat("📊 Creating simulated methylation data...\n")
set.seed(123)

# Simulate mammal40k-style data with realistic CpG sites
n_probes <- 1000   # Using subset for quick example
n_samples <- 8

# Create methylation data: probes as rows, samples as columns
methylation_data <- data.frame(
  CGid = paste0("cg", sprintf("%08d", sample(10000000:99999999, n_probes))),
  matrix(
    # Realistic beta values with some biological variation
    pmax(0, pmin(1, rnorm(n_probes * n_samples, mean = 0.5, sd = 0.2))), 
    nrow = n_probes
  )
)

# Add sample names
sample_names <- paste0("Sample_", sprintf("%02d", 1:n_samples))
names(methylation_data)[-1] <- sample_names

cat("✓ Created methylation data:", dim(methylation_data), "(probes × samples)\n")

# Step 2: Create sample sheet with realistic ages
cat("📋 Creating sample sheet...\n")
sample_sheet <- data.frame(
  Basename = sample_names,
  Age = c(0.5, 1.2, 2.1, 0.8, 1.5, 2.8, 1.1, 0.9),  # Mouse ages in years
  Female = c(1, 0, 1, 1, 0, 1, 0, 1),                # Sex: 1=female, 0=male
  Tissue = rep("Blood", n_samples),                    # All blood samples
  SpeciesLatinName = rep("Mus musculus", n_samples)   # Mouse species
)

cat("✓ Created sample sheet with", nrow(sample_sheet), "samples\n")
print(sample_sheet)

cat("\n🚀 Running EnsembleAge predictions...\n")
cat("======================================\n")

# Step 3: Run predictions using the ultimate one-liner
results <- predictEnsemble(methylation_data, sample_sheet, verbose = TRUE)

cat("\n📈 Results Summary:\n")
cat("===================\n")
cat("Total predictions:", nrow(results), "\n")
cat("Unique samples:", length(unique(results$Basename)), "\n")
cat("Unique clocks:", length(unique(results$epiClock)), "\n")
cat("Clock families:", paste(unique(results$clockFamily), collapse = ", "), "\n")

cat("\n📊 Sample Results (first 10 rows):\n")
print(head(results[, c("Basename", "Age", "epiClock", "epiAge", "AgeAccelation")], 10))

cat("\n📈 Age Predictions by Sample:\n")
sample_summary <- results %>%
  group_by(Basename, Age) %>%
  summarise(
    n_clocks = n(),
    mean_epiAge = round(mean(epiAge, na.rm = TRUE), 2),
    mean_acceleration = round(mean(AgeAccelation, na.rm = TRUE), 2),
    .groups = "drop"
  )
print(sample_summary)

cat("\n💾 Saving results...\n")
write.csv(results, "ensembleage_example_results.csv", row.names = FALSE)
cat("✓ Results saved to: ensembleage_example_results.csv\n")

cat("\n🎉 Example completed successfully!\n")
cat("===================================\n\n")

cat("📚 What happened:\n")
cat("1. ✅ Created simulated mammal methylation data (1000 probes × 8 samples)\n")
cat("2. ✅ Created sample sheet with ages, sex, tissue, and species info\n")
cat("3. ✅ EnsembleAge automatically detected the platform\n")
cat("4. ✅ Selected appropriate clock methods based on data coverage\n")
cat("5. ✅ Generated age predictions using multiple clocks\n")
cat("6. ✅ Saved results in long format (one row per sample-clock combination)\n\n")

cat("🔄 Try different options:\n")
cat("- Long format (default): predictEnsemble(methylation_data, sample_sheet)\n")
cat("- Wide format: predictEnsemble(methylation_data, sample_sheet, output_format = 'wide')\n")
cat("- Both formats: predictEnsemble(methylation_data, sample_sheet, output_format = 'both')\n")
cat("- Custom directory: predictEnsemble(data, samples, output_dir = 'results/')\n")
cat("- Custom prefix: predictEnsemble(data, samples, output_prefix = 'MyStudy')\n")
cat("- Custom filename: predictEnsemble(data, samples, output_filename = 'results.csv')\n")
cat("- No saving: predictEnsemble(data, samples, save_results = FALSE)\n\n")

cat("🚀 Ready to use EnsembleAge with your own data!\n")
