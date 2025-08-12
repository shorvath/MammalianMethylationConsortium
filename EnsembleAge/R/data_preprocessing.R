#' Preprocess methylation data for EnsembleAge prediction
#'
#' This function preprocesses methylation data to ensure it's in the correct format
#' for age prediction. It handles missing values, filters probes, and validates the data.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @param samps Data frame containing sample information
#' @param min_coverage Minimum probe coverage required (default: 0.8)
#' @param handle_missing Character string indicating how to handle missing values:
#'   "remove" (remove samples/probes with too many NAs), "impute" (simple mean imputation),
#'   or "keep" (keep as is). Default: "impute"
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return List containing processed data and metadata
#' @export
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @examples
#' \dontrun{
#' processed <- preprocess_methylation_data(methylation_data, sample_info)
#' prediction_results <- predictAgeAndAgeAcc(processed$data, processed$samples)
#' }
preprocess_methylation_data <- function(dat0sesame, samps, min_coverage = 0.8, 
                                      handle_missing = "impute", verbose = TRUE) {
  
  if (verbose) cat("Starting data preprocessing...\n")
  
  # First, detect and fix data orientation
  if (verbose) cat("Checking data orientation...\n")
  orientation_result <- detect_and_fix_orientation(dat0sesame, samps, verbose = verbose)
  dat0sesame <- orientation_result$data
  
  # Check if this is Mammal320k data that needs complete processing
  if (verbose) cat("Checking for platform-specific processing requirements...\n")
  
  # First detect platform - handle matrices and data frames differently
  if (is.matrix(dat0sesame)) {
    # For matrices, check dimensions and patterns
    n_rows <- nrow(dat0sesame)
    n_cols <- ncol(dat0sesame)
    
    if (n_cols > 250000) {
      col_suffix_pattern <- sum(grepl("_[A-Z]+[0-9]*$", colnames(dat0sesame)[1:min(1000, n_cols)]))
      if (col_suffix_pattern > 100) {
        platform_check <- "Mammal320k"
      } else {
        platform_check <- "Unknown"
      }
    } else if (n_rows > 250000) {
      row_suffix_pattern <- sum(grepl("_[A-Z]+[0-9]*$", rownames(dat0sesame)[1:min(1000, n_rows)]))
      if (row_suffix_pattern > 100) {
        platform_check <- "Mammal320k"
      } else {
        platform_check <- "Unknown"
      }
    } else {
      platform_check <- "Unknown"
    }
  } else {
    # For data frames, use standard detection
    tryCatch({
      platform_check <- detect_platform(dat0sesame)
    }, error = function(e) {
      # If standard detection fails, check for mammal320k patterns
      if ("CGid" %in% names(dat0sesame)) {
        suffix_pattern <- sum(grepl("_[A-Z]+[0-9]*$", dat0sesame$CGid[1:min(100, nrow(dat0sesame))]))
        platform_check <- if (suffix_pattern > 10) "Mammal320k" else "Unknown"
      } else {
        platform_check <- "Unknown"
      }
    })
  }
  
  if (platform_check == "Mammal320k") {
    if (verbose) cat("Detected Mammal320k format. Using complete processing workflow...\n")
    
    # Use the complete mammal320k processing workflow following Amin's exact steps
    # Default to mouse species unless specified otherwise
    species <- "mouse"  # Could be made a parameter in future
    dat0sesame <- process_mammal320k_data(dat0sesame, samps, species = species, verbose = verbose)
  }
  
  # Validate input data after orientation fix
  validation <- validate_data_format(dat0sesame, samps)
  if (!validation$valid) {
    stop("Data validation failed:\n", paste(validation$issues, collapse = "\n"))
  }
  
  if (verbose) {
    cat("Platform detected:", validation$platform, "\n")
    cat("Number of probes:", validation$n_probes, "\n")
    cat("Number of samples:", validation$n_samples, "\n")
  }
  
  # Prepare sample sheet
  samps_processed <- prepare_sample_sheet(samps)
  
  # Filter samples that exist in both datasets
  available_samples <- intersect(samps_processed$Basename, names(dat0sesame))
  if (length(available_samples) == 0) {
    stop("No matching samples found between dat0sesame and samps")
  }
  
  if (length(available_samples) < nrow(samps_processed)) {
    missing_samples <- setdiff(samps_processed$Basename, available_samples)
    if (verbose) {
      cat("Warning: Some samples in samps not found in data:", 
          paste(missing_samples[1:min(3, length(missing_samples))], collapse = ", "), "\n")
    }
  }
  
  # Filter data to matching samples
  dat_filtered <- dat0sesame %>% select(CGid, all_of(available_samples))
  samps_filtered <- samps_processed[samps_processed$Basename %in% available_samples, ]
  
  # Check probe coverage
  coverage <- check_probe_coverage(dat_filtered)
  low_coverage_clocks <- coverage[coverage$Coverage_Percent < min_coverage * 100, ]
  
  if (nrow(low_coverage_clocks) > 0 && verbose) {
    cat("Warning: Some clocks have low probe coverage (<", min_coverage * 100, "%):\n")
    print(low_coverage_clocks[1:min(5, nrow(low_coverage_clocks)), ])
  }
  
  # Handle missing values
  if (handle_missing == "remove") {
    # Remove probes with too many missing values
    probe_na_prop <- rowMeans(is.na(dat_filtered[, -1]))
    probes_to_keep <- probe_na_prop < (1 - min_coverage)
    dat_filtered <- dat_filtered[probes_to_keep, ]
    
    # Remove samples with too many missing values
    sample_na_prop <- colMeans(is.na(dat_filtered[, -1]))
    samples_to_keep <- c(TRUE, sample_na_prop < (1 - min_coverage))
    dat_filtered <- dat_filtered[, samples_to_keep]
    
    # Update sample sheet
    kept_samples <- names(dat_filtered)[-1]
    samps_filtered <- samps_filtered[samps_filtered$Basename %in% kept_samples, ]
    
    if (verbose) {
      cat("Removed", sum(!probes_to_keep), "probes and", 
          sum(!samples_to_keep[-1]), "samples due to missing values\n")
    }
    
  } else if (handle_missing == "impute") {
    # Simple mean imputation for missing values
    if (ncol(dat_filtered) > 1) {
      for (i in 2:ncol(dat_filtered)) {
        missing_idx <- is.na(dat_filtered[, i])
        if (any(missing_idx)) {
          dat_filtered[missing_idx, i] <- mean(dat_filtered[, i], na.rm = TRUE)
        }
      }
      if (verbose) cat("Imputed missing values with column means\n")
    } else {
      if (verbose) cat("No sample columns found for imputation\n")
    }
    
  } else if (handle_missing == "keep") {
    if (verbose) cat("Keeping missing values as is\n")
  }
  
  # Final validation
  final_validation <- validate_data_format(dat_filtered, samps_filtered)
  if (!final_validation$valid) {
    warning("Processed data still has issues:\n", paste(final_validation$issues, collapse = "\n"))
  }
  
  if (verbose) {
    cat("Preprocessing complete!\n")
    cat("Final data dimensions:", nrow(dat_filtered), "probes x", ncol(dat_filtered) - 1, "samples\n")
  }
  
  return(list(
    data = dat_filtered,
    samples = samps_filtered,
    platform = validation$platform,
    coverage = coverage,
    removed_samples = setdiff(samps$Basename, samps_filtered$Basename),
    metadata = list(
      original_probes = validation$n_probes,
      original_samples = validation$n_samples,
      final_probes = nrow(dat_filtered),
      final_samples = ncol(dat_filtered) - 1,
      preprocessing_method = handle_missing,
      min_coverage = min_coverage
    )
  ))
}

#' Quick prediction wrapper with automatic preprocessing
#'
#' This function provides a simple wrapper that automatically preprocesses the data
#' and runs age prediction with sensible defaults.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @param samps Data frame containing sample information
#' @param efficient_loading Logical indicating whether to use efficient loading (default: TRUE)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @param ... Additional arguments passed to preprocessing functions
#' @return Data frame with age predictions and acceleration values
#' @export
#' @examples
#' \dontrun{
#' # Simple usage with efficient loading
#' results <- predict_age_simple(methylation_data, sample_info)
#' 
#' # With custom options
#' results <- predict_age_simple(methylation_data, sample_info, 
#'                              efficient_loading = TRUE, verbose = TRUE)
#' }
predict_age_simple <- function(dat0sesame, samps, efficient_loading = TRUE, verbose = TRUE, ...) {
  
  if (efficient_loading) {
    if (verbose) cat("Using efficient data loading for clock predictions...\n")
    processed <- efficiently_load_clock_data(dat0sesame, samps, verbose = verbose, ...)
    dat_clean <- processed$data
    samps_clean <- prepare_sample_sheet(samps)
    
    if (verbose) {
      cat("Platform:", processed$platform, "\n")
      cat("Compression ratio:", processed$compression_ratio, "(", nrow(dat_clean), "probes retained)\n")
      cat("Using", nrow(dat_clean), "probes and", nrow(samps_clean), "samples\n")
    }
  } else {
    if (verbose) cat("Using standard preprocessing...\n")
    processed <- preprocess_methylation_data(dat0sesame, samps, verbose = verbose, ...)
    dat_clean <- processed$data
    samps_clean <- processed$samples
    
    if (verbose) {
      cat("Platform:", processed$platform, "\n")
      cat("Using", nrow(dat_clean), "probes and", nrow(samps_clean), "samples\n")
    }
  }
  
  if (verbose) cat("Running age predictions...\n")
  
  # Run the main prediction function
  results <- predictAgeAndAgeAcc(dat_clean, samps_clean)
  
  if (verbose) cat("Prediction complete!\n")
  
  return(results)
}

#' Convert M-values to beta values
#'
#' This function converts M-values (log2 ratio of methylated to unmethylated intensities)
#' to beta values (proportion of methylated intensities).
#' 
#' @param mvals Numeric vector or matrix of M-values
#' @return Numeric vector or matrix of beta values
#' @export
#' @examples
#' mvals <- c(-2, -1, 0, 1, 2)
#' bvals <- mvals_to_bvals(mvals)
mvals_to_bvals <- function(mvals) {
  2^mvals / (2^mvals + 1)
}

#' Convert beta values to M-values
#'
#' This function converts beta values to M-values.
#' 
#' @param bvals Numeric vector or matrix of beta values (between 0 and 1)
#' @return Numeric vector or matrix of M-values
#' @export
#' @examples
#' bvals <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' mvals <- bvals_to_mvals(bvals)
bvals_to_mvals <- function(bvals) {
  # Add small epsilon to avoid log(0)
  bvals[bvals <= 0] <- 1e-6
  bvals[bvals >= 1] <- 1 - 1e-6
  log2(bvals / (1 - bvals))
}
