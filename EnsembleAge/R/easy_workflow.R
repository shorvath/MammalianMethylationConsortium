#' EnsembleAge Easy Workflow Functions
#'
#' These functions provide a simple 4-step workflow for users:
#' 1. load_methylation_data() - Load methylation data
#' 2. load_sample_sheet() - Load and prepare sample information  
#' 3. preprocess_data() - Automatically detect platform and preprocess
#' 4. predict_ages() - Run clock predictions
#'
#' @name easy_workflow
NULL

#' Load Methylation Data
#'
#' Load methylation data from various formats (RDS, CSV, etc.)
#' Automatically handles different orientations and formats.
#'
#' @param file_path Character string path to the methylation data file
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Data frame or matrix containing methylation data
#' @export
#' @examples
#' \dontrun{
#' # Load any methylation data format
#' data <- load_methylation_data("path/to/methylation_data.RDS")
#' }
load_methylation_data <- function(file_path, verbose = TRUE) {
  
  if (verbose) cat("Loading methylation data from:", basename(file_path), "\n")
  
  # Check file extension and load accordingly
  file_ext <- tools::file_ext(file_path)
  
  if (file_ext == "RDS" || file_ext == "rds") {
    data <- readRDS(file_path)
  } else if (file_ext == "RData" || file_ext == "rdata") {
    # Load .RData file and find the methylation data object
    env <- new.env()
    load(file_path, envir = env)
    loaded_objects <- ls(env)
    
    if (verbose) cat("  Objects in .RData:", paste(loaded_objects, collapse = ", "), "\n")
    
    # Look for common methylation data object names
    potential_names <- c('betas', 'dat0Noob', 'dat0', 'methylation_data', 'beta_values', 'data')
    data <- NULL
    
    for (name in potential_names) {
      if (name %in% loaded_objects) {
        obj <- get(name, envir = env)
        if ((is.matrix(obj) || is.data.frame(obj)) && nrow(obj) > 100 && ncol(obj) > 5) {
          data <- obj
          if (verbose) cat("  Using object:", name, "\n")
          break
        }
      }
    }
    
    # If no common names found, use the largest object that looks like data
    if (is.null(data)) {
      for (obj_name in loaded_objects) {
        obj <- get(obj_name, envir = env)
        if ((is.matrix(obj) || is.data.frame(obj)) && nrow(obj) > 100 && ncol(obj) > 5) {
          data <- obj
          if (verbose) cat("  Using object:", obj_name, "(auto-detected)\n")
          break
        }
      }
    }
    
    if (is.null(data)) {
      stop("Could not find methylation data in .RData file. Expected a matrix or data.frame with >100 rows and >5 columns.")
    }
    
    # Handle missing CpG IDs for human data
    if ((nrow(data) > 100000) && is.null(rownames(data)) || all(grepl("^[0-9]+$", rownames(data)[1:10]))) {
      if (verbose) cat("  Detected human data missing CpG IDs. Adding rownames as probe indices.\n")
      # For now, add generic probe names - this will need annotation mapping later
      rownames(data) <- paste0("probe_", seq_len(nrow(data)))
    }
    
  } else if (file_ext == "csv") {
    data <- read.csv(file_path, row.names = 1)
  } else if (file_ext == "txt" || file_ext == "tsv") {
    data <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1)
  } else {
    stop("Unsupported file format: ", file_ext, 
         ". Supported formats: RDS, RData, CSV, TXT, TSV")
  }
  
  if (verbose) {
    cat("✓ Data loaded successfully\n")
    cat("  Dimensions:", dim(data), "\n")
    cat("  Data type:", class(data), "\n")
  }
  
  return(data)
}

#' Load and Prepare Sample Sheet
#'
#' Load sample sheet and prepare it with required columns for age prediction.
#' Automatically fills in missing required fields with defaults.
#'
#' @param file_path Character string path to the sample sheet CSV file
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Data frame with prepared sample information
#' @export
#' @examples
#' \dontrun{
#' # Load sample sheet
#' samples <- load_sample_sheet("path/to/sample_sheet.csv")
#' }
load_sample_sheet <- function(file_path, verbose = TRUE) {
  
  if (verbose) cat("Loading sample sheet from:", basename(file_path), "\n")
  
  # Load the sample sheet
  samples <- read.csv(file_path)
  
  if (verbose) {
    cat("✓ Sample sheet loaded\n")
    cat("  Samples:", nrow(samples), "\n")
    cat("  Columns:", ncol(samples), "\n")
  }
  
  # Prepare the sample sheet using existing function
  samples_prepared <- prepare_sample_sheet(samples)
  
  if (verbose) {
    cat("✓ Sample sheet prepared with required columns\n")
    missing_vars <- c()
    if (all(is.na(samples_prepared$Female))) missing_vars <- c(missing_vars, "Sex/Female")
    if (all(is.na(samples_prepared$Age))) missing_vars <- c(missing_vars, "Age")
    if (length(missing_vars) > 0) {
      cat("  Note: Missing variables:", paste(missing_vars, collapse = ", "), "\n")
      cat("  (Defaults applied - provide these for better predictions)\n")
    }
  }
  
  return(samples_prepared)
}

#' Preprocess Methylation Data
#'
#' Automatically detect platform and preprocess methylation data.
#' Handles all platforms: Mammal40k, Mammal320k, HumanEPIC, Human450k.
#' 
#' @param methylation_data Data frame or matrix from load_methylation_data()
#' @param sample_sheet Data frame from load_sample_sheet()
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Preprocessed data frame ready for age prediction
#' @export
#' @examples
#' \dontrun{
#' # Preprocess data (automatic platform detection)
#' processed_data <- preprocess_data(methylation_data, sample_sheet)
#' }
preprocess_data <- function(methylation_data, sample_sheet, verbose = TRUE) {
  
  if (verbose) cat("\n=== Preprocessing Methylation Data ===\n")
  
  # Use the existing preprocessing pipeline
  processed_result <- preprocess_methylation_data(
    dat0sesame = methylation_data,
    samps = sample_sheet,
    verbose = verbose
  )
  
  # Extract the data component from the list
  processed_data <- processed_result$data
  
  if (verbose) {
    cat("✓ Preprocessing complete!\n")
    if ("CGid" %in% names(processed_data)) {
      cat("  Final dimensions:", nrow(processed_data), "CpG sites ×", ncol(processed_data) - 1, "samples\n")
      
      # Check coverage for some popular clocks
      tryCatch({
        coverage <- check_probe_coverage(processed_data)
        avg_coverage <- round(mean(coverage$Coverage_Percent), 1)
        high_coverage <- sum(coverage$Coverage_Percent >= 90)
        cat("  Average clock coverage:", avg_coverage, "%\n")
        cat("  High coverage clocks (≥90%):", high_coverage, "out of", nrow(coverage), "\n")
      }, error = function(e) {
        cat("  Coverage check skipped\n")
      })
    } else {
      cat("  Final dimensions:", dim(processed_data), "\n")
    }
  }
  
  return(processed_data)
}

#' Predict Ages Using Multiple Clocks
#'
#' Run age predictions using ensemble and individual clocks.
#' Automatically selects the best prediction method based on data coverage.
#'
#' @param processed_data Preprocessed data from preprocess_data()
#' @param sample_sheet Sample sheet from load_sample_sheet()
#' @param method Character string specifying prediction method:
#'   - "auto" (default): Automatically select best method
#'   - "ensemble_static": Static ensemble clocks
#'   - "ensemble_dynamic": Dynamic ensemble clocks  
#'   - "ensemble_dual": Dual static ensemble clocks
#'   - "all": All available clocks
#' @param save_results Logical indicating whether to save results to CSV (default: TRUE)
#' @param output_file Character string specifying output file name (default: auto-generated)
#' @param output_format Character string specifying output format:
#'   - "long" (default): One row per sample-clock combination
#'   - "wide": One row per sample, clocks as columns  
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Data frame with age predictions
#' @export
#' @examples
#' \dontrun{
#' # Predict ages (automatic method selection)
#' results <- predict_ages(processed_data, sample_sheet)
#' 
#' # Use specific method
#' results <- predict_ages(processed_data, sample_sheet, method = "ensemble_static")
#' }
predict_ages <- function(processed_data, sample_sheet, method = "auto", 
                        save_results = TRUE, output_file = NULL, output_format = "long", 
                        verbose = TRUE) {
  
  if (verbose) cat("\n=== Predicting Ages ===\n")
  
  # Determine best method if auto
  if (method == "auto") {
    # Check data coverage to select best method
    tryCatch({
      coverage <- check_probe_coverage(processed_data)
      avg_coverage <- mean(coverage$Coverage_Percent)
      
      if (avg_coverage >= 95) {
        method <- "ensemble_static"
        if (verbose) cat("High coverage detected - using ensemble static clocks\n")
      } else if (avg_coverage >= 85) {
        method <- "ensemble_dynamic" 
        if (verbose) cat("Good coverage detected - using ensemble dynamic clocks\n")
      } else if (avg_coverage >= 70) {
        method <- "ensemble_dual"
        if (verbose) cat("Moderate coverage detected - using dual ensemble clocks\n")
      } else {
        method <- "all"
        if (verbose) cat("Lower coverage detected - using all available clocks\n")
      }
    }, error = function(e) {
      method <- "ensemble_static"
      if (verbose) cat("Coverage check failed - defaulting to ensemble static\n")
    })
  }
  
  if (verbose) cat("Using prediction method:", method, "\n")
  
  # Run predictions based on method
  tryCatch({
    results <- switch(method,
      "ensemble_static" = predict_ensemble_static(processed_data, sample_sheet),
      "ensemble_dynamic" = predict_ensemble_dynamic(processed_data, sample_sheet),
      "ensemble_dual" = predict_ensemble_dual_static(processed_data, sample_sheet),
      "all" = predict_all_clocks(processed_data, sample_sheet),
      stop("Unknown method: ", method)
    )
    
    if (verbose) {
      cat("✓ Age predictions successful!\n")
      cat("  Results dimensions:", dim(results), "\n")
      
      # Count age prediction columns
      age_cols <- grep("epiAge", colnames(results), value = TRUE)
      cat("  Age predictions:", length(age_cols), "clocks\n")
      
      # Show sample results
      if (nrow(results) >= 1 && length(age_cols) >= 1) {
        cat("  Sample prediction (first clock):", 
            round(results[1, age_cols[1]], 3), "years\n")
      }
    }
    
    # Convert to requested format
    final_results <- convert_results_format(results, output_format, verbose)
    
    # Save results if requested
    if (save_results) {
      if (is.null(output_file)) {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        output_file <- paste0("EnsembleAge_results_", method, "_", output_format, "_", timestamp, ".csv")
      }
      
      write.csv(final_results, output_file, row.names = FALSE)
      
      if (verbose) {
        file_size <- round(file.info(output_file)$size / 1024, 1)
        cat("✓ Results saved to:", output_file, "(", file_size, "KB)\n")
      }
    }
    
    return(final_results)
    
  }, error = function(e) {
    if (verbose) cat("✗ Age prediction failed:", e$message, "\n")
    
    # Try fallback method
    if (method != "all") {
      if (verbose) cat("Trying fallback method: all clocks...\n")
      return(predict_ages(processed_data, sample_sheet, method = "all", 
                         save_results = save_results, output_file = output_file, 
                         verbose = verbose))
    } else {
      stop("All prediction methods failed: ", e$message)
    }
  })
}

#' Complete EnsembleAge Workflow
#'
#' Run the complete 4-step workflow in one function call.
#' This is the easiest way to use EnsembleAge.
#'
#' @param methylation_file Character string path to methylation data file
#' @param sample_sheet_file Character string path to sample sheet CSV file  
#' @param method Character string specifying prediction method (default: "auto")
#' @param save_results Logical indicating whether to save results (default: TRUE)
#' @param output_file Character string for output file name (default: auto-generated)
#' @param output_format Character string specifying output format: "long" or "wide" (default: "long")
#' @param verbose Logical indicating whether to print progress (default: TRUE)
#' @return Data frame with age predictions
#' @export
#' @examples
#' \dontrun{
#' # Complete workflow in one line!
#' results <- run_ensemble_age(
#'   methylation_file = "path/to/methylation_data.RDS",
#'   sample_sheet_file = "path/to/sample_sheet.csv"
#' )
#' 
#' # With specific method
#' results <- run_ensemble_age(
#'   methylation_file = "methylation_data.RDS",
#'   sample_sheet_file = "sample_sheet.csv", 
#'   method = "ensemble_static"
#' )
#' }
run_ensemble_age <- function(methylation_file, sample_sheet_file, method = "auto",
                            save_results = TRUE, output_file = NULL, output_format = "long", 
                            verbose = TRUE) {
  
  if (verbose) {
    cat("🧬 EnsembleAge: Complete Methylation Age Prediction Workflow 🧬\n")
    cat("================================================================\n\n")
  }
  
  # Step 1: Load methylation data
  if (verbose) cat("📁 Step 1: Loading methylation data...\n")
  methylation_data <- load_methylation_data(methylation_file, verbose = verbose)
  
  # Step 2: Load sample sheet
  if (verbose) cat("\n📋 Step 2: Loading sample sheet...\n")
  sample_sheet <- load_sample_sheet(sample_sheet_file, verbose = verbose)
  
  # Step 3: Preprocess data
  if (verbose) cat("\n🔄 Step 3: Preprocessing data...\n")
  processed_data <- preprocess_data(methylation_data, sample_sheet, verbose = verbose)
  
  # Step 4: Predict ages
  if (verbose) cat("\n🎯 Step 4: Predicting ages...\n")
  results <- predict_ages(processed_data, sample_sheet, method = method,
                         save_results = save_results, output_file = output_file,
                         output_format = output_format, verbose = verbose)
  
  if (verbose) {
    cat("\n🎉 EnsembleAge Workflow Complete! 🎉\n")
    cat("================================\n")
    cat("✅ Data loaded and preprocessed\n")
    cat("✅ Platform automatically detected\n") 
    cat("✅ Age predictions generated\n")
    if (save_results) cat("✅ Results saved to file\n")
    cat("\nReady for analysis! 📊\n")
  }
  
  return(results)
}

#' Convert Results Format
#'
#' Convert prediction results between wide and long formats.
#' 
#' @param results Data frame with prediction results 
#' @param target_format Character string: "wide" or "long"
#' @param verbose Logical indicating whether to print messages
#' @return Data frame in the requested format
#' @keywords internal
convert_results_format <- function(results, target_format, verbose = FALSE) {
  
  # Detect current format
  has_epiclock_col <- "epiClock" %in% names(results)
  has_age_cols <- any(grepl("\\.epiAge$|\\.clock\\.epiAge$", names(results)))
  
  if (has_epiclock_col && !has_age_cols) {
    current_format <- "long"
  } else if (!has_epiclock_col && has_age_cols) {
    current_format <- "wide"
  } else {
    # Try to detect from long format attribute
    long_format_data <- attr(results, "long_format")
    if (!is.null(long_format_data)) {
      current_format <- "wide"
      if (verbose) cat("Using wide format with long format attribute\n")
    } else {
      current_format <- "long"  # Default assumption
    }
  }
  
  if (verbose) cat("Converting from", current_format, "to", target_format, "format\n")
  
  # Return if already in target format
  if (current_format == target_format) {
    return(results)
  }
  
  # Convert between formats
  if (current_format == "wide" && target_format == "long") {
    # Wide to long conversion
    long_format_data <- attr(results, "long_format")
    if (!is.null(long_format_data)) {
      return(long_format_data)
    } else {
      # Manual conversion if no long format attribute
      stop("Cannot convert to long format: no long format data available")
    }
    
  } else if (current_format == "long" && target_format == "wide") {
    # Long to wide conversion
    if (!"epiClock" %in% names(results)) {
      stop("Cannot convert to wide format: missing epiClock column")
    }
    
    # Create exact same naming pattern as core predictAgeAndAgeAcc function
    # Pattern: clockFamily.epiClock.clock.epiAge and clockFamily.epiClock.clock.AgeAcceleration
    results_clean <- results %>%
      dplyr::group_by(Basename, clockFamily, epiClock) %>%
      dplyr::slice(1) %>%  # Take first occurrence if duplicates
      dplyr::ungroup()
    
    # Convert age predictions to wide format (matching core function pattern)
    wide_ages <- results_clean %>%
      dplyr::mutate(clock_age_name = paste(clockFamily, epiClock, "clock", "epiAge", sep = ".")) %>%
      dplyr::select(Basename, clock_age_name, epiAge) %>%
      tidyr::spread(key = clock_age_name, value = epiAge)
    
    # Convert age acceleration to wide format (matching core function pattern)
    wide_accel <- results_clean %>%
      dplyr::mutate(clock_accel_name = paste(clockFamily, epiClock, "clock", "AgeAcceleration", sep = ".")) %>%
      dplyr::select(Basename, clock_accel_name, AgeAccelation) %>%
      tidyr::spread(key = clock_accel_name, value = AgeAccelation)
    
    # Combine with sample info
    sample_info <- results %>%
      dplyr::select(Basename, Age, Female, Tissue) %>%
      dplyr::distinct()
    
    # Merge all together
    wide_results <- sample_info %>%
      dplyr::left_join(wide_ages, by = "Basename") %>%
      dplyr::left_join(wide_accel, by = "Basename")
    
    return(wide_results)
  }
  
  return(results)
}
