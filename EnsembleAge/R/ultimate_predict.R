#' EnsembleAge: Ultimate One-Line Predictor
#'
#' The simplest possible interface - just provide data and sample sheet!
#' Automatically detects platform, preprocesses data, selects best prediction method,
#' and returns results in your preferred format.
#'
#' @param data Either a file path (string) or data object (data.frame/matrix)
#' @param samples Either a file path (string) or sample sheet (data.frame), or NULL to use template
#' @param output_format Character string: "long" (default), "wide", or "both"
#' @param save_results Logical: save results to CSV files (default: TRUE)
#' @param output_dir Character string: directory to save results (default: current directory)
#' @param output_prefix Character string: prefix for output filenames (default: "EnsembleAge")
#' @param output_filename Character string: custom filename (overrides automatic naming, default: NULL)
#' @param verbose Logical: print progress messages (default: TRUE)
#' @return Data frame with age predictions (or list if output_format="both")
#' @export
#' @examples
#' \dontrun{
#' # Ultimate simplicity - just data and samples!
#' results <- predictEnsemble("methylation_data.RDS", "sample_sheet.csv")
#' 
#' # Control output directory and filename
#' results <- predictEnsemble("data.RDS", "samples.csv", 
#'                            output_dir = "results/", 
#'                            output_prefix = "MyStudy")
#' 
#' # Custom filename (overrides automatic naming)
#' results <- predictEnsemble("data.RDS", "samples.csv", 
#'                            output_filename = "my_results.csv")
#' 
#' # Get both formats with custom location
#' results <- predictEnsemble("data.RDS", "samples.csv", 
#'                            output_format = "both",
#'                            output_dir = "output/", 
#'                            output_prefix = "Study1")
#' }
predictEnsemble <- function(data, samples = NULL, output_format = "long", 
                           save_results = TRUE, output_dir = ".", 
                           output_prefix = "EnsembleAge", output_filename = NULL,
                           verbose = TRUE) {
  
  if (verbose) {
    cat("🧬 EnsembleAge: Universal Methylation Age Predictor 🧬\n")
    cat("=====================================================\n\n")
  }
  
  # Step 1: Handle data input (file path or object)
  if (verbose) cat("📁 Loading methylation data...\n")
  if (is.character(data) && length(data) == 1) {
    # It's a file path
    methylation_data <- load_methylation_data(data, verbose = verbose)
    data_name <- tools::file_path_sans_ext(basename(data))
  } else if (is.data.frame(data) || is.matrix(data)) {
    # It's already a data object
    methylation_data <- data
    data_name <- "user_data"
    if (verbose) {
      cat("✓ Data object loaded\n")
      cat("  Dimensions:", dim(methylation_data), "\n")
      cat("  Data type:", class(methylation_data), "\n")
    }
  } else {
    stop("'data' must be a file path (string) or data object (data.frame/matrix)")
  }
  
  # Step 2: Handle sample sheet input
  if (verbose) cat("\n📋 Preparing sample sheet...\n")
  if (is.null(samples)) {
    # Auto-generate sample sheet template
    if (verbose) cat("No sample sheet provided - generating template automatically...\n")
    
    # Determine sample names from data
    if (is.matrix(methylation_data) || is.data.frame(methylation_data)) {
      # Smart detection: If one dimension is much larger than the other, 
      # the larger dimension likely contains probes (especially if >10,000)
      nrows <- nrow(methylation_data)
      ncols <- ncol(methylation_data)
      
      if (ncols > 10000 && nrows < 1000) {
        # Very likely samples as rows, probes as columns (e.g., 5 samples × 326k probes)
        sample_names <- rownames(methylation_data)
        if (is.null(sample_names)) {
          sample_names <- paste0("Sample_", 1:nrows)
        }
      } else if (nrows > 10000 && ncols < 1000) {
        # Very likely probes as rows, samples as columns (e.g., 326k probes × 5 samples)
        sample_names <- colnames(methylation_data)
        if (is.null(sample_names)) {
          sample_names <- paste0("Sample_", 1:ncols)
        }
        # Remove potential CGid column
        sample_names <- sample_names[sample_names != "CGid"]
      } else if (nrows > ncols) {
        # Default: more rows than columns, assume samples as rows
        sample_names <- rownames(methylation_data)
        if (is.null(sample_names)) {
          sample_names <- paste0("Sample_", 1:nrows)
        }
      } else {
        # Default: more columns than rows, assume samples as columns
        sample_names <- colnames(methylation_data)
        if (is.null(sample_names)) {
          sample_names <- paste0("Sample_", 1:ncols)
        }
        # Remove potential CGid column
        sample_names <- sample_names[sample_names != "CGid"]
      }
    }
    
    # Create template sample sheet
    sample_sheet <- data.frame(
      Basename = sample_names,
      Age = 1.0,  # Default age
      Sex = "Unknown",
      Species = "Mus musculus",
      Tissue = "Unknown",
      stringsAsFactors = FALSE
    )
    
    if (verbose) {
      cat("✓ Auto-generated sample sheet with", nrow(sample_sheet), "samples\n")
      cat("  Note: Using default values (Age=1.0, Species='Mus musculus')\n")
      cat("  Provide actual sample sheet for better predictions\n")
    }
    
    # Save template for user
    template_file <- paste0(data_name, "_sample_sheet_template.csv")
    write.csv(sample_sheet, template_file, row.names = FALSE)
    if (verbose) cat("  Template saved as:", template_file, "\n")
    
  } else if (is.character(samples) && length(samples) == 1) {
    # It's a file path
    sample_sheet <- load_sample_sheet(samples, verbose = verbose)
  } else if (is.data.frame(samples)) {
    # It's already a data frame
    sample_sheet <- prepare_sample_sheet(samples)
    if (verbose) {
      cat("✓ Sample sheet object loaded\n")
      cat("  Samples:", nrow(sample_sheet), "\n")
    }
  } else {
    stop("'samples' must be NULL, a file path (string), or data.frame")
  }
  
  # Step 3: Auto-detect platform and preprocess
  if (verbose) cat("\n🔄 Auto-detecting platform and preprocessing...\n")
  processed_data <- preprocess_data(methylation_data, sample_sheet, verbose = verbose)
  
  # Step 4: Auto-select best prediction method
  if (verbose) cat("\n🎯 Auto-selecting optimal prediction method...\n")
  
  # Determine best methods based on coverage
  tryCatch({
    coverage <- check_probe_coverage(processed_data)
    avg_coverage <- mean(coverage$Coverage_Percent)
    high_cov_clocks <- sum(coverage$Coverage_Percent >= 95)
    
    if (verbose) {
      cat("  Platform coverage analysis:\n")
      cat("    Average coverage:", round(avg_coverage, 1), "%\n")
      cat("    High coverage clocks (≥95%):", high_cov_clocks, "\n")
    }
    
    # Select methods based on coverage
    if (avg_coverage >= 95 && high_cov_clocks >= 5) {
      methods <- c("ensemble_static", "ensemble_dynamic", "ensemble_dual")
      if (verbose) cat("  Selected: Comprehensive ensemble methods (high coverage)\n")
    } else if (avg_coverage >= 85) {
      methods <- c("ensemble_static", "ensemble_dynamic")
      if (verbose) cat("  Selected: Core ensemble methods (good coverage)\n")
    } else if (avg_coverage >= 70) {
      methods <- c("ensemble_static", "all")
      if (verbose) cat("  Selected: Static ensemble + all available (moderate coverage)\n")
    } else {
      methods <- c("all")
      if (verbose) cat("  Selected: All available clocks (lower coverage)\n")
    }
  }, error = function(e) {
    methods <- c("ensemble_static")
    if (verbose) cat("  Coverage check failed - using ensemble static\n")
  })
  
  # Step 5: Run multiple prediction methods
  if (verbose) cat("\n🚀 Running age predictions...\n")
  
  all_results <- list()
  
  for (method in methods) {
    if (verbose) cat("  Running", method, "predictions...\n")
    
    tryCatch({
      result <- predict_ages(processed_data, sample_sheet, 
                           method = method,
                           output_format = "long",  # Always get long format first
                           save_results = FALSE,    # Save at the end
                           verbose = FALSE)
      
      all_results[[method]] <- result
      
      if (verbose) {
        cat("    ✓ Success:", nrow(result), "predictions from", 
            length(unique(result$epiClock)), "clocks\n")
      }
      
    }, error = function(e) {
      if (verbose) cat("    ✗ Failed:", e$message, "\n")
    })
  }
  
  # Combine all results
  if (length(all_results) > 0) {
    combined_results <- data.table::rbindlist(all_results, fill = TRUE)
    
    if (verbose) {
      cat("✓ Combined predictions complete!\n")
      cat("  Total predictions:", nrow(combined_results), "\n")
      cat("  Unique clocks:", length(unique(combined_results$epiClock)), "\n")
      cat("  Samples:", length(unique(combined_results$Basename)), "\n")
    }
  } else {
    stop("All prediction methods failed")
  }
  
  # Step 6: Format results as requested
  if (verbose) cat("\n📊 Formatting results...\n")
  
  final_results <- if (output_format == "both") {
    list(
      long = combined_results,
      wide = convert_results_format(combined_results, "wide", verbose = FALSE)
    )
  } else {
    convert_results_format(combined_results, output_format, verbose = verbose)
  }
  
  # Step 7: Save results
  if (save_results) {
    if (verbose) cat("\n💾 Saving results...\n")
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      if (verbose) cat("✓ Created output directory:", output_dir, "\n")
    }
    
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    if (!is.null(output_filename)) {
      # Custom filename provided
      if (output_format == "both") {
        # For "both" format, modify the base filename
        base_name <- tools::file_path_sans_ext(output_filename)
        ext <- tools::file_ext(output_filename)
        if (ext == "") ext <- "csv"
        
        long_file <- file.path(output_dir, paste0(base_name, "_long.", ext))
        wide_file <- file.path(output_dir, paste0(base_name, "_wide.", ext))
        
        write.csv(final_results$long, long_file, row.names = FALSE)
        write.csv(final_results$wide, wide_file, row.names = FALSE)
        
        if (verbose) {
          cat("✓ Long format saved:", basename(long_file), "\n")
          cat("✓ Wide format saved:", basename(wide_file), "\n")
        }
      } else {
        # Single format with custom filename
        output_file <- file.path(output_dir, output_filename)
        write.csv(final_results, output_file, row.names = FALSE)
        
        if (verbose) {
          file_size <- round(file.info(output_file)$size / 1024, 1)
          cat("✓ Results saved:", basename(output_file), "(", file_size, "KB)\n")
        }
      }
    } else {
      # Automatic naming with prefix
      if (output_format == "both") {
        # Save both formats
        long_file <- file.path(output_dir, paste0(output_prefix, "_", data_name, "_long_", timestamp, ".csv"))
        wide_file <- file.path(output_dir, paste0(output_prefix, "_", data_name, "_wide_", timestamp, ".csv"))
        
        write.csv(final_results$long, long_file, row.names = FALSE)
        write.csv(final_results$wide, wide_file, row.names = FALSE)
        
        if (verbose) {
          cat("✓ Long format saved:", basename(long_file), "\n")
          cat("✓ Wide format saved:", basename(wide_file), "\n")
        }
      } else {
        # Save single format
        output_file <- file.path(output_dir, paste0(output_prefix, "_", data_name, "_", output_format, "_", timestamp, ".csv"))
        write.csv(final_results, output_file, row.names = FALSE)
        
        if (verbose) {
          file_size <- round(file.info(output_file)$size / 1024, 1)
          cat("✓ Results saved:", basename(output_file), "(", file_size, "KB)\n")
        }
      }
    }
  }
  
  # Final summary
  if (verbose) {
    cat("\n🎉 EnsembleAge Prediction Complete! 🎉\n")
    cat("====================================\n")
    cat("✅ Platform automatically detected\n")
    cat("✅ Data preprocessed and optimized\n")
    cat("✅ Best prediction methods selected\n")
    cat("✅ Multiple clocks predicted\n")
    if (save_results) cat("✅ Results saved to file\n")
    cat("\n📈 Ready for analysis!\n")
  }
  
  return(final_results)
}

