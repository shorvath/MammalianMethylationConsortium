#' EnsembleAge Static Clock Predictions
#'
#' Predict age using the main EnsembleAge static clocks for mouse data.
#' This is the primary clock for mouse age predictions.
#' 
#' For the simplest interface, consider using \code{\link{predictEnsemble}} which
#' automatically handles data loading, platform detection, and preprocessing.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @param samps Data frame containing sample information (must include Age column)
#' @param efficient_loading Logical indicating whether to use efficient loading (default: TRUE)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Data frame with age predictions and acceleration values for static ensemble clocks
#' @export
#' @references Haghani, A., et al. (2025). GeroScience. DOI: 10.1007/s11357-025-01808-1
#' @examples
#' \donttest{
#' # Main mouse age prediction
#' results <- predict_ensemble_static(methylation_data, sample_info)
#' }
predict_ensemble_static <- function(dat0sesame, samps, efficient_loading = TRUE, verbose = TRUE) {
  
  if (verbose) cat("=== EnsembleAge Static Clock Predictions ===\n")
  
  # Run full prediction
  prediction_results <- predict_age_simple(dat0sesame, samps, 
                                          efficient_loading = efficient_loading, 
                                          verbose = verbose)
  
  # Get the long format results
  all_results <- attr(prediction_results, "long_format")
  if (is.null(all_results)) {
    stop("Long format results not available. Please check the prediction function.")
  }
  
  # Filter to only EnsembleAge.Static results (includes both Ensemble.Static and EnsembleAge.Static)
  static_results <- all_results %>%
    dplyr::filter(grepl("Ensemble.*\\.Static", clockFamily)) %>%
    dplyr::select(Basename, Age, epiClock, epiAge, AgeAccelation, Female, Tissue, clockFamily)
  
  if (verbose) {
    cat("EnsembleAge Static predictions complete!\n")
    cat("Predicted", nrow(static_results), "samples using", 
        length(unique(static_results$epiClock)), "static clocks\n")
  }
  
  return(static_results)
}

#' EnsembleAge Dynamic Clock Predictions
#'
#' Predict age using individual EnsembleAge dynamic clocks for mouse data.
#' These are individual clocks that can be analyzed separately.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @param samps Data frame containing sample information (must include Age column)
#' @param efficient_loading Logical indicating whether to use efficient loading (default: TRUE)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Data frame with age predictions and acceleration values for dynamic ensemble clocks
#' @export
#' @examples
#' \donttest{
#' # Individual mouse clock predictions
#' results <- predict_ensemble_dynamic(methylation_data, sample_info)
#' }
predict_ensemble_dynamic <- function(dat0sesame, samps, efficient_loading = TRUE, verbose = TRUE) {
  
  if (verbose) cat("=== EnsembleAge Dynamic Clock Predictions ===\n")
  
  # Run full prediction
  prediction_results <- predict_age_simple(dat0sesame, samps, 
                                          efficient_loading = efficient_loading, 
                                          verbose = verbose)
  
  # Get the long format results
  all_results <- attr(prediction_results, "long_format")
  if (is.null(all_results)) {
    stop("Long format results not available. Please check the prediction function.")
  }
  
  # Filter to only EnsembleAge.Dynamic results
  dynamic_results <- all_results %>%
    dplyr::filter(grepl("EnsembleAge\\.Dynamic", clockFamily)) %>%
    dplyr::select(Basename, Age, epiClock, epiAge, AgeAccelation, Female, Tissue, clockFamily)
  
  if (verbose) {
    cat("EnsembleAge Dynamic predictions complete!\n")
    cat("Predicted", nrow(dynamic_results), "samples using", 
        length(unique(dynamic_results$epiClock)), "individual dynamic clocks\n")
  }
  
  return(dynamic_results)
}

#' EnsembleDual Static Clock Predictions
#'
#' Predict age using EnsembleDual static clocks for both human and mouse data.
#' These clocks work across species and include relative age transformations.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @param samps Data frame containing sample information (must include Age column)
#' @param efficient_loading Logical indicating whether to use efficient loading (default: TRUE)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Data frame with age predictions and acceleration values for dual-species clocks
#' @export
#' @examples
#' \dontrun{
#' # Cross-species age prediction (human and mouse)
#' results <- predict_ensemble_dual_static(methylation_data, sample_info)
#' }
predict_ensemble_dual_static <- function(dat0sesame, samps, efficient_loading = TRUE, verbose = TRUE) {
  
  if (verbose) cat("=== EnsembleDual Static Clock Predictions ===\n")
  
  # Run full prediction
  prediction_results <- predict_age_simple(dat0sesame, samps, 
                                          efficient_loading = efficient_loading, 
                                          verbose = verbose)
  
  # Get the long format results
  all_results <- attr(prediction_results, "long_format")
  if (is.null(all_results)) {
    stop("Long format results not available. Please check the prediction function.")
  }
  
  # Filter to only EnsembleDualAge.Static results
  dual_results <- all_results %>%
    dplyr::filter(grepl("EnsembleDualAge\\.Static", clockFamily)) %>%
    dplyr::select(Basename, Age, epiClock, epiAge, AgeAccelation, Female, Tissue, clockFamily)
  
  if (verbose) {
    cat("EnsembleDual Static predictions complete!\n")
    cat("Predicted", nrow(dual_results), "samples using", 
        length(unique(dual_results$epiClock)), "dual-species clocks\n")
  }
  
  return(dual_results)
}

#' Original Clock Predictions
#'
#' Predict age using the original clock implementations (Ake clocks, Universal clocks, etc.).
#' These are the pre-ensemble clock methods.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @param samps Data frame containing sample information (must include Age column)
#' @param efficient_loading Logical indicating whether to use efficient loading (default: TRUE)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Data frame with age predictions and acceleration values for original clocks
#' @export
#' @examples
#' \dontrun{
#' # Original clock predictions
#' results <- predict_original_clocks(methylation_data, sample_info)
#' }
predict_original_clocks <- function(dat0sesame, samps, efficient_loading = TRUE, verbose = TRUE) {
  
  if (verbose) cat("=== Original Clock Predictions ===\n")
  
  # Run full prediction
  prediction_results <- predict_age_simple(dat0sesame, samps, 
                                          efficient_loading = efficient_loading, 
                                          verbose = verbose)
  
  # Get the long format results
  all_results <- attr(prediction_results, "long_format")
  if (is.null(all_results)) {
    stop("Long format results not available. Please check the prediction function.")
  }
  
  # Filter to exclude ensemble clocks (get original clocks)
  original_results <- all_results %>%
    dplyr::filter(!grepl("EnsembleAge\\.|EnsembleDualAge\\.", clockFamily)) %>%
    dplyr::select(Basename, Age, epiClock, epiAge, AgeAccelation, Female, Tissue, clockFamily)
  
  if (verbose) {
    cat("Original clock predictions complete!\n")
    cat("Predicted", nrow(original_results), "samples using", 
        length(unique(original_results$epiClock)), "original clocks\n")
    cat("Clock families:", paste(unique(original_results$clockFamily), collapse = ", "), "\n")
  }
  
  return(original_results)
}

#' All Clock Predictions
#'
#' Predict age using all available clocks (ensemble and original).
#' This is equivalent to the main predict_age_simple function but with clearer naming.
#' 
#' For the simplest interface, consider using \code{\link{predictEnsemble}} which
#' automatically handles data loading, platform detection, and preprocessing.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @param samps Data frame containing sample information (must include Age column)
#' @param efficient_loading Logical indicating whether to use efficient loading (default: TRUE)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Data frame with age predictions and acceleration values for all clocks
#' @export
#' @examples
#' \dontrun{
#' # All available clock predictions
#' results <- predict_all_clocks(methylation_data, sample_info)
#' }
predict_all_clocks <- function(dat0sesame, samps, efficient_loading = TRUE, verbose = TRUE) {
  
  if (verbose) cat("=== All Clock Predictions ===\n")
  
  # Run full prediction
  prediction_results <- predict_age_simple(dat0sesame, samps, 
                                          efficient_loading = efficient_loading, 
                                          verbose = verbose)
  
  # Get the long format results
  all_results <- attr(prediction_results, "long_format")
  if (is.null(all_results)) {
    stop("Long format results not available. Please check the prediction function.")
  }
  
  if (verbose) {
    cat("All clock predictions complete!\n")
    cat("Predicted", nrow(all_results), "samples using", 
        length(unique(all_results$epiClock)), "total clocks\n")
    
    # Summary by clock family
    family_summary <- all_results %>%
      dplyr::group_by(clockFamily) %>%
      dplyr::summarise(
        n_clocks = length(unique(epiClock)),
        n_predictions = dplyr::n(),
        .groups = "drop"
      )
    
    cat("Clock family summary:\n")
    print(family_summary)
  }
  
  return(all_results)
}

#' Get Clock Information
#'
#' Get information about available clocks in the package.
#' 
#' @param clock_type Character string to filter by clock type (optional)
#' @param verbose Logical indicating whether to print detailed information (default: TRUE)
#' @return Data frame with clock information
#' @export
#' @examples
#' \dontrun{
#' # Get all clock information
#' clock_info <- get_clock_info()
#' 
#' # Get only ensemble static clocks
#' ensemble_static <- get_clock_info("EnsembleAge.Static")
#' }
get_clock_info <- function(clock_type = NULL, verbose = TRUE) {
  
  # Load clock coefficients
  clock_file_path <- system.file("data", "Clock_coefficients.RDS", package = "EnsembleAge")
  if (clock_file_path == "") {
    clock_file_path <- file.path("data", "Clock_coefficients.RDS")
  }
  
  if (!file.exists(clock_file_path)) {
    stop("Clock coefficients file not found.")
  }
  
  epiclocks <- readRDS(clock_file_path)
  
  # Extract clock information
  clock_info <- data.frame(
    clockFamily = character(0),
    clockName = character(0),
    n_probes = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (family_name in names(epiclocks)) {
    family_clocks <- epiclocks[[family_name]]
    for (clock_name in names(family_clocks)) {
      clock <- family_clocks[[clock_name]]
      n_probes <- length(clock$CGid[clock$CGid != "Intercept"])
      
      clock_info <- rbind(clock_info, data.frame(
        clockFamily = family_name,
        clockName = clock_name,
        n_probes = n_probes,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Filter by clock type if specified
  if (!is.null(clock_type)) {
    clock_info <- clock_info %>%
      dplyr::filter(grepl(clock_type, clockFamily, ignore.case = TRUE))
  }
  
  if (verbose) {
    cat("Available clocks:\n")
    cat("Total families:", length(unique(clock_info$clockFamily)), "\n")
    cat("Total clocks:", nrow(clock_info), "\n")
    cat("Total unique probes:", length(unique(unlist(lapply(epiclocks, function(fam) {
      unlist(lapply(fam, function(clk) clk$CGid[clk$CGid != "Intercept"]))
    })))), "\n\n")
    
    # Summary by family
    family_summary <- clock_info %>%
      dplyr::group_by(clockFamily) %>%
      dplyr::summarise(
        n_clocks = dplyr::n(),
        total_probes = sum(n_probes),
        avg_probes = round(mean(n_probes), 1),
        .groups = "drop"
      )
    
    print(family_summary)
  }
  
  return(clock_info)
}
