#' Get required CpG probes for all clocks
#'
#' This function extracts all CpG probes required by any of the clocks
#' to enable efficient data loading.
#' 
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Character vector of required CpG IDs
#' @export
#' @examples
#' \dontrun{
#' required_probes <- get_required_clock_probes()
#' }
get_required_clock_probes <- function(verbose = TRUE) {
  
  # Load clock coefficients
  clock_file_path <- system.file("data", "Clock_coefficients.RDS", package = "EnsembleAge")
  if (clock_file_path == "") {
    clock_file_path <- file.path("data", "Clock_coefficients.RDS")
  }
  
  if (!file.exists(clock_file_path)) {
    stop("Clock coefficients file not found.")
  }
  
  epiclocks <- readRDS(clock_file_path)
  
  # Extract all required CpG IDs
  all_required_probes <- character(0)
  
  for (family_name in names(epiclocks)) {
    family_clocks <- epiclocks[[family_name]]
    for (clock_name in names(family_clocks)) {
      clock <- family_clocks[[clock_name]]
      clock_probes <- clock$CGid[clock$CGid != "Intercept"]  # Exclude intercept
      all_required_probes <- c(all_required_probes, clock_probes)
    }
  }
  
  # Get unique probes
  unique_probes <- unique(all_required_probes)
  
  if (verbose) {
    cat("Total unique CpG probes required for all clocks:", length(unique_probes), "\n")
  }
  
  return(unique_probes)
}

#' Efficiently load methylation data for clocks
#'
#' This function efficiently loads only the methylation data required for clock predictions,
#' following Amin's approach for Mammal320k to 40k conversion with proper probe filtering,
#' mapping, and missing value imputation.
#' 
#' @param dat0sesame Data frame or matrix containing methylation data
#' @param samps Data frame containing sample information
#' @param required_probes Character vector of required CpG IDs (optional, will be determined automatically)
#' @param impute_missing Logical indicating whether to impute missing probes as 0.5 (default: TRUE)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return List containing filtered data and metadata
#' @export
#' @examples
#' \dontrun{
#' efficient_data <- efficiently_load_clock_data(methylation_data, sample_info)
#' }
efficiently_load_clock_data <- function(dat0sesame, samps, required_probes = NULL, 
                                       impute_missing = TRUE, verbose = TRUE) {
  
  if (verbose) cat("Starting efficient data loading for clocks...\n")
  
  # Get required probes if not provided
  if (is.null(required_probes)) {
    if (verbose) cat("Determining required CpG probes...\n")
    required_probes <- get_required_clock_probes(verbose = verbose)
  }
  
  # First, detect and fix data orientation
  if (verbose) cat("Checking data orientation...\n")
  orientation_result <- detect_and_fix_orientation(dat0sesame, samps, verbose = verbose)
  dat0sesame <- orientation_result$data
  
  # Check if this is Mammal320k data that needs probe mapping
  if (verbose) cat("Checking for probe mapping requirements...\n")
  platform_check <- detect_platform(dat0sesame)
  original_probe_count <- nrow(dat0sesame)
  
  if (platform_check == "Mammal320k" || 
      sum(grepl("_[A-Z]+[0-9]*$", dat0sesame$CGid[1:min(100, nrow(dat0sesame))])) > 10) {
    
    # Mammal320k processing is now handled in preprocess_methylation_data()
    # Skip redundant processing here to avoid conflicts
    if (verbose) cat("Mammal320k detected but processing already done in preprocessing step...\n")
    
    # Just filter to required probes (processing already done in preprocessing)
    available_probes <- intersect(required_probes, dat0sesame$CGid)
    dat0sesame <- dat0sesame %>%
      dplyr::filter(CGid %in% available_probes)
    
    if (verbose) {
      cat("Filtered to", length(available_probes), "required probes out of", length(required_probes), "\n")
    }
  } else {
    # For other platforms (Mammal40k, Human), just filter to required probes
    if (verbose) cat("Filtering data to required CpG probes for", platform_check, "platform...\n")
    available_probes <- intersect(required_probes, dat0sesame$CGid)
    
    if (verbose) {
      cat("Found", length(available_probes), "required probes out of", length(required_probes), "in data\n")
    }
    
    dat0sesame <- dat0sesame %>%
      dplyr::filter(CGid %in% available_probes)
  }
  
  # Check which required probes are missing and impute as 0.5 (following line 228)
  current_probes <- dat0sesame$CGid
  missing_probes <- setdiff(required_probes, current_probes)
  
  if (length(missing_probes) > 0 && impute_missing) {
    if (verbose) {
      cat("Imputing", length(missing_probes), "missing probes with 0.5 (following Amin's approach)...\n")
    }
    
    # Create imputed data for missing probes
    sample_cols <- names(dat0sesame)[names(dat0sesame) != "CGid"]
    n_samples <- length(sample_cols)
    
    missing_data <- data.frame(
      CGid = missing_probes,
      matrix(0.5, nrow = length(missing_probes), ncol = n_samples),
      stringsAsFactors = FALSE
    )
    names(missing_data) <- c("CGid", sample_cols)
    
    # Combine with existing data
    dat0sesame <- rbind(dat0sesame, missing_data)
  }
  
  if (verbose) {
    cat("Efficient loading complete:\n")
    cat("  Original probes:", original_probe_count, "\n")
    cat("  Required probes:", length(required_probes), "\n")
    cat("  Available probes:", length(current_probes), "\n")
    cat("  Missing probes:", length(missing_probes), ifelse(impute_missing, "(imputed as 0.5)", "(not imputed)"), "\n")
    cat("  Final probes:", nrow(dat0sesame), "\n")
    cat("  Compression ratio:", round(nrow(dat0sesame) / original_probe_count, 3), "\n")
  }
  
  return(list(
    data = dat0sesame,
    required_probes = required_probes,
    available_probes = current_probes,
    missing_probes = missing_probes,
    imputed = impute_missing && length(missing_probes) > 0,
    platform = platform_check,
    compression_ratio = round(nrow(dat0sesame) / original_probe_count, 3)
  ))
}

#' Load sample sheet template
#'
#' This function creates and saves a sample sheet template that users can fill out.
#' 
#' @param file_path Path where to save the template (default: "sample_sheet_template.csv")
#' @param sample_names Optional character vector of sample names
#' @param n_samples Number of samples if sample_names not provided (default: 20)
#' @return Data frame with the template (also saved to file)
#' @export
#' @examples
#' # Create and save template
#' template <- load_sample_sheet_template("my_samples_template.csv", n_samples = 50)
load_sample_sheet_template <- function(file_path = "sample_sheet_template.csv", 
                                      sample_names = NULL, n_samples = 20) {
  
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", sprintf("%03d", 1:n_samples))
  } else {
    n_samples <- length(sample_names)
  }
  
  # Create template with only prediction-required columns
  template <- data.frame(
    Basename = sample_names,
    Age = rep(NA_real_, n_samples),
    Female = rep(NA_integer_, n_samples),
    SpeciesLatinName = rep("Mus musculus", n_samples),
    Tissue = rep(NA_character_, n_samples),
    stringsAsFactors = FALSE
  )
  
  # Add instruction rows at the top
  instructions <- data.frame(
    Basename = c(
      "# INSTRUCTIONS: Fill in the information below and delete these instruction rows",
      "# REQUIRED: Basename - Sample identifiers matching your methylation data",
      "# REQUIRED: Age - Chronological age in years (e.g., 1.5 for 18-month-old mouse)",
      "# OPTIONAL: Female - 1 for female, 0 for male (defaults to NA if missing)",
      "# OPTIONAL: SpeciesLatinName - 'Mus musculus', 'Rattus norvegicus', 'Homo sapiens' (defaults to 'Mus musculus')",
      "# OPTIONAL: Tissue - 'Blood', 'Liver', 'Brain', 'Heart', etc. (defaults to 'Unknown' if missing)",
      "# DELETE ALL INSTRUCTION ROWS BEFORE USING THE TEMPLATE"
    ),
    Age = rep("", 7),
    Female = rep("", 7),
    SpeciesLatinName = rep("", 7),
    Tissue = rep("", 7),
    stringsAsFactors = FALSE
  )
  
  # Combine instructions and template
  full_template <- rbind(instructions, template)
  
  # Save to file
  write.csv(full_template, file_path, row.names = FALSE)
  
  cat("Sample sheet template created with", n_samples, "samples\n")
  cat("Template saved to:", file_path, "\n")
  cat("Please:\n")
  cat("1. Fill in the Basename column with your exact sample names\n")
  cat("2. Fill in the Age column with chronological ages in years\n")
  cat("3. Fill in other columns as needed\n") 
  cat("4. Delete the instruction rows (first 9 rows)\n")
  cat("5. Save the file and use it with EnsembleAge\n")
  
  return(full_template)
}
