#' Detect methylation array platform
#'
#' This function detects the methylation array platform based on the number of CpG sites
#' and specific probe patterns in the data.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @return Character string indicating the detected platform: "Mammal320k", "Mammal40k", "Human450k", "HumanEPIC", or "Unknown"
#' @export
#' @examples
#' \dontrun{
#' # Example with methylation data
#' platform <- detect_platform(methylation_data)
#' }
detect_platform <- function(dat0sesame) {
  
  # Handle both raw data matrices and processed data frames
  if (is.matrix(dat0sesame)) {
    # Raw matrix data - check dimensions and probe patterns
    n_rows <- nrow(dat0sesame)
    n_cols <- ncol(dat0sesame)
    
    # Check for mammal320k pattern (many probes as columns with _BC21 suffixes)
    if (n_cols > 250000) {
      col_suffix_pattern <- sum(grepl("_[A-Z]+[0-9]*$", colnames(dat0sesame)[1:min(1000, n_cols)]))
      if (col_suffix_pattern > 100) {
        return("Mammal320k")
      }
    }
    
    # Check for mammal320k pattern (many probes as rows with _BC21 suffixes)  
    if (n_rows > 250000) {
      row_suffix_pattern <- sum(grepl("_[A-Z]+[0-9]*$", rownames(dat0sesame)[1:min(1000, n_rows)]))
      if (row_suffix_pattern > 100) {
        return("Mammal320k")
      }
    }
    
    # For matrices, we can't proceed without CGid column
    stop("Matrix data detected. Please convert to data.frame with CGid column or use preprocess_methylation_data() first.")
  }
  
  if (!"CGid" %in% names(dat0sesame)) {
    # Check if first column might be probe IDs
    first_col <- names(dat0sesame)[1]
    probe_like <- sum(grepl("^cg[0-9]", dat0sesame[[first_col]][1:min(100, nrow(dat0sesame))]))
    
    if (probe_like > 50) {
      warning("No 'CGid' column found, but first column appears to contain probe IDs. Consider renaming to 'CGid'.")
      # Temporarily use first column for detection
      probe_ids <- dat0sesame[[first_col]]
    } else {
      stop("Data must contain a 'CGid' column with CpG identifiers")
    }
  } else {
    probe_ids <- dat0sesame$CGid
  }
  
  n_probes <- nrow(dat0sesame)
  
  # Check for mammal320k suffix pattern first (priority check)
  suffix_pattern <- sum(grepl("_[A-Z]+[0-9]*$", probe_ids[1:min(1000, length(probe_ids))]))
  if (suffix_pattern > 100) {
    return("Mammal320k")
  }
  
  # Standard platform detection based on probe counts
  if (n_probes > 250000) {
    return("Mammal320k")
  } else if (n_probes > 30000 && n_probes < 50000) {
    return("Mammal40k")
  } else if (n_probes > 400000 && n_probes < 500000) {
    return("Human450k")
  } else if (n_probes > 800000) {
    return("HumanEPIC")
  } else {
    # Additional checks based on probe naming patterns
    if (any(grepl("^cg", probe_ids)) && any(grepl("^ch", probe_ids))) {
      if (n_probes > 800000) {
        return("HumanEPIC")
      } else {
        return("Human450k")
      }
    } else if (any(grepl("^cg", probe_ids))) {
      return("Mammal40k")
    } else {
      return("Unknown")
    }
  }
}

#' Validate methylation data format
#'
#' This function validates that the methylation data is in the correct format
#' for EnsembleAge predictions.
#' 
#' @param dat0sesame Data frame containing methylation data
#' @param samps Data frame containing sample information
#' @return List with validation results and suggestions
#' @export
#' @examples
#' \dontrun{
#' validation <- validate_data_format(methylation_data, sample_info)
#' if (!validation$valid) {
#'   stop(validation$message)
#' }
#' }
validate_data_format <- function(dat0sesame, samps) {
  issues <- character(0)
  suggestions <- character(0)
  
  # Check dat0sesame format
  if (!"CGid" %in% names(dat0sesame)) {
    issues <- c(issues, "dat0sesame must contain a 'CGid' column")
    suggestions <- c(suggestions, "Add a 'CGid' column with CpG identifiers")
  }
  
  if (ncol(dat0sesame) < 2) {
    issues <- c(issues, "dat0sesame must contain sample columns in addition to CGid")
  }
  
  # Check for numeric methylation values
  numeric_cols <- sapply(dat0sesame[, !names(dat0sesame) %in% "CGid"], is.numeric)
  if (!all(numeric_cols)) {
    issues <- c(issues, "All sample columns must contain numeric methylation values")
  }
  
  # Check value ranges
  if (any(sapply(dat0sesame[, !names(dat0sesame) %in% "CGid"], function(x) any(x < 0 | x > 1, na.rm = TRUE)))) {
    issues <- c(issues, "Methylation values should be between 0 and 1 (beta values)")
    suggestions <- c(suggestions, "Convert M-values to beta values if necessary")
  }
  
  # Check samps format
  if (!"Basename" %in% names(samps)) {
    issues <- c(issues, "samps must contain a 'Basename' column with sample identifiers")
  }
  
  # Check if sample names match
  sample_cols <- names(dat0sesame)[!names(dat0sesame) %in% "CGid"]
  if ("Basename" %in% names(samps)) {
    missing_samples <- setdiff(samps$Basename, sample_cols)
    if (length(missing_samples) > 0) {
      issues <- c(issues, paste("Some samples in samps$Basename are not found in dat0sesame:", 
                                paste(missing_samples[1:min(5, length(missing_samples))], collapse = ", ")))
    }
  }
  
  # Platform detection
  platform <- detect_platform(dat0sesame)
  
  return(list(
    valid = length(issues) == 0,
    platform = platform,
    issues = issues,
    suggestions = suggestions,
    n_probes = nrow(dat0sesame),
    n_samples = ncol(dat0sesame) - 1
  ))
}

#' Prepare sample sheet with default values
#'
#' This function prepares the sample sheet by adding default values for missing columns
#' required by the EnsembleAge prediction functions.
#' 
#' @param samps Data frame containing sample information
#' @param default_species Character string for default species (default: "Mus musculus")
#' @param default_age Numeric value for default age (default: 0)
#' @return Data frame with standardized sample information
#' @export
#' @importFrom dplyr mutate
#' @examples
#' \dontrun{
#' sample_info <- data.frame(Basename = c("sample1", "sample2"))
#' prepared_samples <- prepare_sample_sheet(sample_info)
#' }
prepare_sample_sheet <- function(samps, verbose = TRUE) {
  
  missing_vars <- character(0)
  defaults_applied <- character(0)
  
  # Ensure Basename column exists
  if (!"Basename" %in% names(samps)) {
    if ("Sample_Name" %in% names(samps)) {
      samps$Basename <- samps$Sample_Name
      if (verbose) cat("Using 'Sample_Name' column as 'Basename'\n")
    } else {
      stop("Sample sheet must contain either 'Basename' or 'Sample_Name' column")
    }
  }
  
  # Ensure Age column exists - this is truly required for clock predictions
  if (!"Age" %in% names(samps)) {
    stop("Sample sheet must contain 'Age' column with chronological ages in years")
  }
  
  # Add missing optional columns with defaults and track what was added
  if (!"Female" %in% names(samps)) {
    samps <- samps %>% mutate(Female = NA)
    missing_vars <- c(missing_vars, "Female")
    defaults_applied <- c(defaults_applied, "Female = NA (sex unknown)")
  }
  
  if (!"Tissue" %in% names(samps)) {
    samps <- samps %>% mutate(Tissue = "Unknown")
    missing_vars <- c(missing_vars, "Tissue")
    defaults_applied <- c(defaults_applied, "Tissue = 'Unknown'")
  }
  
  if (!"SpeciesLatinName" %in% names(samps)) {
    # Try to detect species based on context - for now default to mouse
    # Future enhancement: could detect based on data dimensions or other clues
    samps <- samps %>% mutate(SpeciesLatinName = "Mus musculus")
    missing_vars <- c(missing_vars, "SpeciesLatinName")
    defaults_applied <- c(defaults_applied, "SpeciesLatinName = 'Mus musculus'")
  }
  
  # Notify user about defaults applied
  if (length(missing_vars) > 0 && verbose) {
    cat("Missing variables detected. Applied defaults:\n")
    for (default in defaults_applied) {
      cat("  -", default, "\n")
    }
    cat("For better predictions, consider providing these variables in your sample sheet.\n")
  }
  
  # Convert data types with error handling
  tryCatch({
    samps$Age <- as.numeric(samps$Age)
    if (any(is.na(samps$Age))) {
      warning("Some Age values could not be converted to numeric. Check your data.")
    }
  }, error = function(e) {
    stop("Error converting Age column to numeric: ", e$message)
  })
  
  # Handle Female column conversion
  if (!all(is.na(samps$Female))) {
    tryCatch({
      samps$Female <- as.numeric(samps$Female)
      # Validate Female values (should be 0, 1, or NA)
      invalid_female <- samps$Female[!is.na(samps$Female) & !samps$Female %in% c(0, 1)]
      if (length(invalid_female) > 0) {
        warning("Female column should contain only 0 (male), 1 (female), or NA. Invalid values found.")
      }
    }, error = function(e) {
      warning("Error converting Female column: ", e$message, ". Setting to NA.")
      samps$Female <- NA
    })
  }
  
  # Ensure character columns
  samps$Tissue <- as.character(samps$Tissue)
  samps$SpeciesLatinName <- as.character(samps$SpeciesLatinName)
  samps$Basename <- as.character(samps$Basename)
  
  # Validate that we have samples
  if (nrow(samps) == 0) {
    stop("Sample sheet is empty")
  }
  
  # Check for duplicate Basenames
  if (any(duplicated(samps$Basename))) {
    duplicates <- samps$Basename[duplicated(samps$Basename)]
    warning("Duplicate sample names found: ", paste(unique(duplicates), collapse = ", "))
  }
  
  if (verbose) {
    cat("Sample sheet prepared:", nrow(samps), "samples\n")
  }
  
  return(samps)
}

#' Get available clocks for a platform
#'
#' This function returns the available clocks for a specific methylation platform.
#' 
#' @param platform Character string indicating the platform
#' @return Character vector of available clock names
#' @export
#' @examples
#' available_clocks <- get_available_clocks("Mammal40k")
get_available_clocks <- function(platform) {
  clock_file_path <- system.file("data", "Clock_coefficients.RDS", package = "EnsembleAge")
  if (clock_file_path == "") {
    clock_file_path <- file.path("data", "Clock_coefficients.RDS")
  }
  
  if (!file.exists(clock_file_path)) {
    warning("Clock coefficients file not found. Please ensure the package is properly installed.")
    return(character(0))
  }
  
  epiclocks <- readRDS(clock_file_path)
  
  # Return all available clock families and individual clocks
  all_clocks <- character(0)
  for (family_name in names(epiclocks)) {
    family_clocks <- epiclocks[[family_name]]
    clock_names <- paste(family_name, names(family_clocks), sep = ".")
    all_clocks <- c(all_clocks, clock_names)
  }
  
  return(all_clocks)
}

#' Check probe coverage for platform
#'
#' This function checks what percentage of required probes are available in the input data
#' for each clock type.
#' 
#' @param dat0sesame Data frame containing methylation data with CGid column
#' @return Data frame with clock names and probe coverage percentages
#' @export
#' @examples
#' \dontrun{
#' coverage <- check_probe_coverage(methylation_data)
#' }
check_probe_coverage <- function(dat0sesame) {
  if (!"CGid" %in% names(dat0sesame)) {
    stop("Data must contain a 'CGid' column")
  }
  
  clock_file_path <- system.file("data", "Clock_coefficients.RDS", package = "EnsembleAge")
  if (clock_file_path == "") {
    clock_file_path <- file.path("data", "Clock_coefficients.RDS")
  }
  
  if (!file.exists(clock_file_path)) {
    warning("Clock coefficients file not found.")
    return(data.frame(Clock = character(0), Coverage = numeric(0)))
  }
  
  epiclocks <- readRDS(clock_file_path)
  available_probes <- dat0sesame$CGid
  
  coverage_results <- data.frame(
    Clock = character(0),
    Required_Probes = integer(0),
    Available_Probes = integer(0),
    Coverage_Percent = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (family_name in names(epiclocks)) {
    family_clocks <- epiclocks[[family_name]]
    for (clock_name in names(family_clocks)) {
      clock <- family_clocks[[clock_name]]
      required_probes <- clock$CGid
      available_in_data <- sum(required_probes %in% available_probes)
      coverage_percent <- round(available_in_data / length(required_probes) * 100, 1)
      
      coverage_results <- rbind(coverage_results, data.frame(
        Clock = paste(family_name, clock_name, sep = "."),
        Required_Probes = length(required_probes),
        Available_Probes = available_in_data,
        Coverage_Percent = coverage_percent,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(coverage_results[order(coverage_results$Coverage_Percent, decreasing = TRUE), ])
}

#' Process Mammal320k data with complete workflow
#'
#' This function processes Mammal320k methylation data following the complete workflow:
#' 1. Maps probe IDs to CpG IDs using revised annotation
#' 2. Adds missing probes with median imputation 
#' 3. Handles species-specific missing probe imputation
#' 
#' @param dat0sesame Data frame or matrix containing methylation data with probe IDs
#' @param species Character string indicating species ("mouse", "rat", "human") for missing probe imputation (default: "mouse")
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return Data frame with CpG IDs mapped, missing probes imputed, ready for clock predictions
#' @export
#' @examples
#' \dontrun{
#' processed_data <- process_mammal320k_data(methylation_data, species = "mouse")
#' }
process_mammal320k_data <- function(mammal320Data, sample_sheet, species = "mouse", verbose = TRUE) {
  
  if (verbose) cat("=== Processing Mammal320k Data (Following Amin's Workflow) ===\n")
  
  # Step 1: Load geneMap320k (prioritize local file over package file)
  annotation_file <- file.path("data", "Mus_musculus_Mammalian_320k_mm10_Amin_V10.RDS")
  if (!file.exists(annotation_file)) {
    annotation_file <- system.file("data", "Mus_musculus_Mammalian_320k_mm10_Amin_V10.RDS", 
                                   package = "EnsembleAge")
  }
  
  if (!file.exists(annotation_file)) {
    stop("Mammal320k annotation file not found: ", annotation_file)
  }
  
  geneMap320_local <- readRDS(annotation_file)
  if (verbose) cat("Loaded geneMap320 with", nrow(geneMap320_local), "probes\n")
  
  # Step 2: Load mammalian array reference and median imputation data
  # Load mammalian array (40k reference) - prioritize local file
  mammalian_file <- file.path("data", "Mus_musculus_Mammalian_40k_mm10_Amin_V10.RDS")
  if (!file.exists(mammalian_file)) {
    mammalian_file <- system.file("data", "Mus_musculus_Mammalian_40k_mm10_Amin_V10.RDS", package = "EnsembleAge")
  }
  
  if (!file.exists(mammalian_file)) {
    stop("Mammalian 40k reference file not found: ", mammalian_file)
  }
  
  mammalianArray <- readRDS(mammalian_file)
  if (verbose) cat("Loaded mammalianArray with", nrow(mammalianArray), "probes\n")
  
  # Load median imputation data based on species
  if (species == "mouse") {
    medians_file <- file.path("data", "gold_mouse_medians.csv")
    if (!file.exists(medians_file)) {
      medians_file <- system.file("data", "gold_mouse_medians.csv", package = "EnsembleAge")
    }
    
    if (file.exists(medians_file)) {
      medians <- read.csv(medians_file)
      if (verbose) cat("Loaded mouse median imputation data\n")
    } else {
      warning("Mouse median file not found")
      medians <- NULL
    }
    
  } else if (species == "rat") {
    medians_file <- file.path("data", "gold_Brownrat_medians.csv")
    if (!file.exists(medians_file)) {
      medians_file <- system.file("data", "gold_Brownrat_medians.csv", package = "EnsembleAge")
    }
    
    if (file.exists(medians_file)) {
      medianRat <- read.csv(medians_file)
      if (verbose) cat("Loaded rat median imputation data\n")
    } else {
      warning("Rat median file not found")
      medianRat <- NULL
    }
    
  } else if (species == "human") {
    medians_file <- file.path("data", "gold_human_medians.RDS")
    if (!file.exists(medians_file)) {
      medians_file <- system.file("data", "gold_human_medians.RDS", package = "EnsembleAge")
    }
    
    if (file.exists(medians_file)) {
      medianHuman <- readRDS(medians_file) %>% 
        dplyr::select(CpG, human_median_pantissue) %>% 
        setNames(c("CGid", "mouse_median"))
      if (verbose) cat("Loaded human median imputation data\n")
    } else {
      warning("Human median file not found")
      medianHuman <- NULL
    }
  } else {
    warning("Unknown species: ", species, ". Using default 0.5 imputation")
    medians <- NULL
  }
  
  # Step 3: Following your exact workflow
  normalized_betas_sesame2 <- mammal320Data
  
  if (verbose) cat("Original data dimensions:", dim(normalized_betas_sesame2), "\n")
  
  # Step 4: Your exact workflow - select geneMap320 probes, transpose, map
  # Check data format and find probe IDs
  if ("CGid" %in% names(normalized_betas_sesame2)) {
    # Data already processed with CGid column
    if (verbose) cat("Data already has CGid column, skipping processing...\n")
    return(normalized_betas_sesame2)
  } else if (nrow(normalized_betas_sesame2) > ncol(normalized_betas_sesame2)) {
    # Data is probes x samples (already transposed)
    probe_ids <- rownames(normalized_betas_sesame2)
    if (verbose) cat("Data appears to be probes x samples format\n")
  } else {
    # Data is samples x probes (original format)
    probe_ids <- colnames(normalized_betas_sesame2)
    if (verbose) cat("Data appears to be samples x probes format\n")
  }
  
  # First check which probes are available
  available_probes <- intersect(geneMap320_local$Probe_ID, probe_ids)
  
  if (verbose) {
    cat("Available probes from geneMap320:", length(available_probes), "out of", nrow(geneMap320_local), "\n")
  }
  
  if (length(available_probes) == 0) {
    stop("No matching probes found between geneMap320 and data")
  }
  
  # Process data based on orientation
  if (nrow(normalized_betas_sesame2) > ncol(normalized_betas_sesame2)) {
    # Data is probes x samples (already transposed)
    if (verbose) cat("Processing probes x samples data...\n")
    dat <- as.data.frame(normalized_betas_sesame2)[available_probes, ] %>% 
      tibble::rownames_to_column(var = "Probe_ID") %>% 
      left_join(dplyr::select(.data = geneMap320_local, Probe_ID, CGid), by = "Probe_ID") %>% 
      dplyr::select(-Probe_ID) %>% 
      relocate(CGid, 1)
  } else {
    # Data is samples x probes (original format) - use your exact workflow
    if (verbose) cat("Processing samples x probes data (original format)...\n")
    dat <- as.data.frame(normalized_betas_sesame2) %>% 
      dplyr::select(all_of(available_probes)) %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Probe_ID") %>% 
      left_join(dplyr::select(.data = geneMap320_local, Probe_ID, CGid), by = "Probe_ID") %>% 
      dplyr::select(-Probe_ID) %>% 
      relocate(CGid, 1)
  }
  
  if (verbose) {
    cat("After probe mapping:", nrow(dat), "CpG sites\n")
    cat("Sample columns:", ncol(dat) - 1, "\n")
  }
  
  # Step 5: Add missing probes following Amin's exact code
  if (verbose) cat("Creating missing probes using mammalianArray reference...\n")
  
  # Following Amin's exact code for missing probes
  if (species == "mouse" && !is.null(medians)) {
    missingProbes <- mammalianArray %>% 
      filter(!CGid %in% geneMap320_local$CGid) %>% 
      left_join(medians, by = c("CGid" = "CpG")) %>% 
      filter(!is.na(mouse_median))
    if (verbose) cat("Mouse missing probes:", nrow(missingProbes), "\n")
  }
  
  if (species == "rat" && !is.null(medianRat)) {
    missingProbesRat <- mammalianArray %>% 
      filter(!CGid %in% geneMap320_local$CGid) %>% 
      left_join(medianRat, by = c("CGid" = "CpG")) %>% 
      filter(!is.na(mouse_median))
    if (verbose) cat("Rat missing probes:", nrow(missingProbesRat), "\n")
  }
  
  if (species == "human" && !is.null(medianHuman)) {
    missingProbesHuman <- mammalianArray %>% 
      filter(!CGid %in% geneMap320_local$CGid) %>% 
      left_join(medianHuman, by = c("CGid" = "CGid")) %>% 
      filter(!is.na(mouse_median))
    if (verbose) cat("Human missing probes:", nrow(missingProbesHuman), "\n")
  }
  
  # Check if we have missing probes to add
  has_missing_probes <- FALSE
  if (species == "mouse" && exists("missingProbes") && nrow(missingProbes) > 0) {
    has_missing_probes <- TRUE
  } else if (species == "rat" && exists("missingProbesRat") && nrow(missingProbesRat) > 0) {
    has_missing_probes <- TRUE
  } else if (species == "human" && exists("missingProbesHuman") && nrow(missingProbesHuman) > 0) {
    has_missing_probes <- TRUE
  }
  
  if (has_missing_probes) {
    if (verbose) cat("Adding missing probes with median imputation...\n")
    
    # Following Amin's exact code for creating dat2
    if (ncol(dat) > 1) {
      dat2 <- rbindlist(lapply(2:ncol(dat), function(y) {
        n <- names(dat)[y]
        if (species == "rat" && exists("missingProbesRat")) {
          a <- missingProbesRat %>% mutate(var = n)
        } else if (species == "human" && exists("missingProbesHuman")) {
          a <- missingProbesHuman %>% mutate(var = n)
        } else if (exists("missingProbes")) {
          a <- missingProbes %>% mutate(var = n)
        } else {
          # Fallback - create empty data frame with right structure
          a <- data.frame(CGid = character(0), mouse_median = numeric(0), var = character(0))
        }
        return(a)
      })) %>% 
        spread(key = "var", value = "mouse_median") %>% 
        dplyr::select(names(dat))
      
      dat <- bind_rows(dat, dat2)
      if (verbose) cat("Final data after imputation:", nrow(dat), "CpG sites\n")
    } else {
      if (verbose) cat("Warning: No sample columns for imputation\n")
    }
  } else {
    if (verbose) cat("No missing probes to add\n")
  }
  
  # Step 6: Final transformation (following your exact code)
  # dt <- dat %>% tibble::column_to_rownames("CGid") %>% dplyr::select(sample_sheet$idat_name) %>% as.matrix() %>% t() 
  # dt[is.na(dt)] <- 0.5
  
  # For EnsembleAge compatibility, return as data.frame with CGid column
  # Handle NAs
  sample_cols <- names(dat)[names(dat) != "CGid"]
  for (col in sample_cols) {
    dat[[col]][is.na(dat[[col]])] <- 0.5
  }
  
  if (verbose) {
    cat("=== Mammal320k Processing Complete (Amin's Workflow) ===\n")
    cat("Final dimensions:", dim(dat), "\n")
    cat("Ready for clock predictions\n")
  }
  
  return(dat)
}
