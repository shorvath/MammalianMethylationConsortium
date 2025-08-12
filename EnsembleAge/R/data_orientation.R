#' Detect and fix data orientation
#'
#' This function detects whether methylation data has probes as rows or columns
#' and converts it to the standard format (probes as rows, samples as columns).
#' 
#' @param dat0sesame Data frame or matrix containing methylation data
#' @param samps Optional data frame containing sample information with Basename column
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' @return List containing the corrected data and metadata about the transformation
#' @export
#' @examples
#' \dontrun{
#' # Detect and fix orientation
#' fixed_data <- detect_and_fix_orientation(methylation_data, sample_info)
#' corrected_data <- fixed_data$data
#' }
detect_and_fix_orientation <- function(dat0sesame, samps = NULL, verbose = TRUE) {
  
  original_class <- class(dat0sesame)
  original_dims <- dim(dat0sesame)
  
  if (verbose) {
    cat("Original data class:", class(dat0sesame), "\n")
    cat("Original dimensions:", original_dims[1], "x", original_dims[2], "\n")
  }
  
  # Convert matrix to data frame if needed
  if (is.matrix(dat0sesame)) {
    if (verbose) cat("Converting matrix to data frame...\n")
    
    # Check if row names look like CpG IDs
    row_names <- rownames(dat0sesame)
    col_names <- colnames(dat0sesame)
    
    # CpG pattern matching
    row_cg_pattern <- sum(grepl("^cg|^ch|^rs", row_names[1:min(100, length(row_names))], ignore.case = TRUE))
    col_cg_pattern <- sum(grepl("^cg|^ch|^rs", col_names[1:min(100, length(col_names))], ignore.case = TRUE))
    
    if (verbose) {
      cat("Row names with CpG pattern:", row_cg_pattern, "/", min(100, length(row_names)), "\n")
      cat("Column names with CpG pattern:", col_cg_pattern, "/", min(100, length(col_names)), "\n")
    }
    
    # Decide orientation based on CpG patterns and dimensions
    if (col_cg_pattern > row_cg_pattern) {
      # Columns are CpGs, rows are samples - need to transpose
      if (verbose) cat("Detected: Samples as rows, CpGs as columns. Transposing...\n")
      dat0sesame <- t(dat0sesame)
      transposed <- TRUE
    } else {
      # Rows are CpGs, columns are samples - correct orientation
      if (verbose) cat("Detected: CpGs as rows, samples as columns. No transpose needed.\n")
      transposed <- FALSE
    }
    
    # Convert to data frame with CGid column
    dat0sesame <- dat0sesame %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("CGid")
      
  } else if (is.data.frame(dat0sesame)) {
    transposed <- FALSE
    
    # Check if CGid column exists
    if (!"CGid" %in% names(dat0sesame)) {
      if (verbose) cat("No CGid column found. Checking if row names are CpG IDs...\n")
      
      # Check if row names look like CpG IDs
      row_names <- rownames(dat0sesame)
      if (sum(grepl("^cg|^ch|^rs", row_names[1:min(100, length(row_names))], ignore.case = TRUE)) > 50) {
        if (verbose) cat("Row names appear to be CpG IDs. Adding CGid column...\n")
        dat0sesame <- dat0sesame %>% tibble::rownames_to_column("CGid")
      } else {
        # Check if we need to transpose
        first_col_values <- dat0sesame[1:min(100, nrow(dat0sesame)), 1]
        if (sum(grepl("^cg|^ch|^rs", first_col_values, ignore.case = TRUE)) > 50) {
          if (verbose) cat("First column appears to contain CpG IDs...\n")
          names(dat0sesame)[1] <- "CGid"
        } else {
          # Check for human data pattern: many rows (>100k), few columns
          if (nrow(dat0sesame) > 100000 && ncol(dat0sesame) < 1000) {
            if (verbose) cat("Detected likely human array data (", nrow(dat0sesame), " probes). Creating CGid column from row indices...\n")
            # For human data without proper CpG IDs, create generic IDs
            # This will require probe imputation later, but allows processing to continue
            dat0sesame <- dat0sesame %>% 
              tibble::rownames_to_column("CGid") %>%
              mutate(CGid = ifelse(grepl("^probe_", CGid), CGid, paste0("probe_", CGid)))
            if (verbose) cat("Note: Using generic probe IDs. Most clocks will rely on imputation.\n")
            
          } else if (nrow(dat0sesame) < ncol(dat0sesame) && ncol(dat0sesame) > 10000) {
            if (verbose) cat("Data appears to be transposed (more columns than rows). Transposing...\n")
            dat0sesame <- dat0sesame %>%
              tibble::rownames_to_column("temp_id") %>%
              tibble::column_to_rownames("temp_id") %>%
              t() %>%
              as.data.frame() %>%
              tibble::rownames_to_column("CGid")
            transposed <- TRUE
          } else {
            stop("Cannot determine data orientation. Please ensure data has a 'CGid' column or CpG IDs as row names.")
          }
        }
      }
    } else {
      if (verbose) cat("CGid column found.\n")
    }
  }
  
  # Additional validation if sample info is provided
  if (!is.null(samps) && "Basename" %in% names(samps)) {
    sample_cols <- names(dat0sesame)[!names(dat0sesame) %in% "CGid"]
    matching_samples <- intersect(samps$Basename, sample_cols)
    
    if (verbose) {
      cat("Sample matching: Found", length(matching_samples), "matching samples out of", 
          nrow(samps), "in sample sheet and", length(sample_cols), "in data\n")
    }
    
    # If very few matches, might still be wrong orientation
    if (length(matching_samples) < 0.1 * nrow(samps) && !transposed) {
      if (verbose) cat("Very few matching samples. Data might still need transposing...\n")
      
      # Try checking if sample names match CpG column values
      cgid_values <- dat0sesame$CGid[1:min(100, nrow(dat0sesame))]
      sample_matches_in_cgid <- sum(samps$Basename %in% cgid_values)
      
      if (sample_matches_in_cgid > length(matching_samples)) {
        if (verbose) cat("Found more sample matches in CGid column. Attempting transpose...\n")
        
        # Transpose the data
        dat0sesame_t <- dat0sesame %>%
          tibble::column_to_rownames("CGid") %>%
          t() %>%
          as.data.frame() %>%
          tibble::rownames_to_column("CGid")
        
        # Check matches after transpose
        sample_cols_t <- names(dat0sesame_t)[!names(dat0sesame_t) %in% "CGid"]
        matching_samples_t <- intersect(samps$Basename, sample_cols_t)
        
        if (length(matching_samples_t) > length(matching_samples)) {
          if (verbose) cat("Transpose improved sample matching. Using transposed data.\n")
          dat0sesame <- dat0sesame_t
          transposed <- !transposed
        }
      }
    }
  }
  
  final_dims <- dim(dat0sesame)
  if (verbose) {
    cat("Final dimensions:", final_dims[1], "probes x", final_dims[2] - 1, "samples\n")
    cat("Transposed:", transposed, "\n")
  }
  
  return(list(
    data = dat0sesame,
    transposed = transposed,
    original_class = original_class,
    original_dims = original_dims,
    final_dims = final_dims
  ))
}

#' Create sample sheet template
#'
#' This function creates a template sample sheet that users can fill out
#' with their sample information for use with EnsembleAge predictions.
#' 
#' @param sample_names Character vector of sample names (optional)
#' @param n_samples Number of samples if sample_names not provided (default: 10)
#' @param file_path Path to save the template CSV file (optional)
#' @return Data frame with sample sheet template
#' @export
#' @examples
#' # Create template for specific sample names
#' my_samples <- c("sample1", "sample2", "sample3")
#' template <- create_sample_sheet_template(sample_names = my_samples)
#' 
#' # Create template and save to file
#' create_sample_sheet_template(n_samples = 20, file_path = "sample_sheet_template.csv")
create_sample_sheet_template <- function(sample_names = NULL, n_samples = 10, file_path = NULL) {
  
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", sprintf("%03d", 1:n_samples))
  } else {
    n_samples <- length(sample_names)
  }
  
  template <- data.frame(
    Basename = sample_names,
    Age = rep(NA, n_samples),
    Female = rep(NA, n_samples),
    SpeciesLatinName = rep("Mus musculus", n_samples),
    Tissue = rep(NA, n_samples),
    Group = rep(NA, n_samples),
    Treatment = rep(NA, n_samples),
    stringsAsFactors = FALSE
  )
  
  # Add helpful comments as the first few rows
  comments <- data.frame(
    Basename = c("# REQUIRED: Sample identifiers matching your methylation data column names",
                 "# REQUIRED: Chronological age in years (e.g., 1.5 for 18 months old mouse)",
                 "# OPTIONAL: Sex - 1 for female, 0 for male, NA if unknown",
                 "# OPTIONAL: Species - 'Mus musculus', 'Rattus norvegicus', 'Homo sapiens', etc.",
                 "# OPTIONAL: Tissue type - 'Blood', 'Liver', 'Brain', etc.",
                 "# OPTIONAL: Experimental group identifier",
                 "# OPTIONAL: Treatment condition",
                 "# DELETE THESE COMMENT ROWS BEFORE USING"),
    Age = rep("", 8),
    Female = rep("", 8),
    SpeciesLatinName = rep("", 8),
    Tissue = rep("", 8),
    Group = rep("", 8),
    Treatment = rep("", 8),
    stringsAsFactors = FALSE
  )
  
  # Combine comments and template
  full_template <- rbind(comments, template)
  
  if (!is.null(file_path)) {
    write.csv(full_template, file_path, row.names = FALSE)
    cat("Sample sheet template saved to:", file_path, "\n")
    cat("Please fill in the required information and delete the comment rows.\n")
  }
  
  return(full_template)
}
