#' Predict Age and Age Acceleration using Ensemble Clocks
#'
#' This is the main function for predicting epigenetic age using various clock methods
#' including Universal Clocks, Elastic Epigenetic clocks, and ensemble methods.
#' 
#' @param dat0sesame Data frame containing methylation data with CpG sites as columns
#'   and samples as rows. Must include a 'CGid' column with CpG identifiers.
#' @param samps Data frame containing sample information with columns:
#'   - Basename: Sample identifiers
#'   - Age: Chronological age (optional, defaults to 0)
#'   - SpeciesLatinName: Species name (optional, defaults to "Mus musculus")
#'   - Female: Sex indicator (optional)
#'   - Tissue: Tissue type (optional)
#' @return Data frame with predicted ages and age acceleration values for various clocks
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr spread
#' @importFrom plyr llply
#' @importFrom stringr str_replace
#' @importFrom data.table rbindlist
#' @examples
#' \dontrun{
#' # Load example methylation data
#' methylation_data <- data.frame(
#'   CGid = c("cg12345", "cg67890", "cg11111"),
#'   sample1 = c(0.5, 0.3, 0.8),
#'   sample2 = c(0.4, 0.6, 0.7)
#' )
#' 
#' # Load sample information
#' sample_info <- data.frame(
#'   Basename = c("sample1", "sample2"),
#'   Age = c(0.5, 1.5),
#'   SpeciesLatinName = c("Mus musculus", "Mus musculus")
#' )
#' 
#' # Predict ages
#' results <- predictAgeAndAgeAcc(methylation_data, sample_info)
#' }
predictAgeAndAgeAcc <- function(dat0sesame, samps) {
  
  # Input validation
  if (missing(dat0sesame) || missing(samps)) {
    stop("Both dat0sesame and samps arguments are required")
  }
  
  if (!is.data.frame(dat0sesame) || !is.data.frame(samps)) {
    stop("Both dat0sesame and samps must be data frames")
  }
  
  if (!"CGid" %in% names(dat0sesame)) {
    stop("dat0sesame must contain a 'CGid' column with CpG identifiers")
  }
  
  if (!"Basename" %in% names(samps)) {
    stop("samps must contain a 'Basename' column with sample identifiers")
  }
  
  # Check for matching samples
  sample_cols <- names(dat0sesame)[!names(dat0sesame) %in% "CGid"]
  matching_samples <- intersect(samps$Basename, sample_cols)
  
  if (length(matching_samples) == 0) {
    stop("No matching samples found between dat0sesame columns and samps$Basename.\n",
         "dat0sesame sample columns: ", paste(sample_cols[1:min(5, length(sample_cols))], collapse = ", "), "\n",
         "samps$Basename values: ", paste(samps$Basename[1:min(5, nrow(samps))], collapse = ", "))
  }
  
  if (length(matching_samples) < nrow(samps)) {
    missing_samples <- setdiff(samps$Basename, sample_cols)
    warning("Some samples in samps$Basename not found in dat0sesame: ",
            paste(missing_samples[1:min(3, length(missing_samples))], collapse = ", "))
  }
  
  # Load species data
  species_file_path <- system.file("data", "anAgeUpdatedCaesarVersion51.csv", 
                                   package = "EnsembleAge")
  if (species_file_path == "") {
    # Fallback to local file if package not installed
    species_file_path <- file.path("data", "anAgeUpdatedCaesarVersion51.csv")
  }
  species <- read.csv(species_file_path)
  
  species$SpeciesLatinName <- sapply(1:nrow(species), function(x) {
    if (is.na(species$SpeciesLatinName[x])) {
      species$SpeciesLatinName[x] <- stringr::str_replace(species$labelsCaesar[x], "_", " ")
    } else {
      species$SpeciesLatinName[x] <- as.character(species$SpeciesLatinName[x])
    }
  })
  
  species <- species %>% 
    mutate(GestationTimeInYears = Gestation.Incubation..days. / 365) %>%
    mutate(logGest = log(Gestation.Incubation..days.), 
           logLifespan = log(maxAgeCaesar), 
           logWeight = log(weightCaesar)) %>%
    dplyr::select(SpeciesLatinName, GestationTimeInYears, averagedMaturity.yrs, 
                  maxAgeCaesar, weightCaesar, logGest, logLifespan, logWeight, 
                  MammalNumberHorvath)
  
  # Remove overlapping columns
  n <- which(names(samps) %in% names(species)[-1])
  if (sum(n) > 0) {
    samps <- samps[, -n]
  }
  
  # Add default columns if missing
  if (!"SpeciesLatinName" %in% names(samps)) {
    samps <- samps %>% mutate(SpeciesLatinName = "Mus musculus")
  }
  
  if (!"Age" %in% names(samps)) {
    samps <- samps %>% mutate(Age = 0)
  }
  
  if (!"Female" %in% names(samps)) {
    samps <- samps %>% mutate(Female = NA)
  }
  
  if (!"Tissue" %in% names(samps)) {
    samps <- samps %>% mutate(Tissue = NA)
  }
  
  samps <- samps %>% 
    mutate(Age = ifelse(is.na(Age), 0, Age)) %>% 
    mutate(SpeciesLatinName = ifelse(is.na(SpeciesLatinName), "Mus musculus", SpeciesLatinName))
  
  # Prepare methylation data
  dat0sesame <- dat0sesame %>% 
    tibble::column_to_rownames("CGid") %>% 
    t(.) %>% 
    as.data.frame(.) %>% 
    mutate('Intercept' = 1)
  
  samps <- samps %>% left_join(species, by = "SpeciesLatinName")
  
  # Load clock coefficients
  clock_file_path <- system.file("data", "Clock_coefficients.RDS", 
                                 package = "EnsembleAge")
  if (clock_file_path == "") {
    # Fallback to local file if package not installed
    clock_file_path <- file.path("data", "Clock_coefficients.RDS")
  }
  epiclocks <- readRDS(clock_file_path)
  
  ageResults <- lapply(1:length(epiclocks), function(j) {
    
    clocks <- epiclocks[[j]]
    
    # a loop for the main category of clocks
    ageAccel <- plyr::llply(1:length(clocks), function(i) {
      # cat(paste0("Predicting age for ", names(epiclocks)[j], " clock ", names(clocks)[i], "\n"))
      clock <- clocks[[i]]
      
      dat1 <- dat0sesame %>% dplyr::select(clock$CGid)
      
      if (grepl("AgeTraf", names(epiclocks)[j])) {
        samp <- samps %>%
          # predicting age
          mutate(epiAge = as.numeric(as.matrix(dat1) %*% clock$Coef)) %>% 
          mutate(AgeAccelation = as.vector(residuals(lm(epiAge ~ Age)))) %>% 
          dplyr::select(Basename, Age, AgeAccelation, epiAge, Female, Tissue) %>%
          # age transformation
          mutate(epiAge = anti.trafo(epiAge)) %>% 
          mutate(AgeAccelation = as.vector(residuals(lm(epiAge ~ Age))))
          
      } else if (grepl("uniClocks", names(epiclocks)[j]) & 
                 grepl("(UniClock3)|(UniBloodClock3)|(UniSkinClock3)", names(clocks)[i])) {
        
        samp <- F3_loglifn(samps) # to compute m estimate for tuning point in the log-linear transformation
        samp <- samp %>%
          # predicting age
          mutate(epiAge = as.numeric(as.matrix(dat1) %*% clock$Coef)) %>% 
          mutate(m1 = a_Logli) %>%
          mutate(epiAge = F2_revtrsf(epiAge, m1) * (averagedMaturity.yrs + GestationTimeInYears) - GestationTimeInYears) %>%
          dplyr::select(-m1, -a1_Logli, -a_Logli, -LogliAge) %>% 
          mutate(AgeAccelation = as.vector(residuals(lm(epiAge ~ Age)))) %>% 
          dplyr::select(Basename, epiAge, AgeAccelation, Age, Female, Tissue)
          
      } else if (grepl("uniClocks", names(epiclocks)[j]) & 
                 grepl("(UniClock2)|(Ake.DuoHumanMouse)|(UniBloodClock2)|(UniSkinClock2)|(Ake.DNAmDuoGrimAge411)", names(clocks)[i])) {
        
        samp <- samps %>%
          # predicting age
          mutate(epiAge = as.numeric(as.matrix(dat1) %*% clock$Coef)) %>% 
          mutate(HighmaxAgeCaesar = maxAgeCaesar * 1.3) %>%
          mutate(HighmaxAgeCaesar = ifelse(SpeciesLatinName %in% c("Homo sapiens", "Mus musculus"), 
                                           maxAgeCaesar, HighmaxAgeCaesar)) %>%
          # age transformation
          mutate(epiAge = F2_antitrans(epiAge, y.maxAge = HighmaxAgeCaesar, 
                                       y.gestation = GestationTimeInYears, const = 1)) %>% 
          mutate(AgeAccelation = as.vector(residuals(lm(epiAge ~ Age)))) %>% 
          dplyr::select(Basename, epiAge, AgeAccelation, Age, Female, Tissue)
          
      } else if (grepl("EnsembleAge.Dynamic", names(epiclocks)[j]) && grepl("AgeTransformed", names(clocks)[i])) {
        samp <- samps %>%
          # predicting age
          mutate(epiAge = as.numeric(as.matrix(dat1) %*% clock$Coef)) %>% 
          mutate(AgeAccelation = as.vector(residuals(lm(epiAge ~ Age)))) %>% 
          dplyr::select(Basename, Age, AgeAccelation, epiAge, Female, Tissue) %>%
          # age transformation for EnsembleAge.Dynamic clocks with AgeTransformed
          mutate(epiAge = anti.trafo(epiAge)) %>% 
          mutate(AgeAccelation = as.vector(residuals(lm(epiAge ~ Age))))
      } else if (grepl("EnsembleDualAge.Static", names(epiclocks)[j])) {
        samp <- samps %>%
          # predicting age for EnsembleDualAge.Static (uses RelativeAge transformation)
          mutate(epiAge = as.numeric(as.matrix(dat1) %*% clock$Coef)) %>% 
          # Convert RelativeAge back to chronological age: RelativeAge * maxAgeCaesar
          mutate(epiAge = epiAge * maxAgeCaesar) %>%
          mutate(AgeAccelation = as.vector(residuals(lm(epiAge ~ Age)))) %>% 
          dplyr::select(Basename, Age, AgeAccelation, epiAge, Female, Tissue)
      } else {
        samp <- samps %>%
          # predicting age
          mutate(epiAge = as.numeric(as.matrix(dat1) %*% clock$Coef)) %>% 
          mutate(AgeAccelation = as.vector(residuals(lm(epiAge ~ Age)))) %>% 
          dplyr::select(Basename, Age, AgeAccelation, epiAge, Female, Tissue)
      }
      
      return(samp)
    })
    
    names(ageAccel) <- names(clocks)
    ageAccel <- rbindlist(ageAccel, idcol = "epiClock")
  })
  
  names(ageResults) <- names(epiclocks)
  
  ageResults <- rbindlist(ageResults, idcol = "clockFamily", fill = T) %>% 
    # Clean up EnsembleAge.Dynamic clock names by converting X123.Name to Clock123.Name
    mutate(epiClock = ifelse(grepl("^X[0-9]+\\.", epiClock), 
                            gsub("^X([0-9]+)\\.", "Clock\\1.", epiClock), 
                            epiClock)) %>%
    mutate(clockFamily = ifelse(grepl("(UniClock2)|(UniBloodClock2)|(UniSkinClock2)", epiClock), "UniClock2", 
                                ifelse(grepl("(UniClock3)|(UniBloodClock3)|(UniSkinClock3)", epiClock), "UniClock3", clockFamily))) %>%
    mutate(epiClock = ifelse(grepl("Skin", epiClock), "Skin", 
                             ifelse(grepl("Blood", epiClock), "Blood", 
                                    ifelse(grepl("(UniClock)|(EnsemblDualAge)", epiClock), 
                                           "panTissue", epiClock))))
  
  # convert to wide format
  targetClocks <- c("LifespanUberClock", "DNAmAgeElasticFinal", 
                    "DNAmAgeInterventionFinal", "DNAmAgeDevelopmentFinal", 
                    "UniClock2", "UniClock3", 
                    "Ensemble.Static", "Ensemble.Static.Top", 
                    "EnsembleAge.Dynamic", "EnsembleDualAge.Static")
  newNames <- c("LifespanUberClock", "DNAmAgeElasticFinal", 
                "DNAmAgeInterventionFinal", "DNAmAgeDevelopmentFinal",
                "UniClock2", "UniClock3", 
                "Ensemble.Static", "Ensemble.Static.Top",
                "EnsembleAge.Dynamic", "EnsembleDualAge.Static")
  
  dat1 <- ageResults %>% 
    mutate(clockFamily = factor(clockFamily, levels = targetClocks, labels = newNames)) %>% 
    mutate(epiClock = paste(clockFamily, epiClock, "clock", "epiAge", sep = ".")) %>% 
    dplyr::select(Basename, epiClock, epiAge) %>% 
    spread(key = "epiClock", value = "epiAge") 
  
  dat2 <- ageResults %>% 
    mutate(clockFamily = factor(clockFamily, levels = targetClocks, labels = newNames)) %>% 
    mutate(epiClock = paste(clockFamily, epiClock, "clock", "AgeAcceleration", sep = ".")) %>% 
    dplyr::select(Basename, epiClock, AgeAccelation) %>% 
    spread(key = "epiClock", value = "AgeAccelation") %>% 
    left_join(dat1, by = "Basename")
  
  # For compatibility with user functions, also return the long format
  attr(dat2, "long_format") <- ageResults
  
  return(dat2)
}
