#' EnsembleAge: Epigenetic Age Prediction Using Ensemble Clock Methods
#'
#' @description
#' EnsembleAge provides tools for predicting epigenetic age using various clock methods 
#' including Universal Clocks, Elastic Epigenetic clocks, static ensemble methods, and 
#' dynamic ensemble clocks. The package automatically detects your data platform and 
#' prepares the data accordingly for accurate age prediction across various tissues and species.
#'
#' @details
#' The package supports multiple methylation array platforms:
#' \itemize{
#'   \item \strong{Mammal320k}: Mammalian Methylation Array (~320,000 probes)
#'   \item \strong{Mammal40k}: Reduced Mammalian Array (~40,000 probes)  
#'   \item \strong{Human EPIC}: Illumina EPIC Array (~850,000 probes)
#'   \item \strong{Human 450k}: Illumina 450k Array (~450,000 probes)
#' }
#'
#' \strong{Main Functions:}
#' \itemize{
#'   \item \code{\link{predictEnsemble}}: Ultimate one-liner - automatic detection and prediction
#'   \item \code{\link{predict_all_clocks}}: All available clock predictions
#'   \item \code{\link{predict_ensemble_static}}: Mouse-optimized ensemble clocks
#'   \item \code{\link{predict_ensemble_dual_static}}: Human-optimized dual-species clocks
#'   \item \code{\link{predict_ensemble_dynamic}}: 50 specialized mouse clocks
#'   \item \code{\link{detect_platform}}: Automatic platform detection
#' }
#'
#' \strong{Quick Start:}
#' \preformatted{
#' library(EnsembleAge)
#' 
#' # Ultimate one-liner - automatically detects platform and runs predictions
#' results <- predictEnsemble("path/to/methylation_data.RDS", "path/to/sample_sheet.csv")
#' 
#' # View results
#' head(results)
#' }
#'
#' @references
#' Haghani, A., et al. (2025). EnsembleAge: enhancing epigenetic age assessment 
#' with a multi-clock framework. \emph{GeroScience}. 
#' DOI: 10.1007/s11357-025-01808-1
#' 
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/ahaghani/EnsembleAge}
#'   \item \doi{10.1007/s11357-025-01808-1}
#'   \item Report bugs at \url{https://github.com/ahaghani/EnsembleAge/issues}
#' }
#' 
#' @author Amin Haghani \email{dr.a.haghani@@gmail.com}
#' 
#' @keywords package
"_PACKAGE"
