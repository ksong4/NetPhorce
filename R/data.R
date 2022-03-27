#' A MaxQuant result file of an experiment with one experimental condition and 5 time points.
#'
#' @description The demo file is a MaxQuant result file for an experiment with one experimental condition, 5 time points, and 4 replicates. This file contains 340 detected peptides and 299 columns. Only a selection of those columns will be used for the actual data analysis.
#'
#' @usage data(oneConditionExample)
#' @docType data
#' @format A data frame with 340 samples and 299 columns, of which,
#' @examples
#' data(oneConditionExample)
"oneConditionExample"

#' A MaxQuant result file of an experiment with two experimental conditions and 5 time points.
#'
#' @description The demo file is a MaxQuant result file for an experiment with two experimental conditions (`loss-of-function line` and `control line`), 5 time points, and 4 replicates. The file contains 790 peptides and 673 columns. Only a selection of those columns will be used for the actual data analysis.
#'
#' @usage data(twoConditionsExample)
#' @docType data
#' @format A data frame with 790 samples and 299 columns, of which,
#' @examples
#' data(twoConditionsExample)
"twoConditionsExample"

#' Pre-loaded Kinase/Phosphatase Table for 25 Species
#'
#' @description The table contains previously published kinase or phosphatase proteins for 25 species. Please use \code{\link{ checkPreloadKinaseTable}} to print a summary for the kinase and phosphatase table.
#'
#' @usage data(kinasesPhosphatases)
#' @docType data
#' @format A data frame with 27,534 kinases and phosphatases records across 25 Species, of which,
#' \describe{
#' \item{Abb}{Three-letter abbreviation for the species}
#' \item{Species}{Full species name}
#' \item{ID}{Kinase/Phosphatase ID}
#' \item{Family}{Kinase or Phosphatase classification}
#' }
#' @examples
#' data(kinasesPhosphatases)
"kinasesPhosphatases"


