## Use setClassUnion to define instance where it can be either data.frame or NULL
#' @keywords internal
setClassUnion("data.frameOrNull", c("data.frame", "NULL"))

#' The main class for NetPhorce data
#'
#' Contains the outputs from the phosphoproteome data analysis, statistical analysis, and network inference, and the necessary data for plotting and data extraction.
#' \describe{
#'   \item{Design}{Contains data design information and filtering parameters}
#'   \item{data.filtered}{Contains the data points that passed filtering criteria in a long format.}
#'   \item{data.filtered.aov.summary}{Contains the  anova results in a long format.}
#'   \item{regulationData}{Contains data associated with regulation check. }
#'   \item{netPhorceResults}{Contains the network inference data. }
#'   \item{networkPlotData}{Contains the data needed for network plotting using \code{\link[visNetwork:visNetwork]{visNetwork::visNetwork()}} from the visNetwork package.}
#'   \item{Misc}{Contains accessory data including default plotting colors and FASTA Keys, if present.}
#' }
#'
#' @name netPhorce-class
#' @rdname netPhorce-class
#' @exportClass netPhorce
#' @keywords internal
setClass(Class="netPhorce",
         representation = representation(
           Design = "list",
           Misc = "list",
           data.filtered = "data.frame",
           data.filtered.aov.summary = "data.frame",
           regulationData = "data.frameOrNull",
           netPhorceResult = "data.frameOrNull",
           networkPlotData = "list"
         )
)

#' method extensions to show for phyloseq objects.
#'
#' See the general documentation of \code{\link[methods]{show}} method for
#' expected behavior.
#'
#' @seealso \code{\link[methods]{show}}
#'
#' @inheritParams methods::show
#' @rdname show-methods
#' @keywords internal
#' @examples
#' # data(netPhorceData)
#' # show(netPhorceData)
#' # netPhorceData
setMethod("show", "netPhorce", function(object){
  cat("netPhorce object", fill=TRUE)

  ## MaxQuant Data Description
  cat("--- MaxQuant Data Description ---\n")
  cat("Number of Input Peptides: \t", object@Misc$dataStats$nInputSamples, " \n")
  cat("Number of Conditions Detected: \t", object@Misc$dataStats$nConditions, "\n")
  cat("\t", paste(sort(unique(object@Design$ConDesign$condition)), collapse = "; "), " \n")
  cat("Number of Time Points Detected: ", object@Misc$dataStats$nTimePoints, " \n")
  cat("\t", paste(sort(unique(object@Design$ConDesign$timepoint)), collapse = "; "), " \n")
  cat("Number of Replicates Detected: \t", object@Misc$dataStats$nReplicates, " \n")
  cat("\t", paste(sort(unique(object@Design$ConDesign$replicate)), collapse = "; "), " \n")

  ## Filtering Criteria
  cat("--- Filtering Criteria ---\n")
  cat("Minimum Replicates Required: \t\t", object@Design$minReplication, " \n")
  cat("Minimum Local Probability Required: \t", object@Design$minLocalProb, " \n")
  # cat("Q-Value Cutoff: ", object@Misc$dataStats$qVal, " \n")


  ## Statistics
  cat("--- Statistics ---\n")
  cat("Number of Presence/Absence Peptides Found: \t\t ", object@Misc$dataStats$nPresenceAbsence, "\n")
  cat("Number of Significant Peptides Found with q-value = 0.1: ", object@Misc$dataStats$nStatsSetqCut, "\n")

  ## Kinase/Phosphotase Table
  cat("--- Kinsae/Phosphatase Table ---\n")
  if(object@Misc$kinaseCheck$KinaseCheck == TRUE){
    print(object@Misc$kinsaseCheck$KinaseSummary)
  } else {
    cat("No Kinases and Phosphotase identified. ")
  }

  ## Regulation
  if(!is.na(object@Misc$dataStats$upReg)){
    cat("--- Regulation ---\n")
    cat("Up Regulation Fold Change Cutoff: \t", object@Misc$dataStats$upReg, " \n")
    cat("Down Regulation Fold Change Cutoff: \t", object@Misc$dataStats$downReg, " \n")
    cat("Bottom Percentage Fold Change Cutoff: \t", object@Misc$dataStats$absMinThreshold, " \n")
    cat("Q-Value Cutoff:\t\t\t\t", object@Misc$dataStats$qValueCutOff, " \n")
  }

  ## Network
  if(!is.null(object@netPhorceResult)){
    cat("--- Network Statistics ---\n")
    cat("Conditions \t Number of Connections\n")
    for(con in unique(object@netPhorceResult$condition)){
      cat("  ", con, "\t",
          length(unique(subset(object@netPhorceResult, condition == con))),
          "\n"
      )
    }
  }
})
