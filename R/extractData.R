#' Extract a summary table from processed NetPhorce data
#'
#' @description Obtain the summary table from NetPhorce data after \code{\link{processData}}.
#' @param netPhorceData Processed NetPhorce Object
#' @return A single data.frame
#' @export
#' @examples
#' \dontrun{
#' ## Loading One Condition Data
#' data("oneConditionExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = oneConditionExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = oneConditionExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = oneConditionExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Extract Summary Table
#' summaryTable <- extractSummaryTable(netPhorceData = netPhorceData)
#' }
extractSummaryTable <- function(netPhorceData = netPhorceData){
  if(!is.null(netPhorceData@data.filtered)){
    data.clustered.avg_V2 = netPhorceData@data.filtered %>%
      group_by(set, UniqueID, timepoint, condition) %>% filter(obs_tp >= netPhorceData@Design$minReplication) %>%
      summarize(m = mean(normValue, na.rm=TRUE), .groups = "drop") %>% ungroup() %>%
      complete(UniqueID, condition,fill= list(set= "UniqueSet", m=0),
               timepoint = netPhorceData@Design$uniqueTimes$Origin[1]) %>%
      spread(timepoint, m, fill=0)
    if(is.null(netPhorceData@data.filtered.aov.summary)){
      Output = data.clustered.avg_V2 %>%
        left_join( netPhorceData@data.filtered  %>% mutate(Sample = paste(timepoint, replicate, sep = "_")) %>%
                    dplyr::select(Sample, condition, normValue, UniqueID) %>% spread(Sample, normValue, fill = 0),
                   by = c("UniqueID", "condition")) %>%
        separate(UniqueID, into = c("Part1", "multiplicity"), sep = "_(?=[0-9]+$)") %>%
        separate(Part1, into = c("Protein", "AA"), sep = "_(?=[0-9A-z;]+$)") %>%
        left_join(netPhorceData@Misc$seqWinData, by = c("Protein", "AA")) %>%
        relocate(`Sequence Window`, .after = AA)
    } else {
      Output = data.clustered.avg_V2 %>% left_join(netPhorceData@data.filtered.aov.summary %>%
                                                     dplyr::select(UniqueID, p.value, qvalue)) %>%
        left_join( netPhorceData@data.filtered  %>% mutate(Sample = paste(timepoint, replicate, sep = "_")) %>%
                    dplyr::select(Sample, condition, normValue, UniqueID) %>% spread(Sample, normValue, fill = 0),
                   by = c("UniqueID", "condition")) %>%
        separate(UniqueID, into = c("Part1", "multiplicity"), sep = "_(?=[0-9]+$)") %>%
        separate(Part1, into = c("Protein", "AA"), sep = "_(?=[0-9A-z;]+$)") %>%
        left_join(netPhorceData@Misc$seqWinData, by = c("Protein", "AA"))%>%
        relocate(`Sequence Window`, .after = AA)
    }
    return(Output)
  } else {
    print("NetPhorce data was incomplete. Please process the data again.")
  }
}


#' Extract the network inference results from processed NetPhorce data
#'
#' @description Obtain the network inference from NetPhorce data after \code{\link{networkAnalysis}}.
#' @param netPhorceData Processed NetPhorce Object
#' @return A single data.frame
#' @export
#' @examples
#' \dontrun{
#' ## Loading One Condition Data
#' data("oneConditionExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = oneConditionExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = oneConditionExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = oneConditionExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#'
#' ## Confirm the correct Kinase/Phosphotase information for the data
#' netPhorceData <- validateKinaseTable(netPhorceData = netPhorceData,
#'                                      defaultKinaseTable = TRUE,
#'                                      abbrev = "Ath")
#' ## Evaluate the netPhorce object with regulation thresholds
#' netPhorceData <- regulationCheck(netPhorceData = netPhorceData,
#'                                  upReg = 0.25,
#'                                  downReg = 0.25,
#'                                  absMinThreshold = 0.1,
#'                                  qValueCutOff = 0.05,
#'                                  verbose = TRUE)
#' ## Network Analysis
#' netPhorceData <- networkAnalysis(netPhorceData = netPhorceData,
#'                                requestPlotData = TRUE)
#' ## Extract Network Results
#' networkResults <- extractNetworkResults(netPhorceData = netPhorceData)
#' }
extractNetworkResults <- function(netPhorceData = netPhorceData){
  if(!is.null(netPhorceData@netPhorceResult)){
    return(netPhorceData@netPhorceResult)
  } else {
    print(paste0("The network is missing. Please follow the instruction for",
                 "network analysis and try again after networkAnalysis() step."))
  }
}

#' Load Shiny Plot Object to R
#'
#' @description Load the Shiny Plot Objects to R for more manipulation of colors, fills and layout.
#' @param ShinyRDS RDS object downloaded for a plot in NetPhorce Shiny.
#' @return A NetPhorce Object
#' @export
convertShinyToR <- function(ShinyRDS = ShinyRDS){
  if(is.null(ShinyRDS$PlotType)){
    stop("Please check the RDS data, the current file does not contain any plotting data. ")
  }
  if(ShinyRDS$PlotType %in% (c("DistributionPlot", "Histogram/Boxplot", "PCA_Normalized_Plot", "PCA_Unnormalized_Plot",  "DistributionPlot",
                              "SignificantPeptidesHeatmap", "AbsencePresencePeptidesHeatmap",
                              "IndividualPeptide", "MultiplePeptides_Top","MultiplePeptides_Bottom"))){
    returnData <- new(
      "netPhorce",
      Design = ShinyRDS$Design,
      Misc = ShinyRDS$Misc,
      data.filtered = ShinyRDS$data.filtered,
      data.filtered.aov.summary = ShinyRDS$data.filtered.aov.summary
    )
  } else if(ShinyRDS$PlotType == "RegulationPlot") {
    returnData <- new(
      "netPhorce",
      Design = ShinyRDS$Design,
      Misc = ShinyRDS$Misc,
      data.filtered = ShinyRDS$data.filtered,
      data.filtered.aov.summary = ShinyRDS$data.filtered.aov.summary
    )
    returnData@regulationData = ShinyRDS$regulationData
  } else if(ShinyRDS$PlotType %in% c("VisNetworkPlot_Left", "VisNetworkPlot_Right")) {
    returnData <- new(
      "netPhorce",
      Design = ShinyRDS$Design,
      Misc = ShinyRDS$Misc,
      data.filtered = ShinyRDS$data.filtered,
      data.filtered.aov.summary = ShinyRDS$data.filtered.aov.summary
    )
    returnData@regulationData = ShinyRDS$regulationData
    returnData@networkPlotData = ShinyRDS$networkPlotData
  }
  return(returnData)
}



