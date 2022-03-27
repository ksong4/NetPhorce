#' Network Analysis
#' @description The phosphorylation signaling network shows the inferred causal relations between the kinases/phosphatases and downstream peptides. Network inference is based on dynamic Bayesian principles. For more than one condition, the common phosphopeptides between the conditions are used as input and a network is generated for each condition separately. As such, the network results might differ slightly when running each condition separately through NetPhorce compared to running the conditions together.
#' @param netPhorceData (Required). Processed netPhorceData data
#' @param requestPlotData (Required). If TRUE, data required to visualize the network in R via \code{\link{visNetwork}} package will be saved.
#' @return The netPhorce object with the following slots:
#' \describe{
#'   \item{Design}{Contains data design information and filtering parameters}
#'   \item{data.filtered}{Contains all the data points in a `Long` format that passed filtering criteria.}
#'   \item{data.filtered.aov.summary}{Contains all the data points in a `Long` format with anova results.}
#'   \item{Misc}{Contains accessory  data including default plotting colors and FASTA Keys, if present.}
#'   \item{regulationData}{Contains regulation data calculated through \code{\link{regulationCheck}} function}
#'   \item{networkPhorceResutls}{Contains the network results}
#'   \item{networkPlotData}{Contains data required to generate the network via \code{\link{visNetwork}} package}
#' }
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
#' ## Validating the Kinase/Phosphatase Information
#' netPhorceData <- validateKinaseTable(netPhorceData = netPhorceData,
#'                                      defaultKinaseTable = TRUE,
#'                                      abbrev = "Ath")
#' ## Regulation Validation based on user inputs
#' netPhorceData <- regulationCheck(netPhorceData = netPhorceData,
#'                                  upReg = 0.25,
#'                                  downReg = 0.25,
#'                                  absMinThreshold = 0.1,
#'                                  qValueCutOff = 0.05,
#'                                  verbose = TRUE)
#' ## Network Analysis
#' networkData <- networkAnalysis(netPhorceData = netPhorceData,
#'                                requestPlotData = TRUE)
#' }
networkAnalysis <- function(netPhorceData = netPhorceData,
                            requestPlotData = TRUE){
  ESS = 1e-15 ### Dirichlet distr param.
  data.significant.proteins <- netPhorceData@regulationData
  n_levels = data.significant.proteins %>% filter(!is.na(discreteState)) %>% distinct(discreteState) %>% nrow ### determined by numberOfDiscreteStates
  conditions = data.significant.proteins %>% filter(!is.na(experiment)) %>% distinct(experiment) %>% nrow ### determined by numberOfDiscreteStates
  maxConditionMatching = conditions-1
  percentChangeThreshold = 0.5
  # returnList <- list()

  pb = txtProgressBar(min = 0, max = 9,
                      initial = 0, style = 3)

  data_V4 <-
    netPhorceData@data.filtered %>%
    mutate(timepoint = factor(timepoint, levels = c(as.vector(netPhorceData@Design$uniqueTimes$Origin)))) %>%
    arrange(timepoint) %>%
    dplyr::select(timepoint, replicate) %>%
    group_by(timepoint) %>%
    distinct() %>%
    summarise(ReplicateCount = n())

  # print("1")
  setTxtProgressBar(pb, 1)
  ################  . . 3.3.1 Identify all possible target-regulator combinations ############
  data.scoring = data.significant.proteins %>% group_by(condition) %>% nest()
  # print(head(data.scoring$data))
  # print("1.5")

  data.scoring = data.scoring %>%  ### expand the possible edges from kinases to genes, then filter to keep only data where changes correspond.
    mutate(potentialEdges = map(data, ~ { ### Reg -> target are filtered by sigSignChange, not by change in state
      allEdges = .x %>% dplyr::select(UniqueID, isKinase, experiment, sigChangeSign) %>% filter(!is.na(sigChangeSign)) %>%
        tidyr::expand(nesting(regulator = UniqueID, regulatorCondition = experiment, regulatorChange = sigChangeSign, regulatorisKinase = isKinase),
                      nesting(target = UniqueID, targetCondition = experiment, targetChange = sigChangeSign)) %>%
        filter(regulatorisKinase) %>% dplyr::select(-regulatorisKinase) %>%
        filter(regulator != target) %>%
        filter(regulatorChange != 0 & targetChange != 0) %>%
        filter(!is.na(levels(targetCondition)[as.numeric(targetCondition)-1])) %>% ### The regulators must have the proper timepoints and proper % matching changes
        filter(regulatorCondition == targetCondition | regulatorCondition == levels(targetCondition)[as.numeric(targetCondition)-1]) ### Allow co-regulators

      allEdges.filterd = allEdges %>% group_by(regulator, target) %>%
        add_count(coReg = regulatorCondition == targetCondition, wt = (abs(regulatorChange) == abs(targetChange))/maxConditionMatching, name = "percentChange") %>%
        ungroup() %>% filter(percentChange > percentChangeThreshold)
      return(allEdges.filterd)
    }))
  # print(head(data.scoring))

  # print("2")
  setTxtProgressBar(pb, 2)
  message(paste0("\nUndergoing a time consuming step, please wait... "))

  data.scoring.searches = data.scoring %>% dplyr::select(condition, potentialEdges) %>% unnest(cols = potentialEdges) %>%
    mutate(regulator = paste(regulator, coReg, sep = "**"), .keep = "unused") %>% group_by(condition, target) %>% nest() %>%
    mutate(neighborCombinations = map(data, possibly( ~ { ### possibly for error catching.
      allSearchCombinations = bind_rows( ### Generate combinations of regulators whithin each searchID to be scored.
        .x %>% distinct(regulator) %>% pull(regulator) %>% as.list %>% possibly(~combn(.x, m=1), tibble())() %>% as.data.frame %>% gather(searchID,regulator) %>% mutate(searchSizeK = 1),
        .x %>% distinct(regulator) %>% pull(regulator) %>% as.list %>% possibly(~combn(.x, m=2), tibble())() %>% as.data.frame %>% gather(searchID,regulator) %>% mutate(searchSizeK = 2),
        .x %>% distinct(regulator) %>% pull(regulator) %>% as.list %>% possibly(~combn(.x, m=3), tibble())() %>% as.data.frame %>% gather(searchID,regulator) %>% mutate(searchSizeK = 3)
      ) %>% as_tibble %>% unnest(cols = c(regulator)) %>% tidyr::separate(regulator, c("regulator", "coReg"), sep = "\\*\\*", remove = T)

    },tibble()))) %>% ungroup()
  # print(head(data.scoring.searches))

  # print("3")
  setTxtProgressBar(pb, 3)
  ######### . . 3.3.2 Scoring of the identified target - (co)regulator combinations ###############
  coRegulater = data.scoring.searches %>% dplyr::select(condition, target, neighborCombinations) %>% unnest(neighborCombinations) %>%
    dplyr::select(condition, target, regulator, coReg) %>% distinct()

  allEdgesCompareState.Nj = data.significant.proteins %>% group_by(condition) %>% nest() %>%
    mutate(data = map(data, ~ {.x %>% dplyr::select(UniqueID, isKinase, experiment, discreteState) %>%
        filter(!is.na(discreteState)) %>%
        tidyr::expand(nesting(regulator = UniqueID, regulatorCondition = experiment, regulatorState = discreteState, regulatorisKinase = isKinase),
                      nesting(target = UniqueID, targetCondition = experiment, targetState = discreteState)) %>% filter(regulatorisKinase) })) %>%
    unnest(cols = c(data)) %>% right_join(coRegulater, by = c("condition", "regulator", "target")) %>% group_by(condition) %>% nest() %>%
    mutate(data = map(data, ~ {.x %>% filter((as.numeric(regulatorCondition) == as.numeric(targetCondition) & coReg == TRUE) | (as.numeric(regulatorCondition) == as.numeric(targetCondition)-1) & coReg == FALSE) %>% filter(regulator != target)})) %>% unnest(cols = c(data))

  # print(head(allEdgesCompareState.Nj))
  # print("4")
  setTxtProgressBar(pb, 4)
  allEdgesCompareState.Nj.full = data.scoring.searches %>%
    dplyr::select(condition, target, neighborCombinations) %>%
    unnest(neighborCombinations) %>% left_join(allEdgesCompareState.Nj, by = c("condition", "target", "regulator", "coReg"))

  data.scoring.genist.1reg = allEdgesCompareState.Nj.full %>% filter(searchSizeK == 1) %>% dplyr::select(-regulatorisKinase) %>%
    group_by(condition, target, searchID) %>%
    mutate(Nj = case_when(
      regulatorState == -1 & targetState == -1 ~ "A",
      regulatorState == -1 & targetState == 1 ~ "B",
      regulatorState == -1 & targetState == 5 ~ "C",
      regulatorState == 1 & targetState == -1 ~ "D",
      regulatorState == 1 & targetState == 1 ~ "E",
      regulatorState == 1 & targetState == 5 ~ "F",
      regulatorState == 5 & targetState == 1 ~ "G",
      regulatorState == 5 & targetState == -1 ~ "H",
      regulatorState == 5 & targetState == 5 ~ "I")) %>%
    mutate(qi = case_when(
      regulatorState == -1 ~ "A",
      regulatorState == 1 ~ "B",
      regulatorState == 5 ~ "C")) %>%
    add_count(qi) %>% dplyr::count(condition, target, coReg, searchID, regulator, searchSizeK, Nj, n, name = "Nj.count") %>%
    mutate(BDei =
             log(gamma(Nj.count+ESS/(n_levels*n_levels^searchSizeK))/gamma(ESS/(n_levels*n_levels^searchSizeK)))) %>% ungroup() %>%
    group_by(condition, target, regulator, searchID, searchSizeK, coReg, n) %>% summarize(BDei = sum(BDei)) %>%
    mutate_at(vars(BDei), ~ BDei + log(gamma(ESS/n_levels^searchSizeK)/gamma(n+ESS/n_levels^searchSizeK))) %>%
    group_by(condition, target, regulator, searchID, searchSizeK, coReg) %>% summarize(BDei = sum(BDei)) %>%  ungroup() %>%
    group_by(condition, target) %>% filter(BDei == max(BDei))

  tp0 = as.character(data_V4$timepoint)[1] #QUESTION KUNCHENG; IS DATA_V4 ALWAYS ordered chronologically??

  # print("5")
  setTxtProgressBar(pb, 5)
  # message(paste0("\nUndergoing a time consuming step, please wait... "))

  data.scoring.genist.2reg = allEdgesCompareState.Nj.full %>% filter(searchSizeK == 2) %>%
    dplyr::select(-regulatorisKinase) %>%
    group_by(condition, target, searchID) %>%
    filter(regulatorCondition != tp0 | targetCondition != tp0 | all(coReg == TRUE)) %>%
    ungroup() %>% group_by(condition, target, searchID, targetCondition) %>% mutate(Reg = row_number()) %>%
    pivot_wider(names_from = Reg, values_from= c(regulator, regulatorState, coReg, regulatorCondition)) %>%
    ungroup() %>% group_by(condition, target, searchID) %>%
    mutate(Nj = case_when(
      regulatorState_1 == -1 & regulatorState_2 == -1 & targetState == -1 ~ "A",
      regulatorState_1 == -1 & regulatorState_2 == -1 & targetState == 1 ~ "B",
      regulatorState_1 == -1 & regulatorState_2 == -1 & targetState == 5 ~ "C",
      regulatorState_1 == -1 & regulatorState_2 == 1 & targetState == 1 ~ "D",
      regulatorState_1 == -1 & regulatorState_2 == 1 & targetState == -1 ~ "E",
      regulatorState_1 == -1 & regulatorState_2 == 1 & targetState == 5 ~ "F",
      regulatorState_1 == -1 & regulatorState_2 == 5 & targetState == -1 ~ "G",
      regulatorState_1 == -1 & regulatorState_2 == 5 & targetState == 1 ~ "H",
      regulatorState_1 == -1 & regulatorState_2 == 5 & targetState == 5 ~ "I",
      regulatorState_1 == 1 & regulatorState_2 == -1 & targetState == -1 ~ "J",
      regulatorState_1 == 1 & regulatorState_2 == -1 & targetState == 1 ~ "K",
      regulatorState_1 == 1 & regulatorState_2 == -1 & targetState == 5 ~ "L",
      regulatorState_1 == 1 & regulatorState_2 == 1 & targetState == 1 ~ "M",
      regulatorState_1 == 1 & regulatorState_2 == 1 & targetState == -1 ~ "N",
      regulatorState_1 == 1 & regulatorState_2 == 1 & targetState == 5 ~ "O",
      regulatorState_1 == 1 & regulatorState_2 == 5 & targetState == -1 ~ "P",
      regulatorState_1 == 1 & regulatorState_2 == 5 & targetState == 1 ~ "Q",
      regulatorState_1 == 1 & regulatorState_2 == 5 & targetState == 5 ~ "R",
      regulatorState_1 == 5 & regulatorState_2 == -1 & targetState == -1 ~ "S",
      regulatorState_1 == 5 & regulatorState_2 == -1 & targetState == 1 ~ "T",
      regulatorState_1 == 5 & regulatorState_2 == -1 & targetState == 5 ~ "U",
      regulatorState_1 == 5 & regulatorState_2 == 1 & targetState == 1 ~ "V",
      regulatorState_1 == 5 & regulatorState_2 == 1 & targetState == -1 ~ "W",
      regulatorState_1 == 5 & regulatorState_2 == 1 & targetState == 5 ~ "X",
      regulatorState_1 == 5 & regulatorState_2 == 5 & targetState == -1 ~ "Y",
      regulatorState_1 == 5 & regulatorState_2 == 5 & targetState == 1 ~ "Z",
      regulatorState_1 == 5 & regulatorState_2 == 5 & targetState == 5 ~ "AA")) %>%
    mutate(qi = case_when(
      regulatorState_1 == -1 & regulatorState_2 == -1 ~ "A",
      regulatorState_1 == -1 & regulatorState_2 == 1 ~ "B",
      regulatorState_1 == -1 & regulatorState_2 == 5 ~ "C",
      regulatorState_1 == 1 & regulatorState_2 == -1 ~ "D",
      regulatorState_1 == 1 & regulatorState_2 == 1 ~ "E",
      regulatorState_1 == 1 & regulatorState_2 == 5 ~ "F",
      regulatorState_1 == 5 & regulatorState_2 == -1 ~ "G",
      regulatorState_1 == 5 & regulatorState_2 == 1 ~ "H",
      regulatorState_1 == 5 & regulatorState_2 == 5 ~ "I")) %>% add_count(qi) %>%
    dplyr::count(condition, target, coReg_1, coReg_2, searchID, regulator_1, regulator_2, searchSizeK, n, Nj, name = "Nj.count") %>%
    mutate(BDei =
             log(gamma(Nj.count+ESS/(n_levels*n_levels^searchSizeK))/gamma(ESS/(n_levels*n_levels^searchSizeK)))) %>% ungroup() %>%
    group_by(condition, target, regulator_1, regulator_2, coReg_1, coReg_2, searchID, searchSizeK, n) %>%
    summarize(BDei = sum(BDei)) %>%
    mutate_at(vars(BDei), ~ BDei + log(gamma(ESS/n_levels^searchSizeK)/gamma(n+ESS/n_levels^searchSizeK))) %>%
    group_by(condition, target, regulator_1, regulator_2, coReg_1, coReg_2, searchID, searchSizeK) %>%
    summarize(BDei = sum(BDei)) %>% ungroup() %>% group_by(condition, target) %>% filter(BDei == max(BDei))

  # print("6")
  setTxtProgressBar(pb, 6)
  ##we keep only the max scores across entire dataset instead of grouping per target##
  data.scoring.genist = bind_rows(data.scoring.genist.1reg %>% mutate(set = "1"), data.scoring.genist.2reg %>% pivot_longer(cols = c(ends_with("_1"), ends_with("_2")), names_to = c(".value", "set"), names_pattern = "(.+)_(.+)")) %>%
    distinct(condition, target, regulator, coReg, BDei) %>% group_by(condition, target) %>%
    filter(BDei == max(BDei)) %>% mutate(BDei = exp(BDei), .keep = "unused")

  # print("7")
  setTxtProgressBar(pb, 7)
  ### Sign identification
  data.scoring.sign = data.scoring %>% ### expand with all identified SigChangeSign, count the changes in the same and opposite direction
    dplyr::select(potentialEdges, condition) %>% unnest(c(condition, potentialEdges)) %>% dplyr::select(-percentChange) %>%
    mutate(coReg = as.character(coReg), .keep = "unused") %>%
    right_join(data.scoring.genist, by = c("regulator", "target", "coReg", "condition")) %>%
    group_by(condition, target, regulator, coReg) %>% dplyr::count(regulatorChange == targetChange) %>%
    pivot_wider(names_from = "regulatorChange == targetChange", values_from = n, values_fill = 0) %>%
    mutate(regulation = case_when(
      `FALSE` > `TRUE` ~ "dephosphorylation",
      `TRUE` > `FALSE` ~ "phosphorylation",
      TRUE ~ "undetermined"
    ))

  # print("8")
  setTxtProgressBar(pb, 8)
  network.final = data.scoring.sign %>% dplyr::select(-`FALSE`, -`TRUE`) %>% dplyr::rename(timelapse = coReg) %>%
    mutate_at(vars(timelapse), ~ case_when(timelapse == TRUE ~ "tl.1", timelapse == FALSE ~ "tl.0")) %>%
    left_join(data.scoring.genist %>% dplyr::select(-coReg), by = c("condition", "target", "regulator")) %>% distinct()
  setTxtProgressBar(pb, 9)

  netPhorceData@netPhorceResult <- network.final
  if(requestPlotData == TRUE){
    plotData = list()
    plotData[["data.scoring.sign"]] <- data.scoring.sign
    plotData[["data.scoring.genist"]] <- data.scoring.genist
    netPhorceData@networkPlotData <- plotData
  }
  return(netPhorceData)
}
