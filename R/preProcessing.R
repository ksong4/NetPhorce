#' Identify the Columns needed for Data Processing
#'
#' @description The function selects the correct columns. Please note that the column name is case-sensitive.
#' @param rawMaxQuant (Required). Loaded raw MaxQuant data
#' @param positionCol (Required). List the column name with the phosphorylation position.
#' @param reverseCol (Required). List the column with the reverse information.
#' @param localizationProbCol (Required). List the column name with the localization probability.
#' @param potentialContaminationCol (Required). List the column name with the contamination information.
#' @param aminoAcidCol (Required). List the column name with the amino acid.
#' @param uniqueIDCol (Required). List the column name with the protein  for the peptides.
#' @param seqWindowIDCol (Required). List the column name with the sequence window.
#' @param fastaIDCol (Optional). List the column name with the FASTA header.
#' @return A list of data.frames
#' @export
#' @examples
#' \dontrun{
#' ## Loading One Condition Data
#' data("oneConditionExample")
#' ## Identify the Columns needed for Data Processing
#' identifiedCols <- confirmColumnNames(rawMaxQuant = oneConditionExample,
#'                                 	    positionCol = "Position",
#'                                 	    reverseCol = "Reverse",
#'                                      localizationProbCol = "Localization prob",
#'                                      potentialContaminationCol = "Potential contaminant",
#'                                      aminoAcidCol = "Amino acid",
#'                                 	    uniqueIDCol = "Protein",
#'                                      seqWindowIDCol = "Sequence window",
#'                     	                fastaIDCol = "Fasta headers")
#' }
confirmColumnNames = function(rawMaxQuant = rawMaxQuant,
                              positionCol = positionCol,
                              reverseCol = reverseCol,
                              localizationProbCol = localizationProbCol,
                              potentialContaminationCol = potentialContaminationCol,
                              aminoAcidCol = aminoAcidCol,
                              uniqueIDCol = uniqueID,
                              seqWindowIDCol = seqWindowID,
                              fastaIDCol = NULL){
  ## Text Progress bar
  # pb = txtProgressBar(min = 0, max = 12, initial = 0, style = 3)
  if(is.data.frame(rawMaxQuant)){
    print("MaxQuant Data Matrix Statistics: ")
    print(paste0("Number of Rows (Peptides):     ", nrow(rawMaxQuant)))
    print(paste0("Number of Columns (Variables): ", ncol(rawMaxQuant)))
  }



  ## Name if any of the input column name are not presenting
  if (!(positionCol %in% colnames(rawMaxQuant))) {
    stop("`Position` Column does not present in the input MaxQuant Data. The column name is case-sensitive. ", call. = FALSE)
  }
  if (!(reverseCol %in% colnames(rawMaxQuant))) {
    stop("`Reverse` Column does not present in the input MaxQuant Data. The column name is case-sensitive. ", call. = FALSE)
  }
  if (!(localizationProbCol %in% colnames(rawMaxQuant))) {
    stop("`Localization Probability` Column does not present in the input MaxQuant Data. The column name is case-sensitive. ", call. = FALSE)
  }
  if (!(potentialContaminationCol %in% colnames(rawMaxQuant))) {
    stop("`Potential Contamination` Column does not present in the input MaxQuant Data. The column name is case-sensitive. ", call. = FALSE)
  }
  if (!(aminoAcidCol %in% colnames(rawMaxQuant))) {
    stop("`Amino Acid` Column does not present in the input MaxQuant Data. The column name is case-sensitive. ", call. = FALSE)
  }
  if (!(uniqueIDCol %in% colnames(rawMaxQuant))) {
    stop("`Unique ID` Column does not present in the input MaxQuant Data. The column name is case-sensitive. ", call. = FALSE)
  }
  if (!(seqWindowIDCol %in% colnames(rawMaxQuant))) {
    stop("`Sequence Window` Column does not present in the input MaxQuant Data. The column name is case-sensitive. ", call. = FALSE)
  }
  if(!is.null(fastaIDCol)){
    if (!(fastaIDCol %in% colnames(rawMaxQuant))) {
      stop("`FASTA Header` Column does not present in the input MaxQuant Data. The column name is case-sensitive. ", call. = FALSE)
    }
  }

  print("Required Columns are all Found")

  returnList <- list(
    positionCol = positionCol,
    reverseCol = reverseCol,
    localizationProbCol = localizationProbCol,
    potentialContaminationCol = potentialContaminationCol,
    aminoAcidCol = aminoAcidCol,
    uniqueIDCol = uniqueIDCol,
    seqWindowIDCol = seqWindowIDCol,
    fastaIDCol = fastaIDCol
  )
  return(returnList)

}

#' Identify the Intensity Columns
#'
#' @description The function identifies the intensity columns. To correctly identify the columns and for downstream analysis the intensity column names should consist of a string of characters that includes the Condition, Time Point, and Replication or Independent experiment. These three variables need to be separated by an underscore, "_". The function takes as input CON, TIME, REP, separated by an underscore and in the order that these variables are found in the column names. For example, "CON_TIME_REP" will parse the intensity column named "Condition1_0minutes_1".
#' @param rawMaxQuant (Required). Loaded raw MaxQuant data
#' @param intensityPattern (Required). A string of CON, TIME, and REP separated by an underscore, "_". The condition term cannot start with a numeric value and the only acceptable symbols allowed in these three terms are Full-Stop/Period (.) or minus-dash (-). E.g. "CON_TIME_REP".
#' @param filterCon (Optional). A list of conditions that you want to process. Only conditions from this list will be included in the downstream analyses.. E.g., c("Con1", "Con2)
#' @param filterTime (Optional). A list of time points that you want to process. Only time points within this list will be included in the downstream analyses. E.g., c("0min", "5min", "15min")
#' @param filterRep (Optional). A list of replications that you want to process. Only replications from this list will be included in the downstream analyses. E.g., c("R1", "R3", "R4")
#' @param verbose (Optional). If TRUE, print a summary table of the intensity columns statistics for verification.
#' @return list of data.frames
#' @export
#' @examples
#' \dontrun{
#' ## Loading One Condition Data
#' data("oneConditionExample")
#' ## Identify the Intensity Columns and Condition, Time Point and Replication
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = oneConditionExample,
#'                                      	intensityPattern = "con_time_rep",
#'                                      	verbose = TRUE)
#' }
confirmIntensityColumns = function(rawMaxQuant = rawMaxQuant,
                                   intensityPattern = intensityPattern,
                                   filterCon = NULL,
                                   filterTime = NULL,
                                   filterRep = NULL,
                                   verbose = TRUE){
  ## Intensity Pattern extractor
  intensityCol_Input <- intensityPattern
  curIntensityPattern <- strsplit(intensityCol_Input, "_")[[1]]
  curIntensityPattern <- toupper(curIntensityPattern)
  conLoc <- match("CON", curIntensityPattern)
  timeLoc <- match("TIME", curIntensityPattern)
  repLoc <- match("REP", curIntensityPattern)

  if (length(curIntensityPattern) != 3) {
    stop("Intensity Columns must have 3 strings, con, time or rep", call. = FALSE)
  } else {
    if (!("CON" %in% curIntensityPattern)) {
      stop("`CON` does not present in the intensity pattern submitted. ", call. = FALSE)
    }
    if (!("TIME" %in% curIntensityPattern)) {
      stop("`TIME` does not present in the intensity pattern submitted. ", call. = FALSE)
    }
    if (!("REP" %in% curIntensityPattern)) {
      stop("`REP` does not present in the intensity pattern submitted. ", call. = FALSE)
    }
  }


  intensityColNameOrder <- data.frame(input = c("CON", "TIME", "REP"),
                                      order = c(conLoc, timeLoc, repLoc),
                                      truth = c("condition","timepoint","replicate"))
  intensityColNameOrder <- intensityColNameOrder %>% arrange(order)
  intensityCol <- grep("Intensity.+([[:alnum:]])+_([[:alnum:]])+_([[:alnum:]])+___",
                       colnames(rawMaxQuant), value = TRUE)

  colTable <- data.frame(Full = intensityCol)
  colTable$Info1 <- sapply(strsplit(colTable$Full, "_|\\s"), "[", 2)
  colTable$Info2 <- sapply(strsplit(colTable$Full, "_|\\s"), "[", 3)
  colTable$Info3 <- sapply(strsplit(colTable$Full, "_|\\s"), "[", 4)
  nInfo1 <- length(unique(colTable$Info1))
  nInfo2 <- length(unique(colTable$Info2))
  nInfo3 <- length(unique(colTable$Info3))
  colTable <- colTable[, 2:4]
  tempSplit <- toupper(strsplit(intensityCol_Input, "_")[[1]])
  colnames(colTable) <- c(tempSplit[1], tempSplit[2], tempSplit[3])
  colTable <- colTable %>%
    distinct() %>%
    group_by(.data$CON, .data$TIME) %>%
    summarise(REPs = paste0(.data$REP, collapse = ";")) %>%
    pivot_wider(names_from = CON, values_from = REPs)
  colnames(colTable)[1] <- "Time"

  ## Processing Filtering for filtering
  ### Filtering Condition
  if(!(is.null(filterCon)) |!(is.null(filterTime)) | !(is.null(filterRep))){
    applyFilter = TRUE
  } else {
    applyFilter = FALSE
  }

  if(verbose){
    cat("First Found Column Name is:\n")
    cat(as.character(intensityCol[1]))
    cat("\n")
    cat("---------------------------------------------------------------------------------------\n")
    cat("Summary below is shown as unique Replications per Condition for each unique Time Point:\n")
    print(as.data.frame(colTable), quote = FALSE, row.names = FALSE)
  }

  returnList <- list(
    intensityCol = intensityCol,
    intensityColNameOrder = intensityColNameOrder,
    intensityReport = colTable,
    applyFilter = applyFilter,
    filterCon = filterCon,
    filterTime = filterTime,
    filterRep = filterRep
  )

  return(returnList)
}

#' Validate existing or user-uploaded Kinase/Phosphatase data
#'
#' @description This function is designed to identify kinases and phosphatases by matching the peptide IDs. Identified kinases or phosphatases are required to accurately infer a signaling network.
#'
#' @param netPhorceData (Required). netPhorceData calculated from \code{\link{processData}} function.
#' @param defaultKinaseTable (Required). If TRUE, this function will search the preloaded kinase/phosphatase table, which contains the kinase and/or phosphatase IDs of 25 species. Please use \code{\link{checkPreloadKinaseTable}} function to check details of these 25 species.
#' @param userKinaseTable (Optional). A user-provided kinase/phosphatase table, which is required to include an "ID" and "Family" column (CASE-SENSITIVE column names). The ID column includes the peptide's universal ID that matches the protein ID column from the raw MaxQuant data. The Family column will include information regarding whether the peptide is a "Kinase" or "Phosphatase" (CASE-SENSITIVE).
#' @param abbrev (Optional). A three character string of the species abbreviation from the preloaded kinases/phosphatase table. This parameter will be evaluated if the species parameter is NULL.
#' @param species (Required). A string of the species name as listed in the preloaded kinases/phosphatase table. If the species and abbrev parameter are NULL, the function will assume the user is using customized kinase phosphatase data.
#' @param verbose (Optional). If TRUE, print a summary table of matched/unmatched kinases and phosphatases.
#' @return A list of data.frames
#' @export
#' @examples
#' \dontrun{
#' ## Loading One Condition Data
#' data("oneConditionExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = oneConditionExample,
#'                                 	positionCol = "Position",
#'                                 	reverseCol = "Reverse",
#'                                 	localizationProbCol = "Localization prob",
#'                                 	potentialContaminationCol = "Potential contaminant",
#'                                 	aminoAcidCol = "Amino acid",
#'                	                 uniqueIDCol = "Protein",
#'                                 	seqWindowIDCol = "Sequence window",
#'                                 	fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = oneConditionExample,
#'                                      	intensityPattern = "con_time_rep",
#'                                      	verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = oneConditionExample,
#'                          	processedColNames = identifiedCols,
#'               	           processedIntensity = intensityCols,
#'                          	minReplication = 3,
#'                          	minLocalProb = 0.75)
#' ## Validating the Kinase/Phosphatase Information
#' netPhorceData <- validateKinaseTable(netPhorceData = netPhorceData,
#'                                      defaultKinaseTable = TRUE,
#'                                  	  abbrev = "Ath")
#' }
validateKinaseTable <- function(netPhorceData = netPhorceData,
                                defaultKinaseTable = TRUE,
                                userKinaseTable = NULL,
                                abbrev = NULL,
                                species = NULL,
                                verbose = TRUE){
  if(defaultKinaseTable == TRUE){
    data(kinasesPhosphatases, envir=environment())
    Kinases = kinasesPhosphatases
    Kinases = Kinases %>% mutate(ID, sub("_", "-", ID)) %>% mutate(ID = str_to_upper(ID))
    colnames(Kinases) <- toupper(colnames(Kinases))
    Kinases <- Kinases[,c("ABB", "SPECIES", "ID", "FAMILY")]
    # print(head(Kinases))
    kinaseTableCheck = FALSE

    if(is.null(abbrev) & is.null(species)){
      stop("Please enter a valid abbreviate or full species names from the preload
           kinase/phosphatase dataset. Use `checkPreloadKinaseTable()` to see what species are avaliable. ")
    } else if(is.null(species)){
      if(!(abbrev %in% unique(Kinases$ABB))){
        stop("Please enter a valid abbreviate names from the preload
         kinase/phosphatase dataset. Use `checkPreloadKinaseTable()` to see what species are avaliable. ")
      } else {
        Kinases <- subset(Kinases, ABB == abbrev)
        kinaseTableCheck = TRUE
      }
    } else if(is.null(abbrev)){
      if(!(species %in% unique(Kinases$SPECIES))){
        stop("Please enter a valid abbreviate names from the preload
           kinase/phosphatase dataset. Use `checkPreloadKinaseTable()` to see what species are avaliable. ")
      } else {
        Kinases <- subset(Kinases, SPECIES == species)
        kinaseTableCheck = TRUE
      }
    }
  } else {
    if(is.null(userKinaseTable)){
      stop("Please use either preloaded Kinase or Phosphotase table or upload your own. ")
    } else{
      if(c("ID", "FAMILY") %in% colnames(userKinaseTable)){
        Kinases = userKinaseTable
        kinaseTableCheck = TRUE
      } else {
        stop("Please name sure the Kinase/Phosphatase table you uploaded have 'ID' and 'Family' columns. CASE_SENSITIVE")
      }

    }
  }

  if(kinaseTableCheck == TRUE){
    Kinases$ID <- toupper(Kinases$ID)
    # print(head(Kinases))
    data.clustered.avg <-  netPhorceData@data.filtered %>%
      group_by(set, UniqueID, timepoint, condition) %>% filter(obs_tp >= netPhorceData@Design$minReplication) %>%
      summarize(m = mean(normValue, na.rm=TRUE)) %>% ungroup() %>%
      mutate(timepoint = paste(timepoint, condition, sep = "_"), .keep = "unused") %>%
      spread(timepoint, m, fill=0)
    filteredData <- data.clustered.avg %>% dplyr::select(UniqueID) %>%
      # separate(UniqueID, into = c("Protein", "AA", "multiplicity"), sep = "_", remove=FALSE) %>%
      separate(UniqueID, into = c("Part1", "multiplicity"), sep = "_(?=[0-9]+$)") %>%
      separate(Part1, into = c("Protein", "AA"), sep = "_(?=[0-9A-z;]+$)") %>%
      mutate(Protein = sub("(\\.\\d{1,2})$", "", Protein)) %>%
      dplyr::select(-AA, -multiplicity)
    filteredData <- filteredData %>% dplyr::select(Protein) %>% mutate(Protein = str_to_upper(Protein)) %>% group_by(Protein) %>% dplyr::summarise(Total = n())
    # print(head(filteredData))
    totalTable <- merge(filteredData, Kinases, by.x = "Protein", by.y = "ID", all.x = TRUE)
    # print(head(totalTable))
    totalTable <-
      totalTable %>%
      group_by(FAMILY) %>%
      summarise(`Matched Proteins` = n(),
                `Matched Peptides` = sum(Total))
    totalTable$FAMILY[is.na(totalTable$FAMILY)] <- "Unmatched"
    totalTable <- totalTable %>% column_to_rownames("FAMILY")
    # print("TotalTable 2")
    totalTable <- as.data.frame(totalTable)

    requiredRowNames <- c("Kinase", "Phosphatase", "Unmatched")
    dummyRowName <- setdiff(requiredRowNames, rownames(totalTable))
    kinaseTableCheck = FALSE
    if(length(dummyRowName) == 0){
      message("Kinase and Phosphatase Matching Table:")
      print(totalTable)
      kinaseTableCheck = TRUE
    } else if (length(dummyRowName) == 1){
      message("Kinase and Phosphatase Matching Table:")
      if(length(dummyRowName) > 0){
        for(i in seq(length(dummyRowName))){
          tempRow <- data.frame("Matched Proteins" = 0,
                                "Matched Peptides" = 0)
          rownames(tempRow) <- dummyRowName
          colnames(tempRow) <- c("Matched Proteins", "Matched Peptides")
          # print(tempRow)
          totalTable <- rbind(totalTable, tempRow)
        }
      }
      totalTable <- totalTable[requiredRowNames, ]
      print(totalTable)
      kinaseTableCheck = TRUE
    } else {
      totalTable <- data.frame(WARNING = c("No Kinase and Phosphatase Match Found",
                                           "Please check the foramt of the Kinase and Phosphatase File Format Requirement"))
      totalTable
    }
  }


  netPhorceData@Misc$kinaseCheck = list(KinaseSummary = totalTable,
                                        KinaseTable = Kinases,
                                        KinaseCheck = kinaseTableCheck)
  return(netPhorceData)
}


#' Print a summary kinase/phosphatase table
#' @description Check the preloaded kinase/phosphatase table.
#' @param Abbrev (Optional) The function returns the available IDs of the listed abbreviated species. If left as NULL, the function will output the abbreviation and full species name of the 25 species for which inases and/or phosphatases are available.
#' @return data.frame
#' @export
#' @examples
#' \dontrun{
#' ## Print the 25 Abbreviations and Species.
#' checkPreloadKinaseTable()
#' ## Print the IDs of the provided species
#' checkPreloadKinaseTable(Abbrev = "Ath")
#' }
checkPreloadKinaseTable <- function(Abbrev = NULL){
  data(kinasesPhosphatases, envir=environment())
  Kinases = kinasesPhosphatases
  Kinases = Kinases %>% mutate(ID, sub("_", "-", ID)) %>% mutate(ID = str_to_upper(ID))
  colnames(Kinases) <- toupper(colnames(Kinases))
  Kinases <- Kinases[,c("ABB", "SPECIES", "ID", "FAMILY")]
  colnames(Kinases)[1] <- "ABBREVIATION"
  if(is.null(Abbrev)){
    Kinases_Temp <- Kinases %>%
      group_by(ABBREVIATION, SPECIES) %>%
      filter(row_number() == 1) %>% select(ABBREVIATION, SPECIES, ID) %>%
      rename(`Example ID` = ID)
    Kinases %>%
      # select(ABBREVIATION, SPECIES) %>%
      group_by(ABBREVIATION, SPECIES) %>%
      dplyr::summarise(`Number of Records` = n(), .groups = "drop") %>%
      # distinct() %>%
      left_join(Kinases_Temp, by = c("ABBREVIATION", "SPECIES")) %>%
      as.data.frame
  } else {
    if(Abbrev %in% unique(Kinases$ABBREVIATION)){
      Kinases %>%
        filter(ABBREVIATION == Abbrev) %>%
        as_tibble()
    } else {
      stop("The provided abbreviation does not belong to any speceies we have.
           Please check your spelling, and the Abbreviation is case-senstive. ")
      Kinases %>%
        select(ABBREVIATION, SPECIES) %>%
        distinct() %>%
        as.data.frame
    }
  }

}

#' The regulation check based on the input parameter for the Network Analysis
#'
#' @description Preprocessing the NetPhorce data before \code{\link{networkAnalysis}} by setting the thresholds for the minimal up- and down- fold changes, minimal absolute fold change, and q-value.
#' @param netPhorceData (Required). Processed NetPhorceData
#' @param upReg (Required). The percentage fold change required between timepoints to be treated as an increase in phosphorylation.
#' @param downReg (Required). The percentage fold change required between timepoints to be treated as a decrease in phosphorylation.
#' @param absMinThreshold (Required). The bottom percentage of all fold changes that will be considered as unchanged. Setting this threshold will exclude small intensity fold changes that are considered above the upregulation and downregulation thresholds as a result of small median-centered intensity values.
#' @param qValueCutOff (Required). The q-value threshold. A lower threshold increases stringency for including significant phosphopeptides in the network inference algorithm.
#' @param verbose (Optional). If TRUE, a summary table will be printed.
#' @return
#' The netPhorce object with the following slots:
#' \describe{
#'   \item{Design}{Contains data design information and filtering parameters}
#'   \item{data.filtered}{Contains all the data points in a `Long` format that passed filtering criteria.}
#'   \item{data.filtered.aov.summary}{Contains all the data points in a `Long` format with anova results.}
#'   \item{Misc}{Contains accessory  data including default plotting colors and FASTA Keys, if present.}
#'   \item{regulationData}{Contains regulation data calculated through \code{\link{regulationCheck}} function}
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
#' }
regulationCheck <- function(netPhorceData = netPhorceData,
                            upReg = 0.1,
                            downReg = 0.1,
                            absMinThreshold = 0.1,
                            qValueCutOff = 0.05,
                            verbose = FALSE){

  if(netPhorceData@Misc$kinaseCheck$KinaseCheck == FALSE){
    stop("Please make sure there are identified the Kinase/Phosphatase from your sample. ")
  }
  if(upReg < 0 | upReg > 0.5 | !is.numeric(upReg)){
    stop("Please input a valid up regulation cut off, it needs to be a number between 0 and 0.5.")
  }
  if(downReg < 0 | downReg > 0.5 | !is.numeric(downReg)){
    stop("Please input a valid down regulation cut off, it needs to be a number between 0 and 0.5.")
  }
  if(absMinThreshold < 0 | absMinThreshold > 0.5 | !is.numeric(absMinThreshold)){
    stop("Please input a valid absolute minimum threshold cut off, it needs to be a number between 0 and 0.5.")
  }
  if(qValueCutOff < 0 | qValueCutOff > 0.1 | !is.numeric(qValueCutOff)){
    stop("Please input a valid q value cut off, it needs to be a number less than 0.1.")
  }
  returnList <- list()
  newQVal <- qValueCutOff
  data.filtered.aov.summary <- netPhorceData@data.filtered.aov.summary
  ConDesign <- netPhorceData@Design$ConDesign
  uniqueTimes <- netPhorceData@Design$uniqueTimes
  Threshold <- netPhorceData@Design$minReplication
  Kinases <- netPhorceData@Misc$kinaseCheck$KinaseTable
  ## Generate the data.clustered.avg_V2
  # print(netPhorceData@data.filtered )

  data.clustered.avg_V2 = netPhorceData@data.filtered %>%
    group_by(set, UniqueID, timepoint, condition) %>% filter(obs_tp >= Threshold) %>%
    summarize(m = mean(normValue, na.rm=TRUE)) %>% ungroup() %>%
    complete(UniqueID, condition,fill= list(set= "UniqueSet", m=0), timepoint = uniqueTimes$Origin[1]) %>%
    spread(timepoint, m, fill=0)
  # print("Checks")
  # print(data.clustered.avg_V2)

  upReg <- upReg
  downReg <- downReg
  abs_min_threshold <- absMinThreshold
  numberOfDiscreteStates = 1
  # print(paste(c(upReg, downReg, absMinThreshold, qValueCutOff, Threshold), collapse = ";"))
  # print("1")
  data.significant.proteins = data.clustered.avg_V2 %>%
    left_join(data.filtered.aov.summary %>% dplyr::select(UniqueID, qvalue), by = "UniqueID") %>%
    filter(qvalue < newQVal  | set == "UniqueSet") %>%
    gather(experiment, avgValue, -UniqueID, -qvalue, -condition, -set) %>%
    group_by(condition, UniqueID) %>% add_count(avgValue == 0) %>% filter(any(`avgValue == 0` == "FALSE" & n >= 3)) %>%
    dplyr::select(-set, -n, -`avgValue == 0`) %>%
    mutate(tp = as.numeric(str_extract(experiment,"[:digit:]+"))) %>% ### interpret timepoints by digits pattern
    mutate(experiment = fct_reorder(experiment,tp)) ### experiment is only time for now, consider other variables in the future.
  # print("2")
  data.significant.proteins = data.significant.proteins %>%
    separate(UniqueID, into = c("Model_name", "AA", "multiplicity"), sep = "_", remove=FALSE) %>%
    mutate(Model_name = sub("(\\.\\d{1,2})$", "", Model_name)) %>% mutate(Model_name = str_to_upper(Model_name)) %>%
    left_join(Kinases %>% mutate(isKinase = TRUE), by = c("Model_name"="ID")) %>% tidyr::replace_na(list(isKinase = FALSE))
  # print("3")
  data.significant.proteins = data.significant.proteins %>% arrange(condition, experiment) %>% group_by(condition) %>%
    mutate(conditionMedian = median(avgValue[avgValue > 0],na.rm=TRUE)) %>% ungroup() %>%
    group_by(condition, UniqueID) %>%
    mutate(mvalue = avgValue - conditionMedian) %>%
    mutate(change = mvalue - lag(mvalue)) %>% ungroup() %>% group_by(condition) %>%
    mutate(quantile_threshold = quantile(abs(change), probs = abs_min_threshold, na.rm=TRUE)) %>%
    group_by(condition, UniqueID) %>%
    mutate(sigChange = case_when( ### Asymertric thresholds version
      change > pmax(abs(upReg*lag(mvalue)), quantile_threshold ) ~ TRUE,
      change < pmin(-1*abs(downReg*lag(mvalue)), -quantile_threshold ) ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    mutate(sigChangeSign = sign(change) * sigChange) %>%
    group_by(condition, UniqueID) %>% filter(any(sigChangeSign != 0)) %>% ungroup() %>%
    mutate(discreteState = case_when(
      abs(mvalue) >= abs(conditionMedian) ~ 5,
      mvalue == 0 ~ 1,
      TRUE ~ sign(mvalue)*ntile(abs(mvalue), numberOfDiscreteStates)))

  # print(data.significant.proteins)
  # returnList[["data.significant.proteins"]] <- data.significant.proteins
  # returnList[["colors_Reg"]] <-netPhorceData@Misc$colors_Reg

  ## Verbose
  if(verbose == TRUE){
    filteredData <- data.significant.proteins %>%
      filter(!is.na(mvalue) & !is.na(change))
    outputTable <- as.data.frame(table(filteredData$sigChangeSign))
    # outputTable <- t(outputTable)
    # print(outputTable)
    rowNamesTemp <- data.frame(curCol = outputTable[1,])
    truth <- data.frame(curCol = c("-1", "0", "1"),
                        realCol = c("Dephosphorylation", "Unchanged", "Phosphorylation"))
    rowNamesTemp <- merge(rowNamesTemp, truth)
    rownames(outputTable) <- rowNamesTemp$realCol
    colnames(outputTable) <- c("NA", "# Fold change occurances")
    outputTable <- as.data.frame(outputTable)
    outputTable <- outputTable %>% dplyr::select(`# Fold change occurances`)
    print(outputTable)
  }
  netPhorceData@regulationData = data.significant.proteins
  netPhorceData@Misc$dataStats$upReg = upReg
  netPhorceData@Misc$dataStats$downReg = downReg
  netPhorceData@Misc$dataStats$absMinThreshold = absMinThreshold
  netPhorceData@Misc$dataStats$qValueCutOff = qValueCutOff
  return(netPhorceData)
}



































































