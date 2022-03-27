#' Initial processing and quality control of MaxQuant data including the removal of contaminated and reverse peptides, unnecessary columns, and invalid peptides based on filtering criteria.
#'
#' @description The function combines the raw MaxQuant data, and the outputs from \code{\link{confirmIntensityColumns}} and \code{\link{confirmColumnNames}} to form a \code{\link{netPhorce-class}} object.
#' @param rawMaxQuant (Required). Loaded MaxQuant data
#' @param processedColNames (Required). Processed required columns from \code{\link{confirmColumnNames}} function.
#' @param processedIntensity (Required). Processed required intensity columns from the \code{\link{confirmIntensityColumns}} function.
#' @param minReplication (Required). The minimal number of valid replicate data points or replicate zeros per sample a peptide should contain across the entire time course and experiment. Peptides that do not meet these criteria are filtered out.
#' @param minLocalProb (Required). The minimal localization probability required for a peptide to be included in downstream analyses.
#' @return The netPhorce object with the following slots:
#' \describe{
#'   \item{Design}{Contains data design information and filtering parameters}
#'   \item{data.filtered}{Contains all the data points in a `Long` format that passed filtering criteria.}
#'   \item{data.filtered.aov.summary}{Contains the anova results in a `Long` format .}
#'   \item{Misc}{Contains accessory data including default plotting colors and FASTA Keys, if present.}
#' }
#' @export
#' @examples
#' \dontrun{
#' ## Loading One Condition Data
#' data("oneConditionExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = oneConditionExample,
#'                                 	    positionCol = "Position",
#'                                 	    reverseCol = "Reverse",
#'                                 	    localizationProbCol = "Localization prob",
#'                                     	potentialContaminationCol = "Potential contaminant",
#'                                 	    aminoAcidCol = "Amino acid",
#'                    	                uniqueIDCol = "Protein",
#'                                 	    seqWindowIDCol = "Sequence window",
#'                                 	    fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = oneConditionExample,
#'                                        	intensityPattern = "con_time_rep",
#'                                        	verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = oneConditionExample,
#'                          	  processedColNames = identifiedCols,
#'                          	  processedIntensity = intensityCols,
#'                          	  minReplication = 3,
#'                          	  minLocalProb = 0.75)
#' }
processData <- function(rawMaxQuant = rawMaxQuant,
                        processedColNames = processedColNames,
                        processedIntensity = intensity,
                        minReplication = minReplication,
                        minLocalProb = minLocalProb){
  ## Text Progress bar
  pb = txtProgressBar(min = 0, max = 12,
                      initial = 0, style = 3)
  ## Error Checking
  if (!(is.numeric(minReplication))) {
    stop("`Threshold` needs to be a integer greater than 3", call. = FALSE)
  } else {
    if(minReplication < 3){
      stop("`Threshold` needs to be greater than 3", call. = FALSE)
    }
  }
  setTxtProgressBar(pb, 1)

  df <- rawMaxQuant

  ## Rename the columns
  if(!is.null(processedColNames$fastaIDCol)){
    df2 <- df[, c(processedColNames$uniqueID, processedColNames$localizationProbCol, processedColNames$positionCol,
                  processedColNames$reverseCol, processedColNames$potentialContaminationCol, processedColNames$aminoAcidCol,
                  processedColNames$seqWindowID, processedColNames$fastaIDCol, processedIntensity$intensityCol
    )]
    df2 <- df2 %>%
      dplyr::rename("Protein" = processedColNames$uniqueID,
                    "Amino acid" = processedColNames$aminoAcidCol,
                    "Potential contaminant" = processedColNames$potentialContaminationCol,
                    "Localization Probability" = processedColNames$localizationProbCol,
                    "Reverse" = processedColNames$reverseCol,
                    "Position" = processedColNames$positionCol,
                    "FASTA Header" = processedColNames$fastaIDCol,
                    "Sequence Window" = processedColNames$seqWindowID)
  } else {
    print("No FASTA Header Loaded")
    df2 <- df[, c(processedColNames$uniqueID, processedColNames$localizationProbCol, processedColNames$positionCol,
                    processedColNames$reverseCol, processedColNames$potentialContaminationCol, processedColNames$aminoAcidCol,
                    processedColNames$seqWindowID, processedIntensity$intensityCol
    )]
    df2 <- df2 %>%
      dplyr::rename("Protein" = processedColNames$uniqueID,
                    "Amino acid" = processedColNames$aminoAcidCol,
                    "Potential contaminant" = processedColNames$potentialContaminationCol,
                    "Localization Probability" = processedColNames$localizationProbCol,
                    "Reverse" = processedColNames$reverseCol,
                    "Position" = processedColNames$positionCol,
                    "Sequence Window" = processedColNames$seqWindowID)
  }

  df2[, processedIntensity$intensityCol] <- sapply(df2[, processedIntensity$intensityCol], as.numeric)
  nInputSamples = nrow(df2)
  setTxtProgressBar(pb, 2)

  ## Filter based on the localization prob, reverse, and potential.contaminant
  df2 <- df2  %>%
    filter(is.na(Reverse)) %>%
    filter(!is.na(Position)) %>% ### Require a valid position
    filter(is.na(`Potential contaminant`)) %>% ### Remove possible contaminats
    filter(as.numeric(`Localization Probability`) > as.numeric(minLocalProb))

  seqWinData = df2 %>% mutate(AA = paste(`Amino acid`, Position, sep = "")) %>%
    dplyr::select(Protein, `Sequence Window`, AA)

  data.select <- df2 %>%
    mutate(Protein = gsub("_", "-", Protein)) %>%
    mutate(UniqueID = paste(Protein, "_", `Amino acid`, Position, sep = ""))
  # print(head(data.select))

  if(!is.null(processedColNames$fastaIDCol)){
    FASTA_HeaderKey <- data.select[, c("Protein","FASTA Header")]
    # FASTA_HeaderKey$ID <-
    # returnList[["FASTA_HeaderKey"]] <- FASTA_HeaderKey
  } else {
    # returnList[["FASTA_HeaderKey"]] <- NULL
    FASTA_HeaderKey = NULL
  }



  data.select <- data.select %>%
    dplyr::select(UniqueID, starts_with("Intensity") & matches("___\\d"))

  ## Unique Protein Column Validation
  if(anyDuplicated(as.vector(data.select$UniqueID))>0){
    stop("Detect Duplicates, Error in unique protein ID column! ", call. = FALSE)

  } else {
    # message("\nNo duplicated Protein ID found. Proceed to next step: ")
  }

  # print("3")
  setTxtProgressBar(pb, 3)
  intensityColNameOrder <- processedIntensity$intensityColNameOrder
  newOrder <- c("measure", intensityColNameOrder$truth, "Multiplicity")

  data.t = data.select %>% na_if(0) %>% pivot_longer(-UniqueID, values_drop_na=TRUE) %>%
    separate(name, into = c(newOrder), remove=FALSE)
  # print(head(data.t))
  data.t = data.t %>% unite(SampleID,condition,timepoint,replicate,remove=FALSE) %>%
    unite(UniqueID, UniqueID, Multiplicity) %>% dplyr::select(-name,-measure)
  # print(data.t)
  # print(data.t)
  # if(intensity$applyFilter == TRUE){
  #   if()
  #   data.t <- data.t %>%
  #     filter(timepoint %in% intensity$filterTime &
  #              condition %in% intensity$filterCon &
  #              replicate %in% intensity$filterRep)
  # }
  # print(head(data.t))

  # ## Part 4
  # print("4")
  setTxtProgressBar(pb, 4)
  data.filtered = data.t %>%
    add_count(UniqueID, timepoint, condition, name = "obs_tp") %>% ### count up observations that are not zero for each proteinID
    group_by(timepoint) %>% mutate(Reps = (length(unique(replicate)))) %>% ungroup() %>%
    group_by(UniqueID) %>% mutate(obs_tr = length(unique(condition))) %>%
    filter(any(obs_tp >= minReplication)  & all(obs_tp >= minReplication | obs_tp <= (Reps - minReplication))) %>% ungroup() %>%
    filter(obs_tr == max(obs_tr))
  nPostFilter = length(unique(data.filtered$UniqueID))
  # print(data.filtered)
  # print(nPostFilter)

  ## Reorder time points
  uniqueTimes <- data.frame(Origin = unique(data.filtered$timepoint))
  uniqueTimes$Time <- as.numeric(str_extract(uniqueTimes$Origin, "[0-9]+"))
  uniqueTimes <- uniqueTimes %>% arrange(Time)

  ## FullRep
  uniqueTimesFull <- data.frame(Origin = unique(data.filtered$SampleID))
  uniqueTimesFull$Time <- as.numeric(str_extract(uniqueTimesFull$Origin, "[0-9]+"))
  uniqueTimesFull <- uniqueTimesFull %>% arrange(Time)

  ## Obtain the unique colors, color blind friendly
  uniqueConditions <- length(uniqueTimes$Origin)
  if(uniqueConditions <= 9){
    colors = c(
      "#0077BB", # Blue
      "#EE7733", # Orange
      "#33BBEE", # Cyan
      "#CC3311", # Red
      "#009988", # Teal
      "#EE3377", # Magenta
      "#117733", # Green
      "#DDCC77", # Sand
      "#AA4499"  # Purple
    )
  } else{
    colors <- distinctColorPalette(uniqueConditions)
  }

  regColors <- c("#CC3311", "#2A77D5", "#DDCC77")
  colors_V2 = structure(c(colors), class = "palette") # Color Blind Safe color scheme

  colors_HeatmapV1 = c("#E2E6BD","#F1DE81","#F6C971","#EEAB65","#DA8459","#B9534C","#8E133B")
  colors_HeatmapV2 = c("#EFF3FF","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2271B5","#084594")
  regColors <- c("#CC3311", "#2A77D5", "#DDCC77")

  # ## Part 5
  # print("5")
  setTxtProgressBar(pb, 5)
  data_V4 <-
    data.filtered %>%
    mutate(timepoint = factor(timepoint, levels = c(as.vector(uniqueTimes$Origin)))) %>%
    arrange(timepoint) %>%
    dplyr::select(timepoint, replicate) %>%
    group_by(timepoint) %>%
    distinct() %>%
    summarise(ReplicateCount = n())

  data_V5 <-
    data.filtered %>%
    mutate(timepoint = factor(timepoint, levels = c(as.vector(uniqueTimes$Origin))),
           condition = factor(condition, levels = unique(data.filtered$condition))) %>%
    arrange(condition, timepoint) %>%
    dplyr::select(timepoint, condition, replicate) %>%
    group_by(condition, timepoint) %>%
    distinct() %>%
    summarise(ReplicateCount = n()) %>%
    pivot_wider(names_from = timepoint, values_from = ReplicateCount)

  # Part 6
  # print("6")
  setTxtProgressBar(pb, 6)
  ## Preliminary data for the plotting - Original - For vsn normalization
  ConDesign = data.filtered %>% distinct(SampleID) %>%
    # dplyr::select(-UniqueID) %>%
    separate(SampleID, into = c("condition","timepoint","replicate"), remove=FALSE) %>% ### Sample is created using the order exp, time, replicate
    mutate(experiment = paste(condition, timepoint, sep = "_")) %>%
    mutate(label = SampleID, ID = SampleID) %>%
    arrange(timepoint) %>%
    column_to_rownames("SampleID")
  # print(ConDesign)

  # Re-order the time
  # print("7")
  setTxtProgressBar(pb, 7)
  data.filtered <- data.filtered %>%
    mutate(timepoint = factor(timepoint, levels = c(as.vector(uniqueTimes$Origin)))) %>%
    mutate(SampleID = factor(SampleID, levels = c(as.vector(unique(data.filtered$SampleID))))) %>%
    arrange(timepoint)
  # print(data.filtered)

  ## Vsn normalization
  # print("8")
  setTxtProgressBar(pb, 8) ## This step generates vsn2: 3007 x 112 matrix (1 stratum) warning
  data.filtered.spread = data.filtered %>% dplyr::select(UniqueID,SampleID,value) %>% spread(SampleID,value, fill=NA) %>%
    dplyr::select(UniqueID, ConDesign %>% rownames) %>% column_to_rownames("UniqueID") %>% as.data.frame
  row_data = rownames(data.filtered.spread) %>% enframe(value = "ID") %>% mutate(name = ID) %>% as.data.frame
  se <- suppressMessages(SummarizedExperiment(assays = log2(as.matrix(data.filtered.spread)), colData = ConDesign, rowData = row_data))
  se@metadata$formula <- "~ tp + replicate"
  se_norm = suppressMessages(normalize_vsn(se))
  data.norm = assays(se_norm)[[1]] %>% as.data.frame() %>% rownames_to_column("UniqueID") %>% as_tibble()

  # print("9")
  setTxtProgressBar(pb, 9)
  data.filtered = data.filtered %>% ungroup() %>%
    complete(UniqueID,timepoint,condition,fill=list(obs_tp = 0)) %>% ### Complete allows us to consider timepoints with no data
    left_join(data.norm %>% gather(SampleID, normValue, -UniqueID), by = c("UniqueID", "SampleID")) %>%
    group_by(UniqueID) %>% mutate(set = if_else(any(obs_tp < minReplication),"UniqueSet","StatsSet")) %>% ungroup() %>%
    filter(!is.na(SampleID))

  # print("10")
  setTxtProgressBar(pb, 10)
  data.clustered.avg = data.filtered %>%
    group_by(set, UniqueID, timepoint, condition) %>% filter(obs_tp >= minReplication) %>%
    summarize(m = mean(normValue, na.rm=TRUE)) %>% ungroup() %>%
    mutate(timepoint = paste(timepoint, condition, sep = "_"), .keep = "unused") %>%
    spread(timepoint, m, fill=0)

  data.clustered.avg_V2 = data.filtered %>%
    group_by(set, UniqueID, timepoint, condition) %>% filter(obs_tp >= minReplication) %>%
    summarize(m = mean(normValue, na.rm=TRUE)) %>% ungroup() %>%
    complete(UniqueID, condition,fill= list(set= "UniqueSet", m=0), timepoint = uniqueTimes$Origin[1]) %>%
    spread(timepoint, m, fill=0)
  # print(data.clustered.avg)
  # data.clustered.avg$UniqueID <- paste0(data.clustered.avg$UniqueID, "_", data.clustered.avg$experiment)
  ### This is a bit slow.

  # print("11")
  setTxtProgressBar(pb, 11)
  n = n_distinct(ConDesign$condition)
  t = n_distinct(uniqueTimes$Origin)
  message(paste0("\nNumber of Conditions Found: ", n, ". Undergoing a time consuming step, please wait... "))
  ###### . . 1.3.1 !!!!!! Commented out data.filtered.aov ######
  ###### . . 1.3.0 StatsSet checker ######
  statsSetChecker <- data.filtered %>% filter(set == "StatsSet")

  ###### . . 1.3.1 !!!!!! Commented out data.filtered.aov ######
  if(nrow(statsSetChecker) == 0){
      data.filtered.aov = NULL
      data.filtered.aov.summary = NULL
  } else {
      if(n == 1) {
          # print(head(data.filtered))
          data.filtered.aov = data.filtered %>% filter(set == "StatsSet") %>% group_by(UniqueID) %>% nest() %>% ungroup() %>%
              mutate(AOV = map(data, ~ aov(normValue ~ replicate + timepoint + Error(replicate), data=.x)) ) %>%
              mutate(result = map(AOV,broom::tidy))
          data.filtered.aov.summary = data.filtered.aov %>% dplyr::select(UniqueID,result) %>% unnest(cols = c(result)) %>% filter(!is.na(p.value)) %>%
              mutate(qvalue = qvalue(p.value, lfdr.out = TRUE)[["qvalues"]])
      } else {
          if (n > 1) {
              if (t < 3) {
                  noNetwork <- TRUE
                  cat("There are less than 3 time points detected, the subsequent network analysis portion
                        is disabled. Please use at least 3 time points to enable the network calculation.")
                  data.filtered.aov = data.filtered %>% filter(set == "StatsSet") %>%
                      group_by(UniqueID) %>% nest() %>% ungroup() %>%
                      mutate(AOV = map(data, ~ aov(normValue ~ replicate + condition + Error(replicate), data=.x)) ) %>%
                      mutate(result = map(AOV,broom::tidy))
                  data.filtered.aov.summary = data.filtered.aov %>% dplyr::select(UniqueID,result) %>% unnest(cols = c(result)) %>% filter(!is.na(p.value)) %>%
                      mutate(qvalue = qvalue(p.value, lfdr.out = TRUE)[["qvalues"]])
              } else {
                  data.filtered.aov = data.filtered %>% filter(set == "StatsSet") %>%
                      group_by(UniqueID) %>% nest() %>% ungroup() %>%
                      mutate(AOV = map(data, ~ aov(normValue ~ replicate + timepoint*condition + Error(replicate), data=.x)) ) %>%
                      mutate(result = map(AOV,broom::tidy))
                  data.filtered.aov.summary = data.filtered.aov %>% dplyr::select(UniqueID,result) %>%
                    unnest(cols = c(result)) %>% filter(!is.na(p.value)) %>%
                      filter(term == "condition") %>% mutate(qvalue = qvalue(p.value, lfdr.out = TRUE)[["qvalues"]])

              }
          }
      }
  }

  dataStats = data.frame(
    nConditions = n_distinct(ConDesign$condition),
    nTimePoints = n_distinct(ConDesign$timepoint),
    nReplicates = n_distinct(ConDesign$replicate),
    nInputSamples = nInputSamples,
    nPostFilter = nPostFilter,
    qValueCutOff = 0.1,
    upReg = NA,
    downReg = NA,
    absMinThreshold = NA
  )
  if(nrow(data.filtered %>% filter(set == "StatsSet"))>0){
    # dataStats$nStatsSet = nrow(length(unique(data.filtered %>% filter(set == "StatsSet") %>% pull("UniqueID"))))
    dataStats$nStatsSetqCut = nrow(data.filtered.aov.summary)
  } else {
    # dataStats$nStatsSet = 0
    dataStats$nStatsSetqCut = 0
  }
  if(nrow(data.filtered %>% filter(set == "UniqueSet"))>0){
    dataStats$nPresenceAbsence = length(unique(data.filtered %>% filter(set == "UniqueSet") %>% pull("UniqueID")))
  } else {
    dataStats$nPresenceAbsence = 0
  }

  setTxtProgressBar(pb, 12)
  message("\nComplete.")
  Design = list("ConDesign" = ConDesign,
                "uniqueTimesFull" = uniqueTimesFull,
                "uniqueTimes" = uniqueTimes,
                "intensityColNameOrder" = intensityColNameOrder,
                "minReplication" = minReplication,
                "minLocalProb" = minLocalProb)

  Misc = list("colors" = colors,
              "colors_V2" = colors_V2,
              "colors_HeatmapV1" = colors_HeatmapV1,
              "colors_HeatmapV2" = colors_HeatmapV2,
              "colors_Reg" = regColors,
              "FASTA_HeaderKey" = FASTA_HeaderKey,
              "dataStats" = dataStats,
              "seqWinData" = seqWinData)
  returnData <- new(
    "netPhorce",
    Design = Design,
    Misc = Misc,
    data.filtered = data.filtered,
    data.filtered.aov.summary = data.filtered.aov.summary
  )
  return(returnData)
}

