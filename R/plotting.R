#' Plot the Distribution of the Intensity Data
#'
#' @description Boxplots of the intensity data for each time point, condition, and replicate before and after variance stabilizing normalization (vsn). The function generates a ggplot object by default. An interactive plotly version (used on our Shiny application) is included if \emph{plotly = TRUE} is set.
#' @param netPhorceData (Required). Processed netPhorce Object
#' @param condition (Required). Select a specific condition from your experiment.
#' @param plotly (Required). If TRUE, output an interactive \code{\link{plotly}} version, else output a static \code{\link{ggplot2}} version.
#' @return a ggplot/plotly object
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Generate Distribution Plot – GGPLOT version (Static)
#' plotDistribution(netPhorceData = netPhorceData,
#'                  condition = NULL,
#'                  plotly = FALSE)
#' ## Generate Distribution Plot – PLOTLY version (Interactive)
#' plotDistribution(netPhorceData = netPhorceData,
#'                  condition = NULL,
#'                  plotly = TRUE)
#' }
plotDistribution = function(netPhorceData = netPhorceData,
                            condition = NULL,
                            plotly = FALSE){
  ## Requested Data Objects from NetPhore
  if(!is.null(netPhorceData@data.filtered)){
    data.filtered = netPhorceData@data.filtered
    ConDesign = netPhorceData@Design$ConDesign
    colors = netPhorceData@Misc$colors
    uniqueTimesFull = netPhorceData@Design$uniqueTimesFull

    nConditions = n_distinct(ConDesign$condition)
    if(is.null(condition)){
      if(nConditions == 1){
        curCondition = unique(ConDesign$condition)
      } else {
        cat("Multiple condition detected, please use the 'condition = ' arguemnt to select the specific condition.\n")
        cat("The following plot is based on the first found condition which is: \n")
        curCondition = unique(ConDesign$condition)[1]
        cat(curCondition)
      }
    } else {
      curCondition = condition
      if(!(curCondition %in% ConDesign$condition)){
        stop("Selected condition is not in the current data, the avaliable conditions are: ",
             paste(unique(ConDesign$condition), collapse = "; "), call. = FALSE)
      }
    }

    data.filtered.spread = data.filtered %>%
      dplyr::select(UniqueID,SampleID,value) %>%
      spread(SampleID,value, fill=NA) %>%
      dplyr::select(UniqueID, ConDesign %>% rownames) %>%
      column_to_rownames("UniqueID") %>%
      as.data.frame

    row_data = rownames(data.filtered.spread) %>%
      enframe(value = "ID") %>%
      mutate(name = ID) %>%
      as.data.frame

    se <- suppressMessages(SummarizedExperiment(assays = log2(as.matrix(data.filtered.spread)),
                               colData = ConDesign,
                               rowData = row_data)
    )
    se@metadata$formula <- "~ tp + replicate"
    se_norm = suppressMessages(normalize_vsn(se))

    ConDesign_V2 <- ConDesign %>%
      distinct(timepoint, experiment) %>%
      distinct()
    ConDesign_V2$ConTime <- as.numeric(str_extract(ConDesign_V2$timepoint, "[0-9]+"))
    ConDesign_V2$condition <- sapply(strsplit(ConDesign_V2$experiment, "_|\\s"), "[", 1)
    ConDesign_V2 <- ConDesign_V2 %>%
      arrange(condition, ConTime)
    ConDesign_V2 <- subset(ConDesign_V2, condition == curCondition)

    sePlot <- plot_normalization(se = se, se_norm = se_norm,
                                 filterCondition = curCondition)
    sePlot[["data"]][["ID"]] <- factor(sePlot[["data"]][["ID"]], levels = rev(uniqueTimesFull$Origin))

    orderTime <- unique(ConDesign_V2$experiment)
    sePlot[["data"]][["experiment"]] <- factor(sePlot[["data"]][["experiment"]], levels = orderTime)
    sePlot[["data"]][["var"]] <- gsub("se_norm", "Distribution of Normalized Data", sePlot[["data"]][["var"]])
    sePlot[["data"]][["var"]] <- gsub("se", "Distribution of Raw Data", sePlot[["data"]][["var"]])
    sePlot[["data"]][["var"]] <- factor(sePlot[["data"]][["var"]],
                                        levels = c("Distribution of Raw Data", "Distribution of Normalized Data"))
    sePlot <- sePlot +
      scale_fill_manual(values = colors)

    if(plotly == FALSE){
      sePlot +
        guides(fill=guide_legend(title="Time Points")) +
        labs(title = paste0("Condition: ", curCondition))
    } else {
      sePlot <- sePlot +
        theme(legend.title = element_blank(),
              axis.title.x=element_blank())
      ggplotly(sePlot) %>%
        layout(legend = list(title = list(text = '<b> Time Points </b>')),
               xaxis  = list(title = "<b>log<sub>2</sub> Intensity</b>"),
               title = list(
                 text = paste0("<b>Condition: ", curCondition, "</b>"),
                 x = 0.02,
                 font = list(
                   size = 15
                 )
               ))
    }
  } else {
    print("NetPhorce data was incomplete. Please process the data again")
  }
}


#' PCA Plot
#'
#' @description The principle component analysis (PCA) output before and after variance stabilizing normalization (vsn).
#' @param netPhorceData (Required). Processed NetPhorce Object
#' @param condition (Required). Select a specific condition from your experiment.
#' @param normalized (Required). Use raw or normalized intensity data
#' @param plotly (Required). If TRUE, output an interactive \code{\link{plotly}} version, else output a static \code{\link{ggplot2}} version.
#' @return A ggplot/plotly object
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(
#'   rawMaxQuant = twoConditionsExample,
#'   positionCol = "Position",
#'   reverseCol = "Reverse",
#'   localizationProbCol = "Localization prob",
#'   potentialContaminationCol = "Potential contaminant",
#'   aminoAcidCol = "Amino acid",
#'   uniqueIDCol = "Protein",
#'   seqWindowIDCol = "Sequence window",
#'   fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Generate Distribution Plot – GGPLOT version (Static)
#' plotPCA(netPhorceData = netPhorceData,
#'         condition = "tot3",
#'         plotly = FALSE)
#' ## Generate Distribution Plot – PLOTLY version (Interactive)
#' plotPCA(netPhorceData = netPhorceData,
#'         condition = "Col0",
#'         plotly = TRUE)
#' }
plotPCA <- function(netPhorceData = netPhorceData,
                    condition = NULL,
                    normalized = FALSE,
                    plotly = FALSE){
  ## Requested Data Objects from NetPhorce
  if(!is.null(netPhorceData@data.filtered)){
    data.filtered = netPhorceData@data.filtered
    ConDesign = netPhorceData@Design$ConDesign
    colors = netPhorceData@Misc$colors
    uniqueTimes = netPhorceData@Design$uniqueTimes
    intensityColNameOrder = netPhorceData@Design$intensityColNameOrder

    nConditions = n_distinct(ConDesign$condition)
    if(is.null(condition)){
      if(nConditions == 1){
        curCondition = unique(ConDesign$condition)
      } else {
        cat("Multiple condition detected, please use the 'condition = ' arguemnt to select the specific condition.\n")
        cat("The following plot is based on the first found condition which is: \n")
        curCondition = unique(ConDesign$condition)[1]
        cat(curCondition)
      }
    } else {
      curCondition = condition
      if(!(curCondition %in% ConDesign$condition)){
        stop("Selected condition is not in the current data, the avaliable conditions are: ",
             paste(unique(ConDesign$condition), collapse = "; "), call. = FALSE)
      }
    }

    data.filtered.spread = data.filtered %>%
      dplyr::select(UniqueID,SampleID,value) %>%
      spread(SampleID,value, fill=NA) %>%
      dplyr::select(UniqueID, ConDesign %>% rownames) %>%
      column_to_rownames("UniqueID") %>%
      as.data.frame

    if(normalized == FALSE){
      ##PCA unnorm data
      rawPCAData = as.data.frame(t(data.filtered.spread)) %>% rownames_to_column(var = "SampleID") %>%
        mutate_if(is.numeric, ~log(.,2)) %>%
        separate(SampleID, into = c("condition", "timepoint", "replicate"), remove=FALSE) %>% replace(is.na(.), 0)
      titleText = "Unnormalized"
    } else {
      ##PCA norm data
      row_data = rownames(data.filtered.spread) %>%
        enframe(value = "ID") %>%
        mutate(name = ID) %>%
        as.data.frame
      se <- suppressMessages(SummarizedExperiment(assays = log2(as.matrix(data.filtered.spread)),
                                                  colData = ConDesign,
                                                  rowData = row_data)
      )
      se@metadata$formula <- "~ tp + replicate"
      se_norm = suppressMessages(normalize_vsn(se))
      data.norm = assays(se_norm)[[1]] %>% as.data.frame() %>% rownames_to_column("UniqueID") %>% as_tibble()

      rawPCAData = as.data.frame(t(data.norm %>% column_to_rownames("UniqueID"))) %>%
        rownames_to_column(var = "SampleID") %>%
        separate(SampleID, into = c("condition", "timepoint", "replicate"), remove=FALSE) %>% replace(is.na(.), 0)
      titleText = "Normalized"
    }

    rawPCA_V2 <- subset(rawPCAData, condition == curCondition)

    pcaData <- prcomp(rawPCA_V2[, (5:ncol(rawPCA_V2))], center = TRUE, scale. = TRUE)

    finalPCAData <- as.data.frame(pcaData$x[,1:2]) %>% cbind(rawPCA_V2[1:4])
    finalPCAData$timepoint <- factor(finalPCAData$timepoint, levels = uniqueTimes$Origin)
    finalPCAData_Percentage <- paste0(round(pcaData$sdev[1:2], 2), "%")
    nTime <- length(unique(finalPCAData$timepoint))


    if(plotly == FALSE){
      ggplot(data = finalPCAData, aes(x = PC1, y = PC2, colour = timepoint)) +
        stat_ellipse(aes(x = PC1, y = PC2, fill = timepoint),
                     geom = "polygon", alpha = 0.2,  type = "norm", show.legend = FALSE) +
        geom_point(aes(color = timepoint))  +
        scale_color_manual(values = colors[1:nTime]) +
        scale_fill_manual(values = colors[1:nTime]) +
        labs(x = paste0("PCA1 (", finalPCAData_Percentage[1], ")"),
             y = paste0("PCA1 (", finalPCAData_Percentage[2], ")"),
             title = paste0(titleText, " PCA Plot | Condition: ", curCondition))  +
        guides(colour=guide_legend(title="Time Points")) + theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    } else {
      ellipseData <- NULL
      uniqueTime <- unique(finalPCAData$timepoint)
      for(i in seq_along(uniqueTime)){
        subData <- subset(finalPCAData, timepoint == uniqueTime[i])
        centers <- c(mean(subData[,1]), mean(subData[,2]))
        ellipseDataTemp <- as.data.frame(ellipse::ellipse(cov(subData[,1:2]), centre = centers))
        ellipseDataTemp$timepoint <- uniqueTime[i]
        if(i == 1){
          ellipseData <- ellipseDataTemp
        } else {
          ellipseData <- rbind(ellipseData, ellipseDataTemp)
        }
      }
      finalPCAData$Text <- paste0("<b>Replication</b>: ", finalPCAData$replicate, "<br>",
                                  "<b>Time Point</b>:  ", finalPCAData$timepoint)
      plot_ly() %>%
        add_polygons(data = as.data.frame(ellipseData),
                  x = ~PC1,
                  y = ~PC2,
                  color = ~timepoint,
                  colors = colors[1:nTime],
                  # type = "area",
                  fill = "toself", legendgroup = ~timepoint, showlegend = FALSE,
                  opacity = 0.4) %>%
        add_trace(data = finalPCAData,
                  x = ~PC1,
                  y = ~PC2,
                  color = ~timepoint,
                  colors = colors[1:nTime],
                  type = "scatter",
                  mode = 'markers',
                  hovertext = ~Text,
                  legendgroup = ~timepoint) %>%
        layout(
          title = paste0("<b>", titleText, " PCA Plot | Condition: ", curCondition, "</b>"),
          xaxis = list(
            title = paste0("<b>PC1 (", finalPCAData_Percentage[1], ")</b>")
          ),
          yaxis = list(
            title = paste0("<b>PC2 (", finalPCAData_Percentage[2], ")</b>")
          ),
          legend = list(title = list(text = "<b>Time Points</b>"))
        )
    }
  } else {
    print("NetPhorce data was incomplete. Please process the data again")
  }
}

#' Plot the Distribution of the Intensity Data
#'
#' @description The histogram and/or boxplot of the intensity data for each time point after normalization. The function generates a \code{\link{ggplot}} plot.
#' @param netPhorceData (Required). Processed NetPhorce Object
#' @param condition (Required). Select a specific condition from your experiment.
#' @param histogram (Optional). Enable the histogram Plot. If the boxplot parameter is also TRUE, a side-by-side plot will be generated.
#' @param boxplot (Optional). Enable the boxplot. If the histogram parameter is also TRUE, a side-by-side plot will be generated.
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(
#'   rawMaxQuant = twoConditionsExample,
#'   positionCol = "Position",
#'   reverseCol = "Reverse",
#'   localizationProbCol = "Localization prob",
#'   potentialContaminationCol = "Potential contaminant",
#'   aminoAcidCol = "Amino acid",
#'   uniqueIDCol = "Protein",
#'   seqWindowIDCol = "Sequence window",
#'   fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Side-by-side Plot
#' plotHistBox(netPhorceData = netPhorceData,
#'              condition = "tot3",
#'              histogram = TRUE,
#'              boxplot = TRUE)
#' ## Boxplot Only
#' plotHistBox(netPhorceData = netPhorceData,
#'              condition = "Col0",
#'              histogram = FALSE,
#'              boxplot = TRUE)
#' ## Histogram Only
#' plotHistBox(netPhorceData = netPhorceData,
#'              condition = "tot3",
#'              histogram = TRUE,
#'              boxplot = FALSE)
#' }
plotHistBox <- function(netPhorceData = netPhorceData,
                        condition = NULL,
                        histogram = TRUE,
                        boxplot = TRUE){
  ## Requested Data Objects from NetPhorce
  if(!is.null(netPhorceData@data.filtered)){
    data.filtered = netPhorceData@data.filtered
    ConDesign = netPhorceData@Design$ConDesign
    colors = netPhorceData@Misc$colors
    uniqueTimes = netPhorceData@Design$uniqueTimes
    intensityColNameOrder = netPhorceData@Design$intensityColNameOrder

    if(histogram == FALSE & boxplot == FALSE){
      stop("At least histogram or boxplot needs to be TRUE for plotting", call. = FALSE)
    }

    nConditions = n_distinct(ConDesign$condition)
    if(is.null(condition)){
      if(nConditions == 1){
        curCondition = unique(ConDesign$condition)
      } else {
        cat("Multiple condition detected, please use the 'condition = ' arguemnt to select the specific condition.\n")
        cat("The following plot is based on the first found condition which is: \n")
        curCondition = unique(ConDesign$condition)[1]
        cat(curCondition)
      }
    } else {
      curCondition = condition
      if(!(curCondition %in% ConDesign$condition)){
        stop("Selected condition is not in the current data, the avaliable conditions are: ",
             paste(unique(ConDesign$condition), collapse = "; "), call. = FALSE)
      }
    }

    data.filtered.spread = data.filtered %>%
      dplyr::select(UniqueID,SampleID,value) %>%
      spread(SampleID,value, fill=NA) %>%
      dplyr::select(UniqueID, ConDesign %>% rownames) %>%
      column_to_rownames("UniqueID") %>%
      as.data.frame

    row_data = rownames(data.filtered.spread) %>%
      enframe(value = "ID") %>%
      mutate(name = ID) %>%
      as.data.frame

    se <- suppressMessages(SummarizedExperiment(assays = log2(as.matrix(data.filtered.spread)),
                                                colData = ConDesign,
                                                rowData = row_data)
    )
    se@metadata$formula <- "~ tp + replicate"
    se_norm = suppressMessages(normalize_vsn(se))
    data.norm = assays(se_norm)[[1]] %>% as.data.frame() %>% rownames_to_column("UniqueID") %>% as_tibble()

    histData <- data.norm %>%
      gather(key = "variable",
             value = "Value",
             colnames(data.norm[2:(ncol(data.norm))]))
    histData$TimePoint <- sapply(strsplit(as.vector(histData$variable), "_"), "[", 2)
    histData$TimePoint <- factor(histData$TimePoint,
                                 levels = (uniqueTimes$Origin))
    histData$Replication <- sapply(strsplit(as.vector(histData$variable), "_"), "[", 3)
    histData$condition <- sapply(strsplit(as.vector(histData$variable), "_"), "[", 1)

    histData <- subset(histData, condition == curCondition)

    nTimes <- length(unique(histData$TimePoint))


    if(histogram == TRUE & boxplot == FALSE){
      ggplot(histData) +
        geom_histogram(aes(Value, fill = TimePoint), bins = 200) +
        facet_grid(TimePoint ~ .) +
        labs(x = expression(paste(Log[2]*" Intensity")),
             y = "Frequency",
             title = paste0("Condition: ", curCondition)) +
        theme_minimal() +
        theme(strip.text.y = element_blank(),
              text = element_text(size = 15),
              plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values = (colors[1:nTimes]),
                          breaks = uniqueTimes$Origin)
    } else if(histogram == FALSE & boxplot == TRUE){
      ggplot(histData) +
        geom_boxplot(aes(x = TimePoint, y = Value, fill = TimePoint)) +
        labs(y = expression(paste(Log[2]*" Intensity")),
             title = paste0("Condition: ", curCondition)) +
        coord_flip() +
        facet_grid(TimePoint ~ ., scales = "free") +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              text = element_text(size = 15),
              plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values = (colors[1:nTimes]),
                          breaks = uniqueTimes$Origin)
    } else if(histogram == TRUE & boxplot == TRUE){
      histPlot <- ggplot(histData) +
        geom_histogram(aes(Value, fill = TimePoint), bins = 200) +
        facet_grid(TimePoint ~ .) +
        labs(x = expression(paste(Log[2]*" Intensity")), y = "") +
        theme_minimal() +
        theme(strip.text.y = element_blank(), text = element_text(size = 15)) +
        scale_fill_manual(values = (colors[1:nTimes]),
                          breaks = uniqueTimes$Origin)


      boxPlot <- ggplot(histData) +
        geom_boxplot(aes(x = TimePoint, y = Value, fill = TimePoint)) +
        labs(y = expression(paste(Log[2]*" Intensity"))) +
        coord_flip() +
        facet_grid(TimePoint ~ ., scales = "free") +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(), text = element_text(size = 15)) +
        scale_fill_manual(values = (colors[1:nTimes]),
                          breaks = uniqueTimes$Origin)

      histBoxPlot <- ggarrange(histPlot, boxPlot,
                               widths = c(1,1),
                               common.legend = TRUE, legend = "bottom")
      histBoxPlot <- annotate_figure(histBoxPlot,
                                     top = text_grob(paste0("Condition: ", curCondition),
                                                     face = "bold", size = 14))
      histBoxPlot
    }
  } else {
    print("NetPhorce data was incomplete. Please process the data again")
  }
}


#' Confirm or List Unique IDs from the netPhorceData
#'
#' @description This function intended to use to confirm the user’s input Unique IDs for plotting or list the avaliable Unique IDs.
#' @param netPhorceData (Required). Processed netPhorce Object
#' @param uniqueIDList (Required). Vector of unique IDs. If NULL, this function will return all avaliable uniqueIDs.
#' @param verbose (Required). If TRUE, the function will return the available Unique IDs if the user left uniqueIDList as NULL or none of the IDs on the uniqueIDList matched the IDs in the netPhorce data.
#' @return list of vectors for matched and unmatched IDs
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Search for matching IDs
#' findUniqueIDs(netPhorceData = netPhorceData,
#'               uniqueIDList = c(
#'                 # Significant Set Examples:
#'                 "AT1G13030.1", "AT1G13360.3", "AT1G42550.1",
#'                 # Unique/Abscence Set Examples
#'                 "AT1G17280.9", "AT1G22310.2", "AT1G23890.2"
#'               ),
#'               verbose = TRUE)
#' ## Return all available IDs
#' findUniqueIDs(netPhorceData = netPhorceData,
#'               uniqueIDList = NULL)
#' }
findUniqueIDs <- function(netPhorceData = netPhorceData,
                          uniqueIDList = NULL,
                          verbose = FALSE){
  if(!is.null(netPhorceData@data.filtered)){
    data.filtered = netPhorceData@data.filtered
    data.filtered_UniqueID_Table <- data.frame(Full = data.filtered$UniqueID)
    data.filtered_UniqueID_Table$Short <- sapply(strsplit(as.character(data.filtered_UniqueID_Table$Full), "_"), `[`, 1)

    if(is.null(uniqueIDList)){
      if(verbose){
        print("Print all avaliable Unique ID by default (versboe == TRUE): ")
        print(unique(data.filtered_UniqueID_Table$Short))
      }
      stop("No unique id was provided, returning a list of avaliable unique IDs if verbose = TRUE", call. = FALSE)
    } else {
      if(all(uniqueIDList %in% unique(data.filtered_UniqueID_Table$Short))){
        cat("All ", length(uniqueIDList), " provided unqiue IDs are matched to existing Protein IDs\n")
        return(list("MatchedIDs" = uniqueIDList,
                    "UnMatchedIDs" = NULL))
      } else if (any(uniqueIDList %in% unique(data.filtered_UniqueID_Table$Short))){
        cat(paste0("Only ", length(unique(uniqueIDList)[uniqueIDList %in% unique(data.filtered_UniqueID_Table$Short)]),
                   " out of ", length(uniqueIDList), " provided IDs matched to existing Protein IDs\n"))
        return(list("MatchedIDs" = unique(uniqueIDList)[uniqueIDList %in% unique(data.filtered_UniqueID_Table$Short)],
                    "UnMatchedIDs" = unique(uniqueIDList)[!(uniqueIDList %in% unique(data.filtered_UniqueID_Table$Short))]))
      } else {
        if(verbose){
          print(unique(data.filtered_UniqueID_Table$Short))
        }
        stop("No unique ids matched, returning a list of avaliable unique IDs if verbose = TRUE", call. = FALSE)
      }
    }
  } else {
    print("NetPhorce data was incomplete. Please process the data again")
  }
}

#' Identify and Validate Provided Clusters for plotUniqueIDsHeatmaps
#'
#' @description Return matched clusters that are avliable for plotting using \code{\link{plotClustersHeatmap}} function. The user can leave the clusterIDs as `NULL` to export a full table of all the avliable cluster and its corresponding peptide IDs. After identifying the uniqueID, the user can then perceed with plotting by the unique ID or by the cluster in which the peptide belongs to.
#' @param netPhorceData (Required). Processed NetPhorce Object
#' @param clusterIDs (Required). Select a specific cluster number.
#' @param heatmapType (Required). Either "Significant" or "AbsencePresence" are required. "Significant" for the unique peptide IDs that were significantly differentially phosphorylated peptides. "AbsencePresence" for the unique peptide IDs that were absent at least one-time point.
#' @param minQVal (Optional). A q-value threshold only for heatmapType == `Significant`. A lower threshold increases stringency for displaying significant phosphopeptides into the heatmap. The q-value controls the positive false discovery rate and is estimated based on the p-values obtained from a linear mixed model fitted to the intensities of each peptide.
#' @param verbose If TRUE, the function will return all the peptides along with its cluster number in a table format when the clusterIDs is `NULL` or the input cluster is not found.
#' @return list of data.frames
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Returns all available peptides along with IDs from Significant Set
#' clusterIDs_Sig <-
#'   findClusters(netPhorceData = netPhorceData,
#'                clusterIDs = c(0),
#'                heatmapType = "Significant",
#'                minQVal = 0.1,
#'                verbose = TRUE)
#' ## Returns all available peptides Identify Cluster IDs from “AbsencePresence” Set
#' clusterIDs_AbsPrs <-
#'   findClusters(netPhorceData = netPhorceData,
#'                clusterIDs = c(1, 2),
#'                heatmapType = "AbsencePresence",
#'                verbose = TRUE)
#' }
findClusters <- function(netPhorceData = netPhorceData,
                         clusterIDs = NULL,
                         heatmapType = "Significant",
                         minQVal = 0.05,
                         verbose = TRUE){
  if(!is.null(netPhorceData@data.filtered)){
    data.filtered = netPhorceData@data.filtered
    Threshold = netPhorceData@Design$minReplication
    data.filtered.aov.summary = netPhorceData@data.filtered.aov.summary
    data.clustered.avg = data.filtered %>%
      group_by(set, UniqueID, timepoint, condition) %>% filter(obs_tp >= Threshold) %>%
      summarize(m = mean(normValue, na.rm=TRUE)) %>% ungroup() %>%
      mutate(timepoint = paste(timepoint, condition, sep = "_"), .keep = "unused") %>%
      spread(timepoint, m, fill=0)

    if(length(heatmapType) > 1){
      stop("Please choose a single type for heatmap from 'AbsencePresence' or 'Significant'", call. = FALSE)
    }

    if(any(heatmapType %in% c("AbsencePresence", "Significant"))){
      if(heatmapType == "AbsencePresence"){
        curSetName = "UniqueSet"
      } else if(heatmapType == "Significant"){
        curSetName = "StatsSet"
      }
    } else {
      stop("Please choose a single type for heatmap from 'AbsencePresence' or 'Significant'", call. = FALSE)
    }

    if(curSetName == "StatsSet"){
      if(!any(data.clustered.avg$set == curSetName)){
        stop("No significant peptides are found. ", call. = FALSE)
      } else {
        heatmap_Significant_phosphosites.summary = data.clustered.avg %>%
          inner_join(data.filtered.aov.summary %>% filter(qvalue < minQVal) %>% dplyr::select(UniqueID)) %>%
          dplyr::select(-set) %>% distinct() %>%
          nest(data = everything()) %>% ungroup() %>%
          mutate(dist = map(data, ~ .x %>% column_to_rownames("UniqueID") %>% dist())) %>%
          mutate(hclust = map(dist, hclust,  method = "complete")) %>%
          mutate(order = map(hclust,"order")) %>%
          mutate(cutreeDynamic = map2(hclust, dist, ~ cutreeDynamic(dendro = .x, distM = as.matrix(.y),
                                                                    method = "hybrid", deepSplit = 4)))
        UniqueID_Table <- data.frame(Full = heatmap_Significant_phosphosites.summary$data[[1]]$UniqueID)
        UniqueID_Table$Cluster <- heatmap_Significant_phosphosites.summary$cutreeDynamic[[1]]
        UniqueID_Table$Short <- sapply(strsplit(as.character(UniqueID_Table$Full), "_"), `[`, 1)
        clusters <- sort(unique(UniqueID_Table$Cluster))
        if(all(clusterIDs %in% clusters) & !is.null(clusterIDs)){
          cat("All ", length(clusterIDs), " cluster are matched with provided q-value cut-off: ", minQVal, "\n")
          UniqueID_Table_Filtered = UniqueID_Table %>%
            filter(Cluster %in% clusterIDs) %>% rename(PeptideID = Full, UniqueID = Short)
          return(list("MatchedClusters" = clusterIDs,
                      "UnMatchedClusters" = NULL,
                      "MatchedClusterTable" = UniqueID_Table_Filtered,
                      "heatmapType" = "StatsSet",
                      "minQVal" = minQVal))
        } else if (any(clusterIDs %in% clusters)){
          cat(paste0("Only ", length(unique(clusterIDs)[clusterIDs %in% clusters]),
                     " out of ", length(clusterIDs), " cluster are matched with provided q-value cut-off: ", minQVal, "\n"))
          UniqueID_Table_Filtered = UniqueID_Table %>%
            filter(Cluster %in% unique(clusterIDs)[clusterIDs %in% clusters]) %>%
            rename(PeptideID = Full, UniqueID = Short)

          return(list("MatchedClusters" = unique(clusterIDs)[clusterIDs %in% clusters],
                      "UnMatchedClusters" = unique(clusterIDs)[!(clusterIDs %in% clusters)],
                      "MatchedClusterTable" = UniqueID_Table_Filtered,
                      "heatmapType" = "StatsSet",
                      "minQVal" = minQVal))
        } else {
          if(verbose==TRUE){
            cat("The avaliable Clusters are: \n")
            cat(paste(clusters, collapse = ", "), "\n")
            cat("Please use selected cluster(s) using the clusterIDs arguments,\n")
            cat("With a total of", nrow(UniqueID_Table), " peptides. \n")
            fullTable = UniqueID_Table %>% rename(PeptideID = Full, UniqueID = Short) %>%
              group_by(Cluster, UniqueID) %>% arrange(Cluster) %>%
              as_tibble()

            return(fullTable)
          } else{
            cat("No cluster was provided, returning a list of avaliable clusters if verbose = TRUE")
          }
        }
      }
    }

    if(curSetName == "UniqueSet"){
      if(!any(data.clustered.avg$set == curSetName)){
        stop("No significant peptides are found. ", call. = FALSE)
      } else {
        heatmap_Unique_phosphosites.summary = data.clustered.avg %>% filter(set == "UniqueSet") %>%
          dplyr::select(-set) %>% nest(data = everything()) %>% ungroup() %>%
          distinct() %>%
          mutate(dist = map(data, ~ .x %>% column_to_rownames("UniqueID") %>% dist())) %>%
          mutate(hclust = map(dist, hclust,  method = "complete")) %>%
          mutate(order = map(hclust,"order")) %>%
          mutate(cutreeDynamic = map2(hclust, dist, ~ cutreeDynamic(dendro = .x, distM = as.matrix(.y),
                                                                    method = "hybrid", deepSplit = 4)))

        UniqueID_Table_Unique <- data.frame(Full = heatmap_Unique_phosphosites.summary$data[[1]]$UniqueID)
        UniqueID_Table_Unique$Cluster <- heatmap_Unique_phosphosites.summary$cutreeDynamic[[1]]
        UniqueID_Table_Unique$Short <- sapply(strsplit(as.character(UniqueID_Table_Unique$Full), "_"), `[`, 1)
        clusters <- sort(unique(UniqueID_Table_Unique$Cluster))
        # print(clusters)
        if(all(clusterIDs %in% clusters)){
          cat("All ", length(clusterIDs), " cluster are matched with provided q-value cut-off: ", minQVal, "\n")
          UniqueID_Table_Filtered = UniqueID_Table_Unique %>%
            filter(Cluster %in% unique(clusterIDs)[clusterIDs %in% clusters]) %>%
            rename(PeptideID = Full, UniqueID = Short)
          return(list("MatchedClusters" = clusterIDs,
                      "UnMatchedClusters" = NULL,
                      "UniqueID_Table_Unique" = UniqueID_Table_Filtered,
                      "heatmapType" = "UniqueSet"))
        } else if (any(clusterIDs %in% clusters)){
          cat(paste0("Only ", length(unique(clusterIDs)[clusterIDs %in% clusters]),
                     " out of ", length(clusterIDs), " cluster are matched with provided q-value cut-off: ", minQVal, "\n"))
          UniqueID_Table_Filtered = UniqueID_Table_Unique %>%
            filter(Cluster %in% unique(clusterIDs)[clusterIDs %in% clusters]) %>%
            rename(PeptideID = Full, UniqueID = Short)
          return(list("MatchedClusters" = unique(clusterIDs)[clusterIDs %in% clusters],
                      "UnMatchedClusters" = unique(clusterIDs)[!(clusterIDs %in% clusters)],
                      "UniqueID_Table_Unique" = UniqueID_Table_Filtered,
                      "heatmapType" = "UniqueSet"))
        } else {
          if(verbose==TRUE){
            cat("The avaliable Clusters are: \n")
            cat(paste(clusters, collapse = ", "), "\n")
            cat("With a total of", nrow(UniqueID_Table_Unique), " peptides. \n")
            fullTable = UniqueID_Table_Unique %>% rename(PeptideID = Full, UniqueID = Short) %>%
              group_by(Cluster, UniqueID) %>% arrange(Cluster) %>%
              as_tibble()
            return(fullTable)
          } else{
            cat("No cluster was provided, returning a list of avaliable clusters if verbose = TRUE")
          }
        }
      }
    }

  } else {
    print("NetPhorce data was incomplete. Please process the data again")
  }
}

#' Plot Heatmap Based on the Unique IDs Verified from the \code{\link{findUniqueIDs}} function
#'
#' @description generate Heatmap Based on the Unique IDs matched to the netPhorce Data from findUniqueIDs() function
#' @param netPhorceData (Required). Processed netPhorce Object
#' @param foundUniqueIDs (Required). Output from the \code{\link{findUniqueIDs}} function. If nothing in the “MatchedIDs” vector, no plot will be generated
#' @param heatmapType (Required). Either “Significant” or “AbsencePresence” are required. “Significant” for the unique peptide IDs that were significant differentially phosphorylated peptides. “AbsencePresence” for the unique peptide IDs that were absent at at least one time point.
#' @param minQVal (Required). A q-value threshold. A lower threshold increases stringency for displaying significant phosphopeptides into the heatmap. The q-value controls the positive false discovery rate and is estimated based on the p-values obtained from a linear mixed model fitted to the intensities of each peptide.
#' @param plotly (Required). If TRUE, output an interactive \code{\link{plotly}} version, else output a static \code{\link{ggplot2}} version.
#' @return a ggplot/plotly object
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Provided IDs all matched the IDs in the netPhorce Data
#' uniqueIDs <-
#'   findUniqueIDs(netPhorceData = netPhorceData,
#'                 uniqueIDList = c(
#'                   # Significant Set Examples:
#'                   "AT1G13030.1", "AT1G13360.3", "AT1G42550.1",
#'                   # Unique/Abscence Set Examples
#'                   "AT1G17280.9", "AT1G22310.2", "AT1G23890.2"
#'                 ),
#'                 verbose = TRUE)
#' ## Plotting Significant Heatmap - PLOTLY Version
#' plotUniqueIDsHeatmaps(netPhorceData = netPhorceData,
#'                       foundUniqueIDs = uniqueIDs,
#'                       heatmapType = "Significant",
#'                       minQVal = 0.1,
#'                       plotly = TRUE)
#' ## Plotting Absence/Prescence Heatmap - GGPLOT version
#' plotUniqueIDsHeatmaps(netPhorceData = netPhorceData,
#'                       foundUniqueIDs = uniqueIDs,
#'                       heatmapType = "AbsencePresence",
#'                       plotly = FALSE)
#' }
plotUniqueIDsHeatmaps <- function(netPhorceData = netPhorceData,
                                  foundUniqueIDs = NULL,
                                  heatmapType = "Significant",
                                  minQVal = 0.05,
                                  plotly = TRUE){
  if(is.null(foundUniqueIDs)){
    stop("Please run findUniqueIDs() to identify the peptides you want for plotting. ",
         call. = FALSE)
  }

  if(!is.null(netPhorceData@data.filtered)){
    data.filtered = netPhorceData@data.filtered
    Threshold = netPhorceData@Design$minReplication
    data.filtered.aov.summary = netPhorceData@data.filtered.aov.summary
    uniqueTimes = netPhorceData@Design$uniqueTimes
    colors_HeatmapV1 = netPhorceData@Misc$colors_HeatmapV1
    colors_HeatmapV2 = netPhorceData@Misc$colors_HeatmapV2

    data.clustered.avg = data.filtered %>%
      group_by(set, UniqueID, timepoint, condition) %>% filter(obs_tp >= Threshold) %>%
      summarize(m = mean(normValue, na.rm=TRUE)) %>% ungroup() %>%
      mutate(timepoint = paste(timepoint, condition, sep = "_"), .keep = "unused") %>%
      spread(timepoint, m, fill=0)

    if(length(heatmapType) > 1){
      stop("Please choose a single type for heatmap from 'AbsencePresence' or 'Significant'", call. = FALSE)
    }
    if(!is.null(foundUniqueIDs$UnMatchedIDs)){
      cat("Warning: Not all IDs provided are matched, ", length(foundUniqueIDs$UnMatchedIDs),
          " IDs are omitted from plotting\n")
      selectedIDs = foundUniqueIDs$MatchedIDs
    } else {
      selectedIDs = foundUniqueIDs$MatchedIDs
    }
    # print(selectedIDs)

    if(any(heatmapType %in% c("AbsencePresence", "Significant"))){
      if(heatmapType == "AbsencePresence"){
        curSetName = "UniqueSet"
      } else if(heatmapType == "Significant"){
        curSetName = "StatsSet"
      }
    } else {
      stop("Please choose a single type for heatmap from 'AbsencePresence' or 'Significant'", call. = FALSE)
    }

    if(curSetName == "StatsSet"){
      if(!any(data.clustered.avg$set == curSetName)){
        stop("No significant peptides are found. ", call. = FALSE)
      } else {
        heatmap_Significant_phosphosites.summary = data.clustered.avg %>%
          inner_join(data.filtered.aov.summary %>% filter(qvalue < minQVal) %>% dplyr::select(UniqueID)) %>%
          dplyr::select(-set) %>% distinct() %>%
          nest(data = everything()) %>% ungroup() %>%
          mutate(dist = map(data, ~ .x %>% column_to_rownames("UniqueID") %>% dist())) %>%
          mutate(hclust = map(dist, hclust,  method = "complete")) %>%
          mutate(order = map(hclust,"order")) %>%
          mutate(cutreeDynamic = map2(hclust, dist, ~ cutreeDynamic(dendro = .x, distM = as.matrix(.y),
                                                                    method = "hybrid", deepSplit = 4)))
        UniqueID_Table <- data.frame(Full = heatmap_Significant_phosphosites.summary$data[[1]]$UniqueID)
        UniqueID_Table$Cluster <- heatmap_Significant_phosphosites.summary$cutreeDynamic[[1]]
        UniqueID_Table$Short <- sapply(strsplit(as.character(UniqueID_Table$Full), "_"), `[`, 1)
        IDs <- as.vector(subset(UniqueID_Table, Short %in% selectedIDs)$Full)
        clusters <- unique(UniqueID_Table$Cluster)
        # print(UniqueID_Table)

        if(length(IDs) == 0){
          stop(paste0("No significant peptides are found at the maximum allowed q value cut-off at ",
                      minQVal, " with the selected Unique IDs"), call. = FALSE)
        } else {
          # newClusters <- as.vector(subset(UniqueID_Table, Short %in% IDs)$Cluster)

          # return(heatmap_Significant_phosphosites.summary)
          heatmap_temp <-
            heatmap_Significant_phosphosites.summary %>%
            dplyr::select(data, order, cutreeDynamic) %>% unnest(cols = everything()) %>%
            filter(cutreeDynamic %in% clusters) %>%
            gather(experiment,value, -UniqueID, -order, -cutreeDynamic) %>%
            mutate(timepoint = stringr::str_split(experiment, "_") %>% map_chr(., 1),
                   condition = stringr::str_split(experiment, "_") %>% map_chr(., 2)) %>%
            # filter(value != 0 & UniqueID %in% IDs)  %>%
            mutate(value = replace(value, value==0, NA)) %>%
            filter(UniqueID %in% IDs) %>%
            mutate_at(.vars = vars(dplyr::matches("timepoint")),
                      .funs = funs(factor(., levels = uniqueTimes$Origin))) %>%
            ggplot(., aes(x = timepoint, y = UniqueID, fill = value)) + geom_tile() +
            scale_fill_gradientn(colors = colors_HeatmapV1, na.value = "#FFFFFF") +
            ggtitle(paste0("Significant Phosphopeptides (q-value < ", minQVal, ")")) +
            facet_wrap(cutreeDynamic ~ condition , scales = "free_y", dir = "v", ncol=1, strip.position = "right") +
            scale_x_discrete(limits = as.vector(uniqueTimes$Origin)) +
            theme(axis.ticks = element_blank(),
                  axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 8),
                  # axis.text.y = element_blank(),
                  legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"),
                  panel.background = element_rect(fill = "white",
                                                  color = "white", size = 0.5, linetype = "solid"),
                  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                  colour = "gray"),
                  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                  colour = "gray"),
                  plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                  strip.text.x = element_text(size = 6, face = "bold"),
                  strip.text.y = element_text(size = 10, face = "bold"))
          if(plotly == TRUE){
            heatmap_plotly <- heatmap_temp +
              labs(x = "Time Points", y = "", fill = "Average \nLog<sub>2</sub>\nIntensity")
            heatmapPlotly <- ggplotly(heatmap_plotly)
            return(heatmapPlotly)
          } else {
            heatmap <- heatmap_temp +
              labs(x = "Time Points", y = "", fill = expression(paste("Average "*Log[2]*" Intensity")))
            return(heatmap)
          }
        }
      }
    }


    if(curSetName == "UniqueSet"){
      if(!any(data.clustered.avg$set == curSetName)){
        stop("No absence/presence peptides are found. ", call. = FALSE)
      } else {
        heatmap_Unique_phosphosites.summary = data.clustered.avg %>% filter(set == "UniqueSet") %>%
          dplyr::select(-set) %>% nest(data = everything()) %>% ungroup() %>%
          distinct() %>%
          mutate(dist = map(data, ~ .x %>% column_to_rownames("UniqueID") %>% dist())) %>%
          mutate(hclust = map(dist, hclust,  method = "complete")) %>%
          mutate(order = map(hclust,"order")) %>%
          mutate(cutreeDynamic = map2(hclust, dist, ~ cutreeDynamic(dendro = .x, distM = as.matrix(.y),
                                                                    method = "hybrid", deepSplit = 4)))

        UniqueID_Table_Unique <- data.frame(Full = heatmap_Unique_phosphosites.summary$data[[1]]$UniqueID)
        UniqueID_Table_Unique$Cluster <- heatmap_Unique_phosphosites.summary$cutreeDynamic[[1]]
        UniqueID_Table_Unique$Short <- sapply(strsplit(as.character(UniqueID_Table_Unique$Full), "_"), `[`, 1)
        IDs <- as.vector(subset(UniqueID_Table_Unique, Short %in% selectedIDs)$Full)
        clusters <- unique(UniqueID_Table_Unique$Cluster)

        if(length(IDs) == 0){
          stop(paste0("No absence/presence peptides are found with the selected Unique IDs"), call. = FALSE)
        } else {
          # clusters <- unique(heatmap_Significant_phosphosites.summary$cutreeDynamic[[1]])

          heatmap_temp <-
            heatmap_Unique_phosphosites.summary %>%
            dplyr::select(data, order, cutreeDynamic) %>% unnest(cols = everything()) %>%
            filter(cutreeDynamic %in% clusters) %>%
            gather(experiment,value, -UniqueID, -order, -cutreeDynamic) %>%
            mutate(timepoint = stringr::str_split(experiment, "_") %>% map_chr(., 1),
                   condition = stringr::str_split(experiment, "_") %>% map_chr(., 2)) %>%
            # filter(value != 0 & UniqueID %in% IDs)  %>%
            mutate(value = replace(value, value==0, NA)) %>%
            filter(UniqueID %in% IDs) %>%
            mutate_at(.vars = vars(dplyr::matches("timepoint")),
                      .funs = funs(factor(., levels = uniqueTimes$Origin))) %>%
            ggplot(., aes(x = timepoint, y = UniqueID, fill = value)) + geom_tile() +
            scale_fill_gradientn(colors = colors_HeatmapV2, na.value = "#FFFFFF") +
            ggtitle("Absence/presence Phosphopeptides") +
            facet_wrap(cutreeDynamic ~ condition , scales = "free_y", dir = "v", ncol=1, strip.position = "right") +
            scale_x_discrete(limits = as.vector(uniqueTimes$Origin)) +
            theme(axis.ticks = element_blank(),
                  axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 8),
                  # axis.text.y = element_blank(),
                  legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"),
                  panel.background = element_rect(fill = "white",
                                                  color = "white", size = 0.5, linetype = "solid"),
                  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                  colour = "gray"),
                  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                  colour = "gray"),
                  plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                  strip.text.x = element_text(size = 6, face = "bold"),
                  strip.text.y = element_text(size = 10, face = "bold"))
          if(plotly == TRUE){
            heatmap_plotly <- heatmap_temp +
              labs(x = "Time Points", y = "", fill = "Average \nLog<sub>2</sub>\nIntensity")
            heatmapPlotly <- ggplotly(heatmap_plotly)
            return(heatmapPlotly)
          } else {
            heatmap <- heatmap_temp +
              labs(x = "Time Points", y = "", fill = expression(paste("Average "*Log[2]*" Intensity")))
            return(heatmap)
          }
        }
      }
    }
  } else {
    print("NetPhorce data was incomplete. Please process the data again")
  }
}

#' Plot Heatmap Based on the Clusters Verified from the \code{\link{findClusters}} function
#'
#' @description generate Heatmap Based on the clusters matched to the netPhorce Data from findUniqueIDs() function
#' @param netPhorceData (Required). Processed netPhorce Object
#' @param foundClusterIDs (Required). Output from the \code{\link{findClusters}} function. If nothing in the “MatchedClusters” vector, no plot will be generated.
#' @param plotly (Required). If TRUE, output an interactive \code{\link{plotly}} version, else output a static \code{\link{ggplot2}} version.
#' @return a ggplot/plotly object
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Returns all available peptides along with IDs from Significant Set
#' clusterIDs_Sig <-
#'   findClusters(netPhorceData = netPhorceData,
#'                clusterIDs = c(0),
#'                heatmapType = "Significant",
#'                minQVal = 0.1,
#'                verbose = TRUE)
#' ## Returns all available peptides Identify Cluster IDs from “AbsencePresence” Set
#' clusterIDs_AbsPrs <-
#'   findClusters(netPhorceData = netPhorceData,
#'                clusterIDs = c(1, 2),
#'                heatmapType = "AbsencePresence",
#'                verbose = TRUE)
#' ## Plotting Significant Heatmap - PLOTLY Version
#' plotClustersHeatmap(netPhorceData = netPhorceData,
#'                     foundClusterIDs = clusterIDs_Sig,
#'                     plotly = TRUE)
#' ## Plotting Significant Heatmap - GGPLOT version
#' plotClustersHeatmap(netPhorceData = netPhorceData,
#'                     foundClusterIDs = clusterIDs_AbsPrs,
#'                     plotly = FALSE)
#' }
plotClustersHeatmap <- function(netPhorceData = netPhorceData,
                                foundClusterIDs = NULL,
                                plotly = TRUE){
  if(is.null(foundClusterIDs)){
    stop("Please run findClusters to identify the clusters you want for plotting. ",
         call. = FALSE)
  }
  if(is.null(suppressWarnings(foundClusterIDs$heatmapType))){
    stop("Please select the correct clusters using findClusters() function's clusterID argument.")
  }


  if(!is.null(netPhorceData@data.filtered)){
    data.filtered = netPhorceData@data.filtered
    Threshold = netPhorceData@Design$minReplication
    data.filtered.aov.summary = netPhorceData@data.filtered.aov.summary
    uniqueTimes = netPhorceData@Design$uniqueTimes
    colors_HeatmapV1 = netPhorceData@Misc$colors_HeatmapV1
    colors_HeatmapV2 = netPhorceData@Misc$colors_HeatmapV2


    data.clustered.avg = data.filtered %>%
      group_by(set, UniqueID, timepoint, condition) %>% filter(obs_tp >= Threshold) %>%
      summarize(m = mean(normValue, na.rm=TRUE)) %>% ungroup() %>%
      mutate(timepoint = paste(timepoint, condition, sep = "_"), .keep = "unused") %>%
      spread(timepoint, m, fill=0)

    if(foundClusterIDs$heatmapType == "StatsSet"){
        heatmap_Significant_phosphosites.summary = data.clustered.avg %>%
          inner_join(data.filtered.aov.summary %>% filter(qvalue < foundClusterIDs$minQVal) %>% dplyr::select(UniqueID)) %>%
          dplyr::select(-set) %>% distinct() %>%
          nest(data = everything()) %>% ungroup() %>%
          mutate(dist = map(data, ~ .x %>% column_to_rownames("UniqueID") %>% dist())) %>%
          mutate(hclust = map(dist, hclust,  method = "complete")) %>%
          mutate(order = map(hclust,"order")) %>%
          mutate(cutreeDynamic = map2(hclust, dist, ~ cutreeDynamic(dendro = .x, distM = as.matrix(.y),
                                                                    method = "hybrid", deepSplit = 4)))
        UniqueID_Table <- data.frame(Full = heatmap_Significant_phosphosites.summary$data[[1]]$UniqueID)
        UniqueID_Table$Cluster <- heatmap_Significant_phosphosites.summary$cutreeDynamic[[1]]
        UniqueID_Table$Short <- sapply(strsplit(as.character(UniqueID_Table$Full), "_"), `[`, 1)
        IDs <- UniqueID_Table$Full
        clusters <- foundClusterIDs$MatchedClusters
        # print(UniqueID_Table)
        # return(heatmap_Significant_phosphosites.summary)
        heatmap_temp <-
          heatmap_Significant_phosphosites.summary %>%
          dplyr::select(data, order, cutreeDynamic) %>% unnest(cols = everything()) %>%
          filter(cutreeDynamic %in% clusters) %>%
          gather(experiment,value, -UniqueID, -order, -cutreeDynamic) %>%
          mutate(timepoint = stringr::str_split(experiment, "_") %>% map_chr(., 1),
                 condition = stringr::str_split(experiment, "_") %>% map_chr(., 2)) %>%
          # filter(value != 0 & UniqueID %in% IDs)  %>%
          mutate(value = replace(value, value==0, NA)) %>%
          filter(UniqueID %in% IDs) %>%
          mutate_at(.vars = vars(dplyr::matches("timepoint")),
                    .funs = funs(factor(., levels = uniqueTimes$Origin))) %>%
          ggplot(., aes(x = timepoint, y = UniqueID, fill = value)) + geom_tile() +
          scale_fill_gradientn(colors = colors_HeatmapV1, na.value = "#FFFFFF") +
          ggtitle(paste0("Significant Phosphopeptides (q-value < ", foundClusterIDs$minQVal, ")")) +
          facet_wrap(cutreeDynamic ~ condition , scales = "free_y", dir = "v", ncol=1, strip.position = "right") +
          scale_x_discrete(limits = as.vector(uniqueTimes$Origin)) +
          theme(axis.ticks = element_blank(),
                axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 8),
                # axis.text.y = element_blank(),
                legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"),
                panel.background = element_rect(fill = "white",
                                                color = "white", size = 0.5, linetype = "solid"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                colour = "gray"),
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                colour = "gray"),
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                strip.text.x = element_text(size = 6, face = "bold"),
                strip.text.y = element_text(size = 10, face = "bold"))
        if(plotly == TRUE){
          heatmap_plotly <- heatmap_temp +
            labs(x = "Time Points", y = "", fill = "Average \nLog<sub>2</sub>\nIntensity")
          heatmapPlotly <- ggplotly(heatmap_plotly)
          return(heatmapPlotly)
        } else {
          heatmap <- heatmap_temp +
            labs(x = "Time Points", y = "", fill = expression(paste("Average "*Log[2]*" Intensity")))
          return(heatmap)
        }

    }

    if(foundClusterIDs$heatmapType == "UniqueSet"){
      heatmap_Unique_phosphosites.summary = data.clustered.avg %>% filter(set == "UniqueSet") %>%
        dplyr::select(-set) %>% nest(data = everything()) %>% ungroup() %>%
        distinct() %>%
        mutate(dist = map(data, ~ .x %>% column_to_rownames("UniqueID") %>% dist())) %>%
        mutate(hclust = map(dist, hclust,  method = "complete")) %>%
        mutate(order = map(hclust,"order")) %>%
        mutate(cutreeDynamic = map2(hclust, dist, ~ cutreeDynamic(dendro = .x, distM = as.matrix(.y),
                                                                  method = "hybrid", deepSplit = 4)))

      UniqueID_Table_Unique <- data.frame(Full = heatmap_Unique_phosphosites.summary$data[[1]]$UniqueID)
      UniqueID_Table_Unique$Cluster <- heatmap_Unique_phosphosites.summary$cutreeDynamic[[1]]
      UniqueID_Table_Unique$Short <- sapply(strsplit(as.character(UniqueID_Table_Unique$Full), "_"), `[`, 1)
      clusters <-foundClusterIDs$MatchedClusters
      IDs <- UniqueID_Table_Unique$Full

      heatmap_temp <-
        heatmap_Unique_phosphosites.summary %>%
        dplyr::select(data, order, cutreeDynamic) %>% unnest(cols = everything()) %>%
        filter(cutreeDynamic %in% clusters) %>%
        gather(experiment,value, -UniqueID, -order, -cutreeDynamic) %>%
        mutate(timepoint = stringr::str_split(experiment, "_") %>% map_chr(., 1),
               condition = stringr::str_split(experiment, "_") %>% map_chr(., 2)) %>%
        # filter(value != 0 & UniqueID %in% IDs)  %>%
        mutate(value = replace(value, value==0, NA)) %>%
        filter(UniqueID %in% IDs) %>%
        mutate_at(.vars = vars(dplyr::matches("timepoint")),
                  .funs = funs(factor(., levels = uniqueTimes$Origin))) %>%
        ggplot(., aes(x = timepoint, y = UniqueID, fill = value)) + geom_tile() +
        scale_fill_gradientn(colors = colors_HeatmapV2, na.value = "#FFFFFF") +
        ggtitle("Absence/presence Phosphopeptides") +
        facet_wrap(cutreeDynamic ~ condition , scales = "free_y", dir = "v", ncol=1, strip.position = "right") +
        scale_x_discrete(limits = as.vector(uniqueTimes$Origin)) +
        theme(axis.ticks = element_blank(),
              axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 8),
              # axis.text.y = element_blank(),
              legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"),
              panel.background = element_rect(fill = "white",
                                              color = "white", size = 0.5, linetype = "solid"),
              panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                              colour = "gray"),
              panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                              colour = "gray"),
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
              strip.text.x = element_text(size = 6, face = "bold"),
              strip.text.y = element_text(size = 10, face = "bold"))
      if(plotly == TRUE){
        heatmap_plotly <- heatmap_temp +
          labs(x = "Time Points", y = "", fill = "Average \nLog<sub>2</sub>\nIntensity")
        heatmapPlotly <- ggplotly(heatmap_plotly)
        return(heatmapPlotly)
      } else {
        heatmap <- heatmap_temp +
          labs(x = "Time Points", y = "", fill = expression(paste("Average "*Log[2]*" Intensity")))
        return(heatmap)
      }
    }
  } else {
    print("NetPhorce data was incomplete. Please process the data again")
  }
}


#' Identify provided peptide IDs for plotting
#'
#' @description Return the input ID when avaiable in the dataset or output all available IDs
#' @param netPhorceData Processed NetPhorce Object
#' @param peptideIDList A single or multiple peptide IDs. The peptide ID consists of a combination of the protein ID, amino acid, and multiplicity. If left as `NULL`, the full table of available IDs will be printed.
#' @return list of data.frames
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Extract peptides from the NetPhorce Data
#' findPeptideIDs(netPhorceData = netPhorceData,
#'                peptideIDList = c("AT1G01320.2_S1349_1",
#'                                  "AT1G05560.1_S153_1",
#'                                  "AT1G01320.2_S149_1"))
#' }
findPeptideIDs <- function(netPhorceData = netPhorceData,
                           peptideIDList = NULL){
  if(!is.null(netPhorceData@data.filtered)){
    data.filtered = netPhorceData@data.filtered
    data.filtered_UniqueID_Table <- data.frame(PeptideID = data.filtered$UniqueID)
    data.filtered_UniqueID_Table$Protein = data.filtered_UniqueID_Table$PeptideID
    data.filtered_UniqueID_Table <- data.filtered_UniqueID_Table %>%
      separate(Protein, into = c("Part1", "multiplicity"), sep = "_(?=[0-9]+$)") %>%
      separate(Part1, into = c("UniqueID", "AA"), sep = "_(?=[0-9A-z;]+$)")
    # data.filtered_UniqueID_Table$Short <- sapply(strsplit(as.character(data.filtered_UniqueID_Table$Full), "_"), `[`, 1)
    # print(head(data.filtered_UniqueID_Table))
    if(is.null(peptideIDList)){
      cat("No phosphopeptide id was provided, returning a table of avaliable unique IDs. Please use the ID under the column 'PeptideID'\n")
      fulltable = data.filtered_UniqueID_Table %>% select(PeptideID, UniqueID) %>% distinct() %>% as_tibble()
      return(fulltable)
    } else {
      if(all(peptideIDList %in% unique(data.filtered_UniqueID_Table$PeptideID))){
        cat("All ", length(peptideIDList), " provided phosphopeptide IDs are matched to existing phosphopeptide IDs\n")
        return(list("MatchedIDs" = peptideIDList,
                    "UnMatchedIDs" = NULL))
      } else if (any(peptideIDList %in% unique(data.filtered_UniqueID_Table$PeptideID))){
        cat(paste0("Only ", length(unique(peptideIDList)[peptideIDList %in% unique(data.filtered_UniqueID_Table$PeptideID)]),
                   " out of ", length(peptideIDList), " provided IDs matched to existing phosphopeptide IDs\n"))
        return(list("MatchedIDs" = unique(peptideIDList)[peptideIDList %in% unique(data.filtered_UniqueID_Table$PeptideID)],
                    "UnMatchedIDs" = unique(peptideIDList)[!(peptideIDList %in% unique(data.filtered_UniqueID_Table$PeptideID))]))
      } else {
        cat("No phosphopeptide id was matched, to return a table with avaliable peptide IDs if peptideIDList = NULL")
      }
    }
  } else {
    print("NetPhorce data was incomplete. Please process the data again")
  }
}

#' Dotplot displaying the log2 intensity of a single phosphopeptide across the different time points and for each condition.
#'
#' @description Return if the input unque ID is in the dataset and filter out current results related to this
#' @param netPhorceData (Required). Processed NetPhorce Object
#' @param foundPepetidesIDs (Required). Single or Multiple Peptides ID, you can use \code{foundPepetidesIDs} with ` peptideIDList = NULL to extract a full table of avaliable peptides.
#' @param plotAll (Optional). If TRUE, all matched peptides provided will be plotted in their individual plots.
#' @param plotly (Required). If TRUE, output an interactive \code{\link{plotly}} version, else output a static \code{\link{ggplot2}} version.
#' @return one or multiple ggplot/plotly object(s)
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Provided IDs all matched the IDs in the netPhorce Data
#' peptideIDs <-
#'   findPeptideIDs(netPhorceData = netPhorceData,
#'                  peptideIDList = c("AT1G01320.2_S1349_1",
#'                                    "AT1G05560.1_S153_1",
#'                                    "AT1G01320.2_S149_1"))
#' ## Plot the first Phosphopeptide ID - GGPLOT (Static) Version
#' plotSinglePeptide(netPhorceData = netPhorceData,
#'                    foundPepetidesIDs = peptideIDs,
#'                    plotAll = FALSE,
#'                    plotly = TRUE)
#' ## Plot All Provided Phosphopeptide IDs - PLOTLY (Interactive) Version
#' plotSinglePeptide(netPhorceData = netPhorceData,
#'                    foundPepetidesIDs = peptideIDs,
#'                    plotAll = TRUE,
#'                    plotly = FALSE)
#' }
plotSinglePeptide <- function(netPhorceData = netPhorceData,
                              foundPepetidesIDs = NULL,
                              plotAll = FALSE,
                              plotly = FALSE){
  plotData <- netPhorceData@data.filtered
  # FASTA_HeaderKey <- netPhorceData@FASTA_HeaderKey
  FASTA_HeaderKey = NULL
  colors_V2 <- netPhorceData@Misc$colors_V2

  if(is.null(foundPepetidesIDs)){
    stop("Please run findPeptideIDs() function to identify peptides first. ",
         call. = FALSE)
  }

  if(plotAll == TRUE){
    returnList = c()
    for(curID in foundPepetidesIDs$MatchedIDs){
      # print(curID)
      plot1Data <- subset(plotData, UniqueID == curID)
      plot1Data$Time <- as.numeric(str_extract(as.character(plot1Data$timepoint), "[0-9]+"))
      # print(head(plot1Data))
      ## Add FASTA Name for the hover effect
      if(!is.null(FASTA_HeaderKey)){
        curFASTA <- curID
        curFASTA <- paste0(sapply(strsplit(curFASTA, "_"), "[", 1), "_",
                           sapply(strsplit(curFASTA, "_"), "[", 2))
        txtListNew = FASTA_HeaderKey$`FASTA Header`[match(x = curFASTA,
                                                          table = FASTA_HeaderKey$UniqueID)]
        txtListNew <- paste0(strwrap(txtListNew, width = 50), collapse = "<br>")
      } else {
        txtListNew = curID
      }

      ## Extract the FASTA Header if avaliable
      plot1Data$`Time Point` <- plot1Data$timepoint
      plot1Data$`Normalized Value` <- plot1Data$normValue
      plot1Data$Description <- txtListNew
      plot1Data$Replication <- plot1Data$replicate
      ## Base Plot
      plot1_ggplot <-
        ggplot(data = plot1Data, aes(
          x = Time,
          y = `Normalized Value`, label = Description)) +
        geom_point(aes(color = Replication),
                   position = position_jitter(width = 0.15), size = 2, alpha = 0.8) +
        stat_smooth(data = plot1Data %>%
                      group_by(condition, timepoint) %>%
                      filter(n() > 2),
                    method = "gam",
                    formula = y ~ s(x, k=3),
                    se = TRUE) +
        scale_x_continuous(breaks = as.vector(unique(plot1Data$Time)),
                           labels = as.vector(unique(plot1Data$`Time Point`))) +
        facet_grid(condition ~ ., scales = "fixed") +
        theme(plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = "white",
                                              colour = "white",
                                              size = 0.5, linetype = "solid"),
              panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                              colour = "#D1D1D1"),
              panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                              colour = "#D1D1D1")) +
        scale_color_manual(values = colors_V2)

      if(plotly == TRUE){
        plot1_ggplotly <- plot1_ggplot +
          labs(x = " ",
               y = " ",
               # title = " "
               title = " "
          ) +
          theme(legend.title = element_blank())

        plot1_ggplotly <-
          ggplotly(plot1_ggplotly)
        annotation <- list(yref = 'paper', xref = "paper", y = 0.5, x = -0.075,
                           text = "<b>Log<sub>2</sub> Intensity</b>",
                           showarrow = FALSE,
                           textangle = 270)

        plot1_ggplotly = plot1_ggplotly %>%
          layout(
            legend = list(
              x = 1.1, y = 0.5,
              title = list(text = "<b>Replicate</b>")
            ),
            title = curID,
            xaxis = list(title = "<br>Time Points</b>"),
            annotations= list(annotation)
          )
        returnList[[curID]] <- (plot1_ggplotly)
      } else {
        plot1_ggplot = plot1_ggplot +
          labs(x = "Time Points",
               y = expression(paste(Log[2]*" Intensity")),
               title = curID)
        returnList[[curID]] <- (plot1_ggplot)
      }
    }
    return(returnList)
  } else {
    cat("Multiple Phosphopeptide IDs provided, only plotting the first one. Please
        use plotAll if you want to plot all the found phosphopeptide IDS. \n")
    plot1Data <- subset(plotData, UniqueID == foundPepetidesIDs$MatchedIDs[1])
    plot1Data$Time <- as.numeric(str_extract(as.character(plot1Data$timepoint), "[0-9]+"))

    ## Add FASTA Name for the hover effect
    if(!is.null(FASTA_HeaderKey)){
      curFASTA <- foundPepetidesIDs$MatchedIDs[1]
      curFASTA <- paste0(sapply(strsplit(curFASTA, "_"), "[", 1), "_",
                         sapply(strsplit(curFASTA, "_"), "[", 2))
      txtListNew = FASTA_HeaderKey$`FASTA Header`[match(x = curFASTA,
                                                        table = FASTA_HeaderKey$UniqueID)]
      txtListNew <- paste0(strwrap(txtListNew, width = 50), collapse = "<br>")
    } else {
      txtListNew = foundPepetidesIDs$MatchedIDs[1]
    }

    ## Extract the FASTA Header if avaliable
    plot1Data$`Time Point` <- plot1Data$timepoint
    plot1Data$`Normalized Value` <- plot1Data$normValue
    plot1Data$Description <- txtListNew
    plot1Data$Replication <- plot1Data$replicate

    ## Base Plot
    plot1_ggplot <-
      ggplot(data = plot1Data, aes(
        x = Time,
        y = `Normalized Value`, label = Description)) +
      geom_point(aes(color = Replication),
                 position = position_jitter(width = 0.15), size = 2, alpha = 0.8) +
      stat_smooth(data = plot1Data %>%
                    group_by(condition, timepoint) %>%
                    filter(n() > 2),
                  method = "gam",
                  formula = y ~ s(x, k=3),
                  se = TRUE) +
      scale_x_continuous(breaks = as.vector(unique(plot1Data$Time)),
                         labels = as.vector(unique(plot1Data$`Time Point`))) +
      facet_grid(condition ~ ., scales = "fixed") +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "#D1D1D1"),
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "#D1D1D1")) +
      scale_color_manual(values = colors_V2)

    if(plotly == TRUE){
      plot1_ggplotly <- plot1_ggplot +
        labs(x = " ",
             y = " ",
             # title = " "
             title = " "
        ) +
        theme(legend.title = element_blank())

      plot1_ggplotly <-
        ggplotly(plot1_ggplotly)
      annotation <- list(yref = 'paper', xref = "paper", y = 0.5, x = -0.075,
                         text = "<b>Log<sub>2</sub> Intensity</b>",
                         showarrow = FALSE,
                         textangle = 270)

      plot1_ggplotly = plot1_ggplotly %>%
        layout(
          legend = list(
            x = 1.1, y = 0.5,
            title = list(text = "<b>Replicate</b>")
          ),
          title = foundPepetidesIDs$MatchedIDs[1],
          xaxis = list(title = "<br>Time Points</b>"),
          annotations= list(annotation)
        )
      return(plot1_ggplotly)
    } else {
      plot1_ggplot = plot1_ggplot +
        labs(x = "Time Points",
             y = expression(paste(Log[2]*" Intensity")),
             title = foundPepetidesIDs$MatchedIDs[1])
      return(plot1_ggplot)
    }
  }

}

#' Intensity pattern of a single or multiple phosphopeptide across the time course and for each condition.
#'
#' @description Return the intensity pattern across multiple time points for the provided phosphopeptide IDs.
#' @param netPhorceData (Required). Processed NetPhorce Object
#' @param condition (Required). Select a specific condition from your experiment.
#' @param foundPepetidesIDs (Required). Single or multiple peptide IDs. Use \code{foundPepetidesIDs} with ` peptideIDList = NULL ` to extract a full table of available peptide IDs.
#' @param plotly (Required). If TRUE, output an interactive \code{\link{plotly}} version, else output a static \code{\link{ggplot2}} version.
#' @return list of dataframes
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
#' ## Provided IDs all matched the IDs in the netPhorce Data
#' peptideIDs <-
#'   findPeptideIDs(netPhorceData = netPhorceData,
#'                  peptideIDList = c("AT1G01320.2_S1349_1",
#'                                    "AT1G05560.1_S153_1",
#'                                    "AT1G01320.2_S149_1"))
#' ## Plot Mutliple Phosphopeptide IDs - PLOTLY (Interactive) Version
#' plotMultiPeptides(netPhorceData = netPhorceData,
#'                   foundPepetidesIDs = peptideIDs,
#'                   condition = "Col0",
#'                   plotly = TRUE)
#' ## Plot Mutliple Phosphopeptide IDs - GGPLOT (Static) Version
#' plotMultiPeptides(netPhorceData = netPhorceData,
#'                   foundPepetidesIDs = peptideIDs,
#'                   condition = "tot3",
#'                   plotly = FALSE)
#' }
plotMultiPeptides <- function(netPhorceData = netPhorceData,
                              foundPepetidesIDs = NULL,
                              condition = NULL,
                              plotly = TRUE){
  plotData <- netPhorceData@data.filtered
  # FASTA_HeaderKey <- netPhorceData@FASTA_HeaderKey
  FASTA_HeaderKey = NULL
  colors <- netPhorceData@Misc$colors

  if(is.null(foundPepetidesIDs)){
    stop("Please run findPeptideIDs() function to identify peptides first. ",
         call. = FALSE)
  }
  plot1Data <- subset(plotData, UniqueID == foundPepetidesIDs$MatchedIDs[1])
  plot1Data$Time <- as.numeric(str_extract(as.character(plot1Data$timepoint), "[0-9]+"))

  ## Check Condition
  ConDesign <- netPhorceData@Design$ConDesign
  nConditions = n_distinct(ConDesign$condition)
  if(is.null(condition)){
    if(nConditions == 1){
      curCondition = unique(ConDesign$condition)
    } else {
      cat("Multiple condition detected, please use the 'condition = ' arguemnt to select the specific condition.\n")
      cat("The following plot is based on the first found condition which is: \n")
      curCondition = unique(ConDesign$condition)[1]
      cat(curCondition)
      cat("\n")
    }
  } else {
    curCondition = condition
    if(!(curCondition %in% ConDesign$condition)){
      stop("Selected condition is not in the current data, the avaliable conditions are: ",
           paste(unique(ConDesign$condition), collapse = "; "), call. = FALSE)
    }
  }

  ## Add FASTA Name for the hover effect
  if(!is.null(FASTA_HeaderKey)){
    curFASTA <- foundPepetidesIDs$MatchedIDs[1]
    curFASTA <- paste0(sapply(strsplit(curFASTA, "_"), "[", 1), "_",
                       sapply(strsplit(curFASTA, "_"), "[", 2))
    txtListNew = FASTA_HeaderKey$`FASTA Header`[match(x = curFASTA,
                                                      table = FASTA_HeaderKey$UniqueID)]
    txtListNew <- paste0(strwrap(txtListNew, width = 50), collapse = "<br>")
  } else {
    txtListNew = foundPepetidesIDs$MatchedIDs[1]
  }

  plot1Data <- subset(plotData, UniqueID %in% foundPepetidesIDs$MatchedIDs)
  plot1Data <- subset(plot1Data, condition == curCondition)

  plot1Data$`Time Point` <- plot1Data$timepoint
  plot1Data$`Normalized Value` <- plot1Data$normValue
  plot1Data$Description <- txtListNew
  plot1Data$Replication <- plot1Data$replicate
  normalizedData <- plot1Data

  normalizedData <-
    normalizedData %>%
    group_by(UniqueID, timepoint) %>%
    mutate(nTimepointsCount = sum(!is.na(normValue)))
  normalizedData$text <- paste0( "<b>Condition</b>: ", normalizedData$condition, "<br>",
                                 "<b>Time point</b>:   ", normalizedData$timepoint, "<br>",
                                 "<b>Replicate</b>: ", normalizedData$replicate, "<br>",
                                 "<b># Samples per Time point</b>: ", normalizedData$nTimepointsCount)

  colors_darken = darken(colors, 0.2)
  if(plotly == TRUE){
    plot2_ggplotly <-
      plot_ly(data = normalizedData,
              x = ~UniqueID,
              y = ~normValue,
              jitter = 0,
              pointpos = 0,
              color = ~timepoint,
              colors = colors_darken[1:length(unique(normalizedData$timepoint))],
              text  = ~text,
              line = list(color = 'rgb(7,40,89)'),
              marker = list(line = list(color = 'rgb(7,40,89)',
                                        width = 2)
              ),
              boxpoints = 'all',
              type = "box") %>%
      layout(boxmode = "group",
             xaxis = list(title="<b>Samples</b>"),
             yaxis = list(title="<b>Log<sub>2</sub> Intensity</b>"),
             title = paste0("<b>Condition: ", curCondition, "</b>"),
             legend = list(
               # x = 1.1,
               y = 0.5,
               title = list(text = "<b>Time Points</b>")
             ))
    return(plot2_ggplotly)
  } else {
    plot2_ggplot <-
      ggplot(normalizedData, aes(x = UniqueID, y = normValue, fill = timepoint)) +
      geom_boxplot() +
      # geom_point(position=position_dodge(width=0.75), aes(group = Timepoint)) +
      geom_point(position=position_dodge(width=0.75), color = "black", shape = 23) +
      facet_grid(condition ~ .) +
      theme_bw() +
      theme(axis.text.x = element_text(vjust = 0.5),
            legend.title = element_text("Time Points"),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = "Sample ID",
           y = expression(paste("Log"[2]*" Intensity")),
           title = paste0("Condition: ", curCondition),
           fill = "Time Point") +
      scale_fill_manual(values = colors[1:length(unique(normalizedData$timepoint))])
    return(plot2_ggplot)
  }
}


#' The regulation plot allows for the users to visually evaluate the chosen thresholds.
#'
#' @description Return a dot plot where each dot represents the intensity log2 fold change of a time point relative to the immediately preceding time point. This value is plotted onto the median log2 intensity of that time point.
#' @param netPhorceData Processed NetPhorce Object with \code{\link{regulationCheck}} completed.
#' @param condition (Required). Select a specific condition if multiple conditions are present.
#' @param plotly (Required). If TRUE, output an interactive \code{\link{plotly}} version, else output a static \code{\link{ggplot2}} version.
#' @return A ggplot/plotly object
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#'
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
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
#' ## Plot the regulations - GGPLOT (Static) Version
#' plotRegulation(netPhorceData = netPhorceData ,
#'                condition = "Col0",
#'                plotly = FALSE)
#' ## Plot the regulations - PLOTLY (Interactive) Version
#' plotRegulation(netPhorceData = netPhorceData ,
#'                condition = "tot3",
#'                plotly = TRUE)
#' }
plotRegulation <- function(netPhorceData = netPhorceData,
                           condition = NULL,
                           plotly = TRUE){
  data.significant.proteins = netPhorceData@regulationData
  regColors = netPhorceData@Misc$colors_Reg
  filteredData <- data.significant.proteins %>% filter(avgValue != 0) %>% filter(avgValue != change)
  filteredData$sigChangeSign <- as.character(filteredData$sigChangeSign)
  ## Check Condition
  ConDesign <- netPhorceData@Design$ConDesign
  nConditions = n_distinct(ConDesign$condition)
  if(is.null(condition)){
    if(nConditions == 1){
      curCondition = unique(ConDesign$condition)
    } else {
      cat("Multiple condition detected, please use the 'condition = ' arguemnt to select the specific condition.\n")
      cat("The following plot is based on the first found condition which is: \n")
      curCondition = unique(ConDesign$condition)[1]
      cat(curCondition)
      cat("\n")
    }
  } else {
    curCondition = condition
    if(!(curCondition %in% ConDesign$condition)){
      stop("Selected condition is not in the current data, the avaliable conditions are: ",
           paste(unique(ConDesign$condition), collapse = "; "), call. = FALSE)
    }
  }

  ### display qualifying changes
  filteredData <- filteredData %>%
    filter(!is.na(mvalue) & !is.na(change) & condition == curCondition)

  ## Plotly version
  filteredData <- filteredData %>%
    mutate(signChangeText = case_when(
      sigChangeSign == -1 ~ "Dephosphorylation",
      sigChangeSign == 0 ~ "Unchanged",
      sigChangeSign == 1 ~ "Phosphorylation"
    ))

  filteredData$text <- paste0("<b>Peptide ID:</b> ", filteredData$UniqueID, "<br>",
                              "<b>Median Centered Log<sub>2</sub> Intensity:</b> ", filteredData$mvalue, "<br>",
                              "<b>Log<sub>2</sub> Intensity Fold Change:</b> ", filteredData$change, "<br>",
                              "<b>Regulation Sign:</b> ", filteredData$signChangeText)
  # print(head(filteredData))
  if(plotly == TRUE){
    plot_ly(data = filteredData,
            x = ~mvalue,
            y = ~change,
            color = ~signChangeText,
            colors = c(regColors[1], regColors[2], regColors[3]),
            hovertemplate = ~text,
            type = "scatter",
            mode = "markers",
            marker = list(
              line = list(
                color = "#000000",
                width = 1
              )
            )
    ) %>%
      layout(
        xaxis = list(title = "<b>Median Centered Log<sub>2</sub> Intensity</b>"),
        yaxis = list(title = "<b>Log<sub>2</sub> Intensity Fold Change</b>"),
        legend = list(title = list(text = '<b> Regulation Sign</b>')),
        title = paste0("<b>Condition: ",curCondition, "</b>")
      )
  } else{
    filteredData %>%
      ggplot(.,aes(x=mvalue, y = change)) +
      geom_point(aes(fill = signChangeText), pch=21,  color = "black", size = 3) +
      labs(fill = "Regulation Sign",
           title = paste0("Condition: ", curCondition),
           x = expression(paste("Median Centered Log"[2]*" Intensity")),
           y = expression(paste("Log"[2]*" Intensity Fold Change"))) +
      scale_fill_manual(values = c(regColors[1], regColors[2], regColors[3]),
                        labels = c("Dephosphorylation", "Unchanged", "Phosphorylation")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  }
}

#' A VisNetwork plot of the inferred network
#'
#' @description Return the phosphorylation signaling network with the inferred causal relations between the kinases/phosphatases and downstream peptides.
#' @param netPhorceData Processed NetPhorce Object
#' @param condition (Required). Select a specific condition if multiple conditions are present.
#' @param FASTADescription (Optional). If TRUE, the hover information of the nodes will return the FASTA header column information for better identification of the phosphopeptides.
#' @return A visNetwork object
#' @export
#' @examples
#' \dontrun{
#' ## Loading Two Conditions Example
#' data("twoConditionsExample")
#' ## Identify the Key Columns
#' identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
#'                                     positionCol = "Position",
#'                                     reverseCol = "Reverse",
#'                                     localizationProbCol = "Localization prob",
#'                                     potentialContaminationCol = "Potential contaminant",
#'                                     aminoAcidCol = "Amino acid",
#'                                     uniqueIDCol = "Protein",
#'                                     seqWindowIDCol = "Sequence window",
#'                                     fastaIDCol = "Fasta headers")
#' ## Identify the Intensity Columns with Condition, Time Point and Replication Information
#' intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
#'                                          intensityPattern = "con_time_rep",
#'                                          verbose = TRUE)
#' ## Process the data based on the identified columns
#' netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
#'                              processedColNames = identifiedCols,
#'                              processedIntensity = intensityCols,
#'                              minReplication = 3,
#'                              minLocalProb = 0.75)
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
#' ## VisNetwork Plot with FASTA header information included
#' plotNetPhorce(netPhorceData = netPhorceData,
#'               condition = "tot3",
#'               FASTADescription = TRUE)
#' }
plotNetPhorce <- function(netPhorceData = netPhorceData,
                          condition = NULL,
                          FASTADescription = FALSE
                          ){
  ## Check for necessary network calculation data
  if(length(netPhorceData@networkPlotData) == 0){
    stop("No network data detected, please make sure you use requestPlotData = TRUE for networkAnalysis function.\n", call. = FALSE)
  }

  ## Check Condition
  ConDesign <- netPhorceData@Design$ConDesign
  nConditions = n_distinct(ConDesign$condition)
  if(is.null(condition)){
    if(nConditions == 1){
      curCondition = unique(ConDesign$condition)
    } else {
      cat("Multiple condition detected, please use the 'condition = ' arguemnt to select the specific condition.\n")
      cat("The following plot is based on the first found condition which is: \n")
      curCondition = unique(ConDesign$condition)[1]
      cat(curCondition)
      cat("\n")
    }
  } else {
    curCondition = condition
    if(!(curCondition %in% ConDesign$condition)){
      stop("Selected condition is not in the current data, the avaliable conditions are: ",
           paste(unique(ConDesign$condition), collapse = "; "), call. = FALSE)
    }
  }

  ## Check Presence of FASTA Header
  if(FASTADescription == TRUE){
    if(is.null(netPhorceData@Misc$FASTA_HeaderKey)){
      cat("FASTA Header column was not included during data processing step, please re-do the previous function to
          correctly identify the FASTA Header Column.\n")
      cat("No FASTA information is included in the following plot.")
    }else {
      FASTA_HeaderKey = netPhorceData@Misc$FASTA_HeaderKey
    }
  } else{
    FASTA_HeaderKey = NULL
  }

  data.scoring.genist <- netPhorceData@networkPlotData$data.scoring.genist
  data.scoring.sign <- netPhorceData@networkPlotData$data.scoring.sign
  Kinases <- netPhorceData@Misc$kinaseCheck$KinaseTable

  ## Part 1
  data.scoring.sign <- merge(data.scoring.sign, data.scoring.genist,
                             by = c("regulator", "target", "condition", "coReg"), all.x = TRUE)
  ## Filtering
  data.scoring.sign <- subset(data.scoring.sign, condition == curCondition)
  # print(data.scoring.sign)
  # print(head(data.scoring.sign))

  data.scoring.sign <-
    data.scoring.sign %>%
    mutate(Function = case_when(
      regulation == "phosphorylation" ~ 1,
      regulation == "dephosphorylation" ~ -1,
      regulation == "undetermined" ~ 0
    ))
  graphFile <- data.scoring.sign[, c("regulator", "target", "BDei", "Function", "coReg")]
  graphFile$StartGene <- graphFile$regulator
  graphFile$StopGene <- graphFile$target

  colnames(graphFile) <- c("regulator", "target", "Value", "Function", "coReg", "StartGene",  "StopGene")
  ## networkDS, threeejs, visNetowrk, ndtv
  graphFile <- graphFile[, c("StartGene", "StopGene", "Function", "Value", "coReg")]
  g <- graph.data.frame(graphFile[, c("StartGene", "StopGene")], directed = TRUE)

  # uniqueGenes <- unique(as.vector(unlist(graphFile[, c("StartGene", "StopGene")])))

  ## Part 3
  vertices <- data.frame(
    name = V(g)$name,
    group = walktrap.community(g)$membership,
    betweenness = as.vector(betweenness(g, directed=TRUE, normalized = TRUE))
  )

  # vertices <- merge(vertices, extractGenes, by.x = "name", by.y = "ensembl_transcript_id", all.x = TRUE)
  # print(head(vertices))

  # create indices (indexing needs to be JS format)
  graphFile$source.index = match(graphFile$StartGene, vertices$name) - 1
  graphFile$target.index = match(graphFile$StopGene, vertices$name) - 1
  # print(graphFile)

  ## Part 4
  nodes <- as.data.frame(vertices)
  colnames(nodes)[1] <- "id"
  nodes$label <- nodes$id
  nodes$ID <- sapply(strsplit(nodes$id, "\\."), "[", 1)

  edges <- as.data.frame(graphFile)
  edges$arrows.from.type <-  rep(NA, nrow(edges))
  edges$Function <- as.character(edges$Function)
  edges <-
    edges %>%
    mutate(arrows.to.type = case_when(
      Function == "1" ~ "arrow",
      Function == "-1" ~ "arrow",
      Function == "0" ~ "arrow"
    ),
    color = case_when(
      Function == "1" ~ "#2A77D5",
      Function == "-1" ~ "#D52A47",
      Function == "0" ~ "#656565"
    ),
    Functions = case_when(
      Function == "1" ~ "Phosphorylation",
      Function == "-1" ~ "Dephosphorylation",
      Function == "0" ~ "Undertemined"
    ),
    dashes = case_when(
      coReg == TRUE ~ TRUE,
      coReg == FALSE ~ FALSE
    )
    )

  colnames(edges)[1] <- "from"
  colnames(edges)[2] <- "to"
  edges$arrows.to.type <- ifelse(edges$arrows.to.type == "NA", NA, edges$arrows.to.type)
  # nodes$idNew <- sapply(strsplit(nodes$id, "_|\\."), "[", 1) #@Kuncheng: Is this correct ???
  #I am not sure which ID you want to achieve here, but for the Glycine you get "Glyma"
  #for all the peptides. And I don't think you will be able to extract the fasta header like this
  #Previously, I separate the name as follow to get the "unique ID" and remove
  #the splice variant : separate(UniqueID, into = c("Model_name", "AA", "multiplicity"), sep = "_", remove=FALSE) %>%
  #mutate(Model_name = sub("(\\.\\d{1,2})$", "", Model_name))
  nodes$UniqueID_V2 <- nodes$id
  nodes <- nodes %>%
    separate(UniqueID_V2, into = c("Part1", "multiplicity"), sep = "_(?=[0-9]+$)") %>%
    separate(Part1, into = c("Protein", "AA"), sep = "_(?=[0-9A-z;]+$)") %>%
    mutate(idNew = sub("(\\.\\d{1,2})$", "", Protein))
  # print("TESTSDGFSGSETT")
  # print(head(nodes))
  nodes$value <- nodes$betweenness

  # FASTA_HeaderKey = NULL
  if (!is.null(FASTA_HeaderKey)){ ## if the header column is present
    FASTA_HeaderKey <- FASTA_HeaderKey %>%
      mutate(idNew = sub("(\\.\\d{1,2})$", "", Protein))
    nodes <- merge(nodes, FASTA_HeaderKey, by.x = "idNew", by.y = "idNew", all.x = TRUE)
    nodes$Desc <- gsub(pattern = "(.{50})", "\\1<br>", nodes$`FASTA Header`)
    nodes$title <- paste0("<b>Peptide Name:</b> ", nodes$id, "<br>",
                          "<b>Description:</b> ", nodes$Desc, "<br>")
  } else {
    nodes$Desc <- NA
    nodes$Desc <- ifelse(is.na(nodes$Desc), nodes$label, nodes$Desc)
    nodes$title <- nodes$Desc
  }


  # keyData <- readRDS("./data/Kinase_Phosphatase.rds")
  # View(Kinases)
  nodes$idNew <- toupper(nodes$idNew)
  nodes <- merge(x = nodes, y = Kinases,
                 by.x = "idNew", by.y = "ID", all.x = TRUE)
  # View(kinaseTable)

  #@Kuncheng so here this IdNew for sure does not match the kinase ID.
  #I propose you add revisit the IdNew as I proposed above and change it temporaly to uppercase
  #e.g. nodes = nodes %>% mutate(ID = str_to_upper("IdNew)) %>% left_join(Kinases, by = c("ID))
  # print(head(nodes))
  ## Add shape by the
  nodes$shape = ifelse(nodes$FAMILY == "Kinase", "triangle",
                       ifelse(nodes$FAMILY == "Phosphatase", "diamond", NA))
  nodes$shape[is.na(nodes$shape)] = "dot"
  nodes <- nodes %>%
    group_by(id) %>%
    filter(row_number() == 1)
  edges$width <- rescale(edges$Value, to = c(0, 2))
  edges$title <- paste0("<p><b>", str_to_title(edges$Functions), "</b><br>",
                        "<b>Width: </b>", edges$width, "</p>")

  ## Legend
  ledges <- data.frame(
    label = c("Dephosphorylation", "Phosphorylation", "Undertemined", "Time lapse 1", "Time lapse 0"),
    color = c( "#D52A47", "#2A77D5", "#656565", "#000000", "#000000"),
    dashes = c(FALSE, FALSE, FALSE, TRUE, FALSE),
    arrows.from.type = c(NA, NA, NA, NA, NA),
    arrows.to.type = c("arrow", "arrow", "arrow", NA, NA),
    shadow = c(TRUE,TRUE,TRUE,TRUE,TRUE),
    font.align = "top"
    # arrows = c("to", "from",)
  )
  lnodes <- data.frame(
    label = c("Kinase", "Phosphatase", "Others"),
    shape = c("triangle", "diamond", "dot"),
    shadow = c(TRUE, TRUE, TRUE)
  )

  # View(nodes)
  nodes$value <- rescale(nodes$betweenness, to = c(0, 5))

  ## Part 5
  plot <- visNetwork(nodes, edges, height = "750px", width = "100%",
                     main = paste0("Condition: ", curCondition)) %>%
    visNodes(
      shadow = list(enabled = FALSE)
    ) %>%
    visOptions(highlightNearest = TRUE) %>%
    visLegend(width=0.1, addEdges = ledges,
              useGroups = FALSE, addNodes = lnodes,
              zoom = FALSE, main = "Legend") %>%
    visOptions(highlightNearest = TRUE
               # selectedBy = "group"
    ) %>%
    visInteraction(navigationButtons = TRUE) %>%
    visEdges(smooth = FALSE) %>%
    visPhysics(solver = "forceAtlas2Based", repulsion = "damping", minVelocity = 50)

  plot
}

































