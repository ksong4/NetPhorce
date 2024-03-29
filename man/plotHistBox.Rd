% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotHistBox}
\alias{plotHistBox}
\title{Plot the Distribution of the Intensity Data}
\usage{
plotHistBox(
  netPhorceData = netPhorceData,
  condition = NULL,
  histogram = TRUE,
  boxplot = TRUE
)
}
\arguments{
\item{netPhorceData}{(Required). Processed NetPhorce Object}

\item{condition}{(Required). Select a specific condition from your experiment.}

\item{histogram}{(Optional). Enable the histogram Plot. If the boxplot parameter is also TRUE, a side-by-side plot will be generated.}

\item{boxplot}{(Optional). Enable the boxplot. If the histogram parameter is also TRUE, a side-by-side plot will be generated.}
}
\value{
A ggplot object
}
\description{
The histogram and/or boxplot of the intensity data for each time point after normalization. The function generates a \code{\link{ggplot}} plot.
}
\examples{
\dontrun{
## Loading Two Conditions Example
data("twoConditionsExample")
## Identify the Key Columns
identifiedCols <- confirmColumnNames(
  rawMaxQuant = twoConditionsExample,
  positionCol = "Position",
  reverseCol = "Reverse",
  localizationProbCol = "Localization prob",
  potentialContaminationCol = "Potential contaminant",
  aminoAcidCol = "Amino acid",
  uniqueIDCol = "Protein",
  seqWindowIDCol = "Sequence window",
  fastaIDCol = "Fasta headers")
## Identify the Intensity Columns with Condition, Time Point and Replication Information
intensityCols <- confirmIntensityColumns(rawMaxQuant = twoConditionsExample,
                                         intensityPattern = "con_time_rep",
                                         verbose = TRUE)
## Process the data based on the identified columns
netPhorceData <- processData(rawMaxQuant = twoConditionsExample,
                             processedColNames = identifiedCols,
                             processedIntensity = intensityCols,
                             minReplication = 3,
                             minLocalProb = 0.75)
## Side-by-side Plot
plotHistBox(netPhorceData = netPhorceData,
             condition = "tot3",
             histogram = TRUE,
             boxplot = TRUE)
## Boxplot Only
plotHistBox(netPhorceData = netPhorceData,
             condition = "Col0",
             histogram = FALSE,
             boxplot = TRUE)
## Histogram Only
plotHistBox(netPhorceData = netPhorceData,
             condition = "tot3",
             histogram = TRUE,
             boxplot = FALSE)
}
}
