% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotDistribution}
\alias{plotDistribution}
\title{Plot the Distribution of the Intensity Data}
\usage{
plotDistribution(
  netPhorceData = netPhorceData,
  condition = NULL,
  plotly = FALSE
)
}
\arguments{
\item{netPhorceData}{(Required). Processed netPhorce Object}

\item{condition}{(Required). Select a specific condition from your experiment.}

\item{plotly}{(Required). If TRUE, output an interactive \code{\link{plotly}} version, else output a static \code{\link{ggplot2}} version.}
}
\value{
a ggplot/plotly object
}
\description{
Boxplots of the intensity data for each time point, condition, and replicate before and after variance stabilizing normalization (vsn). The function generates a ggplot object by default. An interactive plotly version (used on our Shiny application) is included if \emph{plotly = TRUE} is set.
}
\examples{
\dontrun{
## Loading Two Conditions Example
data("twoConditionsExample")
## Identify the Key Columns
identifiedCols <- confirmColumnNames(rawMaxQuant = twoConditionsExample,
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
## Generate Distribution Plot – GGPLOT version (Static)
plotDistribution(netPhorceData = netPhorceData,
                 condition = NULL,
                 plotly = FALSE)
## Generate Distribution Plot – PLOTLY version (Interactive)
plotDistribution(netPhorceData = netPhorceData,
                 condition = NULL,
                 plotly = TRUE)
}
}
