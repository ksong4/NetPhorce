#' Mofidy DEP function to fit the format of Processed MaxQuant Data
#'
#' @description The original functions from DEP pacakge can be found at:
#' https://bioconductor.org/packages/release/bioc/html/DEP.html
#' The only modification made is replace "condition" with "experiement" to match
#' the NetPhorce definition of each of these two terms. normalize_vsn and theme_DEP1
#' also listed here as necessary functions for plot_normalization.
#'
#' @param se (Required). SummarizedExperiment object from the NetPhorce Object, log2-transformed raw count.
#' @param se_norm (Required). SummarizedExperiment object from the NetPhorce Object. normalized log2-transformed raw count.
#' @param filterCondition (Required). Experimental condition.
#' @importFrom magrittr %>%
#' @import dplyr tidyr purrr tibble methods SummarizedExperiment vsn stringr utils ggplot2 ggpubr ellipse visNetwork dynamicTreeCut forcats
#' @rawNamespace import(plotly, except = c(last_plot))
#' @rawNamespace import(assertthat, except = c(has_name))
#' @import qvalue qvalue
#' @importFrom vsn vsnMatrix predict
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom stats median
#' @importFrom stats quantile
#' @importFrom stats prcomp
#' @importFrom stats cov
#' @importFrom colorspace darken
#' @importFrom DT datatable
#' @importFrom igraph V
#' @importFrom igraph walktrap.community
#' @importFrom igraph betweenness
#' @importFrom igraph graph.data.frame
#' @importFrom scales rescale
#' @keywords internal
#' @return plot object
#' @examples
#' NULL
#### Misc Function (Outside the server class) ####
##### Replacing plot_normalization from DEP package #####
plot_normalization <- function(se, se_norm, filterCondition = NULL) {
  # Get arguments from call
  call <- match.call()
  # print(call)

  arglist <- lapply(call[-1], function(x) x)

  arglist[[3]] <- NULL
  # print(arglist)
  # View(arglist)

  var.names <- vapply(arglist, deparse, character(1))
  arglist <- lapply(arglist, eval.parent, n = 2)
  names(arglist) <- var.names
  # print(arglist)

  # Show error if inputs are not the required classes
  lapply(arglist, function(x) {
    assertthat::assert_that(inherits(x,
                                     "SummarizedExperiment"),
                            msg = "input objects need to be of class 'SummarizedExperiment'")
    if (any(!c("label", "ID", "experiment", "replicate") %in% colnames(colData(x)))) {
      stop("'label', 'ID', 'experiment' and/or 'replicate' ",
           "columns are not present in (one of) the input object(s)",
           "\nRun make_se() or make_se_parse() to obtain the required columns",
           call. = FALSE)
    }
  })

  # Function to get a long data.frame of the assay data
  # annotated with sample info
  gather_join <- function(se) {
    assay(se) %>%
      data.frame() %>%
      gather(ID, val) %>%
      left_join(., data.frame(colData(se)), by = "ID")
  }

  df <- map_df(arglist, gather_join, .id = "var") %>%
    mutate(var = factor(var, levels = names(arglist)))
  # print(head(df))
  if(!is.null(filterCondition)){
    df <- subset(df, condition == filterCondition)
  }
  # Boxplots for conditions with facet_wrap
  # for the original and normalized values
  # return(df)
  ggplot(df, aes(x = ID, y = val, fill = experiment)) +
    geom_boxplot(notch = TRUE, na.rm = TRUE) +
    coord_flip() +
    facet_wrap(~var, ncol = 1) +
    labs(x = "", y = expression(log[2]~"Intensity")) +
    theme_DEP1()
}
##### Replacing normalize_vsn from DEP package #####
normalize_vsn <- function(se) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

  # Variance stabilization transformation on assay data
  se_vsn <- se
  vsn.fit <- vsn::vsnMatrix(2 ^ assay(se_vsn))
  assay(se_vsn) <- vsn::predict(vsn.fit, 2 ^ assay(se_vsn))
  return(se_vsn)
}

##### Replcing them_DEP1() from DEP Package #####
theme_DEP1 <- function() {
  # Use theme_bw() as default
  basesize <- 12
  theme <- ggplot2::theme_bw(base_size = basesize)

  # Change plot title appearance
  theme$plot.title$face <- "bold"
  theme$plot.title$size <- basesize + 2
  theme$plot.title$hjust <- 0.5

  # Change axis title appearance
  theme$axis.title.x$size <- basesize + 2

  theme$axis.title.y$size <- basesize + 2


  theme$axis.text$size <- basesize
  theme$axis.text$colour <- "black"

  # Change legend title appearance
  theme$legend.title$size <- basesize + 2

  # Change legend text appearance
  theme$legend.text$size <- basesize

  # Change strip text (facet headers) appearance
  theme$strip.text$face <- "bold"
  theme$strip.text$size <- basesize + 2
  theme$strip.text$colour <- "black"

  return(theme)
}
