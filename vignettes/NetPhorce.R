## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(NetPhorce)

## -----------------------------------------------------------------------------
## Loading One Condition Data
data("oneConditionExample", package = "NetPhorce")
DT::datatable(oneConditionExample, 
              rownames = FALSE, 
              options = list(
                pageLength = 5, 
                scrollX = TRUE,
                rownames = FALSE,
                autoWidth = TRUE, 
                columnDefs = list(list(targets=c(7), width = "700px"),
                                  list(className = 'dt.center', targets = "_all"))
              )) 

