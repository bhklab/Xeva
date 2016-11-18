library(ggplot2)
#' Given a PDX object this will plot drug response curve
#' i.e. Time vs Volume
#' @param PDX a PDX object
#' @return ggplot object
#' @examples
#' plotDrugResponse(PDX)
#' @export
plotDrugResponse <- function(expSlot, expList)
{
  ## read the PDX object
  toPlt = readRDS("toPlot.Rda")
  expSlot=toPlt$expSlot
  expDesX=toPlt$expDesX




}