library(ggplot2)
#' Given a PDX object this will plot drug response curve
#' i.e. Time vs Volume
#' @param PDX a PDX object
#' @return ggplot object
#' @examples
#' plotDrugResponse(PDX)
#' @export
plotDrugResponse <- function(PDX)
{
  ## read the PDX object
  PDX =readRDS("data/PDX_test_object.Rda")

}