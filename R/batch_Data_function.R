if(1==2) ##----- do we need this ?
{
#' get data for a batch
#'
#' Get data for a batch id or experiment design
#'
#' @param object Xeva object
#' @param batchName batch name
#' @param expDig Experiment design list
#' @param treatment.only \code{FALSE}. If TRUE only in treatment data will be considered
#' @param vol.normal \code{FALSE} . If TRUE volume will ne normalised
#' @param impute.value \code{TRUE}, will impute values where missing
#'
#' @return A list with dataframs of control, treatment and control.mean, treatment.mean
#'
#' @examples
#' data(pdxe)
#' btdata <- batchData(pdxe, batchName="X-1228.CKX620", vol.normal=TRUE)
#'
#' @export
batchData <- function(object, batchName=NULL, expDig =NULL, treatment.only=FALSE,
                      vol.normal=FALSE, impute.value=TRUE)
{
  if(is.null(batchName) & is.null(expDig))
  { stop("please provide 'batchName' or 'expDig'") }

  if(!is.null(batchName))
  { expDig <- expDesign(object, batchName) }

  dfp <- list()
  if(length(expDig$control)>0)
  {
    dfp$control <- .getExperimentMultipalIDs(object, mids=expDig$control,
                                             treatment.only=treatment.only,
                                             vol.normal=vol.normal)

    dfp$control.mean <- .collapseRplicate(dfp$control, impute.value=impute.value)
  }

  if(length(expDig$treatment)>0)
  {
    dfp$treatment <- .getExperimentMultipalIDs(object, mids=expDig$treatment,
                                               treatment.only=treatment.only,
                                               vol.normal=vol.normal)

    dfp$treatment.mean <- .collapseRplicate(dfp$treatment, impute.value=impute.value)
  }

  return(dfp)
}



}
