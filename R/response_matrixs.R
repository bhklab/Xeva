
model_response_class <- function(name, value=NA, fit=NA)
{
  mr <- structure(list(name=name, value=value, fit=fit),
                  class = "modelResponse")
  return(mr)
}

#' Print the model response
#' @param x modelResponse object
#' @param ... Other arguments
#' @export
print.modelResponse <- function(x, ...)
{
  z <- sprintf("%s = %f\n", x$name, x$value)
  cat(z)
}


batch_response_class <- function(name, value=NA, control=NULL, treatment=NULL, ...)
{
  br <- structure(list(name=name, value=value, control=control, treatment=treatment),
                 class = "batchResponse")
  return(br)
}

#' Print the batch response
#' @param x batchResponse object
#' @param ... Other arguments
#' @export
print.batchResponse <- function(x, ...)
{
  cat(sprintf("%s = %f\n", x$name, x$value))

  if(!is.null(x$control))
  { cat(sprintf("control = %f\n", x$control$value)) }

  if(!is.null(x$treatment))
  { cat(sprintf("treatment = %f\n", x$treatment$value)) }
}



#' \code{setResponse} sets response of an Xeva object.
#'
#' @param object Xeva object.
#' @param res.measure Response measure; multiple measures are allowed.
#' @param min.time Default \strong{10} days. Used for \emph{mRECIST} computation.
#' @param treatment.only Default \code{FALSE}. If \code{TRUE}, give data for non-zero dose periods only (if dose data are available).
#' @param max.time Maximum time for data.
#' @param vol.normal Default \code{TRUE} will use
#' @param impute.value Default \code{FALSE}. If \code{TRUE}, impute the missing values.
#' @param concurrent.time Default \code{FALSE}. If \code{TRUE}, cut the batch data such that control and treatment will end at same time point.
#' @param verbose Default \code{TRUE} will print information
#'
#' @return  Returns updated Xeva object.
#'
#' @examples
#' data(brca)
#' \dontrun{ brca  <- setResponse(brca, res.measure = c("mRECIST")) }
#' @export
setResponse <- function(object, res.measure=c("mRECIST", "slope", "AUC", "angle", "abc", "TGI"),
                        min.time=10, treatment.only=FALSE, max.time=NULL,
                        vol.normal=TRUE, impute.value=TRUE, concurrent.time =TRUE,
                        verbose=TRUE)
{
  sen <- slot(object, "sensitivity")

  ###--------compute mRECIST ---------------------------------------------------
  if(any(c("mRECIST", "best.response", "best.average.response") %in% res.measure))
  {
    vl2compute <- c("mRECIST", "best.response", "best.response.time",
                    "best.average.response", "best.average.response.time")
    sen$model[, vl2compute[!(vl2compute%in%colnames(sen$model))]] <- NA

    for(mid in modelInfo(object)$model.id)
    {
      mr <- response(object, model.id=mid, res.measure="mRECIST",
                     treatment.only=treatment.only, max.time=max.time,
                     impute.value=impute.value, min.time=min.time,
                     concurrent.time=FALSE, verbose=verbose)
      for(si in vl2compute)
      { sen$model[mid, si] <- mr[[si]] }
    }
  }

  ###--------compute slope -----------------------------------------------------
  if("slope" %in% res.measure)
  {
    sen$model[, "slope"] <- NA
    for(mid in modelInfo(object)$model.id)
    {
      sl <- response(object, model.id=mid, res.measure="slope",
                     treatment.only=treatment.only, max.time=max.time,
                     impute.value=impute.value, min.time=min.time,
                     concurrent.time=FALSE, verbose=verbose)
      sen$model[mid, "slope"] <- sl$value #$slope
    }
  }

  ###--------compute AUC -------------------------------------------------------

  if("AUC" %in% res.measure)
  {
    sen$model[, "AUC"] <- NA
    for(mid in modelInfo(object)$model.id)
    {
      auc <- response(object, model.id=mid, res.measure="AUC",
                     treatment.only=treatment.only, max.time=max.time,
                     impute.value=impute.value, min.time=min.time,
                     concurrent.time=FALSE, verbose=verbose)
      sen$model[mid, "AUC"] <- auc$value
    }
  }

  ###--------compute doubling time ---------------------------------------------


  ##----------------------------------------------------------------------------
  ##-----------------for batch -------------------------------------------------

  ###--------compute angle for batch -------------------------------------------
  if("angle" %in% res.measure)
  {
    sen$batch[, c("slope.control", "slope.treatment", "angle")] <- NA
    #for(bid in batchNames(object))
    for(bid in batchInfo(object))
    {
      sl <- response(object, batch = bid, res.measure="angle",
                     treatment.only=treatment.only, max.time=max.time,
                     impute.value=impute.value, min.time=min.time,
                     concurrent.time=concurrent.time, verbose=verbose)
      sen$batch[bid, c("slope.control", "slope.treatment", "angle")] <-
        c(sl$control$value, sl$treatment$value, sl$value)
    }
  }

  ###--------compute abc for batch ---------------------------------------------
  if("abc" %in% res.measure)
  {
    sen$batch[, c("auc.control", "auc.treatment", "abc")] <- NA
    #for(bid in batchNames(object))
    for(bid in batchInfo(object))
    {
      sl <- response(object, batch = bid, res.measure="abc",
                     treatment.only=treatment.only, max.time=max.time,
                     impute.value=impute.value, min.time=min.time,
                     concurrent.time=concurrent.time, verbose=verbose)
      sen$batch[bid, c("auc.control", "auc.treatment", "abc")] <-
        c(sl$control$value, sl$treatment$value, sl$value)
    }
  }

###--------compute TGI for batch ---------------------------------------------

  if("TGI" %in% res.measure)
  {
    sen$batch[, c("TGI")] <- NA
    for(bid in batchInfo(object))
    {
      sl <- response(object, batch = bid, res.measure="TGI",
                     treatment.only=treatment.only, max.time=max.time,
                     impute.value=impute.value, min.time=min.time,
                     concurrent.time=concurrent.time, verbose=verbose)
      sen$batch[bid, c("TGI")] <- sl$value
    }
  }

  ##--------------code for batch level mR --------------------------------------
  ##----------------------------------------------------------------------------

  slot(object, "sensitivity") <- sen
  return(object)
}



#' compute response
#'
#' \code{response} Computes the drug response of a PDX model or batch.
#'
#' @param object Xeva object.
#' @param res.measure Drug response measure.
#' @param model.id \code{model.id} for which the durg response is to be computed.
#' @param batch \code{batch.id} or experiment design for which the drug response is to be computed.
#' @param treatment.only Default \code{FALSE}. If \code{TRUE}, give data for non-zero dose periods only (if dose data are available).
#' @param min.time Default \strong{10} days. Used for \emph{mRECIST} computation.
#' @param max.time Maximum time for data.
#' @param vol.normal Default \code{FALSE} will use
#' @param impute.value Default \code{FALSE}. If \code{TRUE}, impute the missing values.
#' @param concurrent.time Default \code{FALSE}. If \code{TRUE}, cut the batch data such that control and treatment will end at same time point.
#' @param verbose Default \code{TRUE} will print information.
#'
#' @return  Returns model or batch drug response object.
#'
#' @examples
#' data(brca)
#' response(brca, model.id="X.1004.BG98", res.measure="mRECIST")
#'
#' response(brca, batch="X-6047.paclitaxel", res.measure="angle")
#'
#' ed <- list(batch.name="myBatch", treatment=c("X.6047.LJ16","X.6047.LJ16.trab"),
#'              control=c("X.6047.uned"))
#' response(brca, batch=ed, res.measure="angle")
#'
#' @export
response <- function(object, model.id=NULL,
                     batch=NULL,
                     #res.measure=c("angle", "mRECIST", "AUC", "angle", "abc"),
                     res.measure=c("mRECIST", "slope", "AUC", "angle", "abc", "TGI"),
                     treatment.only=FALSE, max.time=NULL, impute.value=TRUE,
                     min.time=10, concurrent.time =TRUE, vol.normal=FALSE,
                     verbose=TRUE)
{
  if(is.null(model.id) & is.null(batch)) #Name) & is.null(expDig))
  { stop("'model.id', 'batch' all NULL") }

  model.measure <- c("mRECIST", "best.response", "best.response.time",
                     "best.average.response", "best.average.response.time",
                     "slope", "AUC")

  batch.measure <- c("angle", "abc", "TGI")

  ##------------- for model ----------------------------------------------------
  if(!is.null(model.id))
  {
    dl <- getExperiment(object, model.id=model.id[1], treatment.only=treatment.only,
                        max.time=max.time, vol.normal=vol.normal,
                        impute.value=impute.value)

    ###--------compute mRECIST -------------------------------------------------
    if(any(c("mRECIST", "best.response", "best.average.response") %in% res.measure))
    {
      if(verbose==TRUE) {cat(sprintf("computing mRECIST for %s\n", model.id))}
      mr <- mRECIST(dl$time, dl$volume, min.time=min.time, return.detail=TRUE)
      return(mr)
    }

    ###--------compute slope -----------------------------------------------------
    if(res.measure=="slope")
    {
      if(verbose==TRUE) {cat(sprintf("computing slope for %s\n", model.id))}
      return(slope(dl$time, dl$volume, degree=TRUE))
    }

    ###--------compute AUC -------------------------------------------------------
    if(res.measure=="AUC")
    {
      if(verbose==TRUE) {cat(sprintf("computing AUC for %s\n", model.id))}
      return( AUC(dl$time, dl$volume))
    }

  }

  ##-----------------for batch -------------------------------------------------
  if(is.null(model.id))
  {
    dl <- getExperiment(object, batch=batch,
                        treatment.only=treatment.only, max.time=max.time,
                        vol.normal=vol.normal, impute.value=impute.value,
                        concurrent.time=concurrent.time)

    cInd <- dl$batch$exp.type == "control"
    tInd <- dl$batch$exp.type == "treatment"

    contr.time <- contr.volume <- treat.time <- treat.volume <- NULL
    if(sum(cInd)>1)
    { contr.time <- dl$batch$time[cInd]; contr.volume <- dl$batch$mean[cInd] }

    if(sum(tInd)>1)
    { treat.time <- dl$batch$time[tInd]; treat.volume <- dl$batch$mean[tInd]}

    if(verbose==TRUE){
      if(is.character(batch))
        {bName <- batch} else {bName <- batch$batch.name}
      cat(sprintf("computing %s for batch %s\n",res.measure, bName))

      }
    ###--------compute angle for batch -----------------------------------------
    if(res.measure =="angle")
    {
      rtx <- angle(contr.time, contr.volume, treat.time,treat.volume, degree=TRUE)
      return(rtx)
    }

    ###--------compute abc for batch ---------------------------------------------
    if(res.measure=="abc")
    {
      rtx <- ABC(contr.time, contr.volume, treat.time, treat.volume)
      return(rtx)
    }

    ###--------compute abc for batch ---------------------------------------------
    if(res.measure=="TGI")
    {
      rtx <- TGI(contr.volume, treat.volume)
      return(rtx)
    }
  }
}
