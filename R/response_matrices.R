
model_response_class <- function(name, value=NA, fit=NA)
{
  mr <- structure(list(name=name, value=value, fit=fit),
                  class = "modelResponse")
  return(mr)
}

#' Print the model response
#' @param x modelResponse object
#' @param ... Other arguments
#' @return prints the modelResponse
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
#' @return prints the batchResponse
#' @export
print.batchResponse <- function(x, ...)
{
  cat(sprintf("%s = %f\n", x$name, x$value))

  if(!is.null(x$control))
  { cat(sprintf("control = %f\n", x$control$value)) }

  if(!is.null(x$treatment))
  { cat(sprintf("treatment = %f\n", x$treatment$value)) }
}



#' set PDX response
#'
#' \code{setResponse} sets response of all PDXs in an Xeva object.
#'
#' @param object Xeva object.
#' @param res.measure Response measure, multiple measures are allowed. See `Details` below
#' @param min.time Minimum number of days for \emph{mRECIST} computation. Default \strong{10} days.
#' @param treatment.only Default \code{FALSE}. If \code{TRUE}, give data for non-zero dose periods only (if dose data are available).
#' @param max.time Maximum number of days to consider for analysis. Data byond this will be discarded. Default \code{NULL} takes full data.
#' @param vol.normal If TRUE it will will normalize the volume. Default \code{FALSE}
#' @param impute.value Default \code{FALSE}. If \code{TRUE}, impute the missing volume values.
#' @param concurrent.time Default \code{FALSE}. If \code{TRUE}, cut the batch data such that control and treatment will end at same time point.
#' @param log.volume If TRUE log of the volume will be used for response calculation. Default \code{FALSE}
#' @param verbose Default \code{TRUE} will print information.
#'
#' @return  Returns updated Xeva object.
#'
#' @details At present fellowing response measure are implemented
#' * mRECIST Computes mRECIST for indivial PDX model
#' * slope Computes slope of the fitted indivial PDX curve
#' * AUC  Computes area under a PDX curve for indivial PDX model
#' * angle Computes angle between treatment and control PDX curves
#' * abc Computes area between the treatment and control PDX curves
#' * TGI Computes  tumor growth inhibition using treatment and control PDX curves
#' * lmm Computes linear mixed model (lmm) statistics for a PDX batch
#' @md
#'
#' @examples
#' data(brca)
#' brca  <- setResponse(brca, res.measure = c("mRECIST"), verbose=FALSE)
#' @export
setResponse <- function(object,
                        res.measure=c("mRECIST", "slope", "AUC", "angle", "abc", "TGI", "lmm"),
                        min.time=10, treatment.only=FALSE, max.time=NULL,
                        vol.normal=FALSE, impute.value=TRUE, concurrent.time =TRUE,
                        log.volume=FALSE, verbose=TRUE)
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
                     concurrent.time=FALSE,
                     vol.normal=vol.normal,
                     log.volume=log.volume, verbose=verbose)
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
                     concurrent.time=FALSE,
                     vol.normal=vol.normal, log.volume=log.volume,
                     verbose=verbose)
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
                     concurrent.time=FALSE,
                     vol.normal=vol.normal, log.volume=log.volume,
                     verbose=verbose)
      sen$model[mid, "AUC"] <- auc$value
    }
  }

  ##----------------------------------------------------------------------------
  ##-----------------for batch -------------------------------------------------

  ###--------compute angle for batch -------------------------------------------
  if("angle" %in% res.measure)
  {
    sen$batch[, c("slope.control", "slope.treatment", "angle")] <- NA
    for(bid in batchInfo(object))
    {
      sl <- response(object, batch = bid, res.measure="angle",
                     treatment.only=treatment.only, max.time=max.time,
                     impute.value=impute.value, min.time=min.time,
                     concurrent.time=concurrent.time,
                     vol.normal=vol.normal, log.volume=log.volume,
                     verbose=verbose)
      sen$batch[bid, c("slope.control", "slope.treatment", "angle")] <-
        c(sl$control$value, sl$treatment$value, sl$value)
    }
  }

  ###--------compute abc for batch ---------------------------------------------
  if("abc" %in% res.measure)
  {
    sen$batch[, c("auc.control", "auc.treatment", "abc")] <- NA
    for(bid in batchInfo(object))
    {
      sl <- response(object, batch = bid, res.measure="abc",
                     treatment.only=treatment.only, max.time=max.time,
                     impute.value=impute.value, min.time=min.time,
                     concurrent.time=concurrent.time,
                     vol.normal=vol.normal, log.volume=log.volume,
                     verbose=verbose)
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
                     concurrent.time=concurrent.time,
                     vol.normal=vol.normal, log.volume=log.volume,
                     verbose=verbose)
      sen$batch[bid, c("TGI")] <- sl$value
    }
  }

  ###--------compute lmm for batch ---------------------------------------------
  if("lmm" %in% res.measure)
  {
    sen$batch[, c("lmm")] <- NA
    for(bid in batchInfo(object))
    {
      sl <- response(object, batch = bid, res.measure="lmm",
                     treatment.only=treatment.only, max.time=max.time,
                     impute.value=impute.value, min.time=min.time,
                     concurrent.time=concurrent.time,
                     vol.normal=vol.normal,
                     log.volume=log.volume, verbose=verbose)
      sen$batch[bid, c("lmm")] <- sl$value
    }
  }
  ##--------------code for batch level mR --------------------------------------
  ##----------------------------------------------------------------------------

  slot(object, "sensitivity") <- sen
  return(object)
}



#' compute PDX response
#'
#' \code{response} Computes the drug response of an individual PDX model or batch.
#'
#' @param object Xeva object.
#' @param model.id \code{model.id} for which the durg response is to be computed.
#' @param batch \code{batch.id} or experiment design for which the drug response is to be computed.
#' @param res.measure Drug response measure. See `Details` below
#' @param treatment.only Default \code{FALSE}. If \code{TRUE}, give data for non-zero dose periods only (if dose data are available).
#' @param max.time Maximum time for data.
#' @param impute.value Default \code{FALSE}. If \code{TRUE}, impute the missing values.
#' @param min.time Default \strong{10} days. Used for \emph{mRECIST} computation.
#' @param concurrent.time Default \code{FALSE}. If \code{TRUE}, cut the batch data such that control and treatment will end at same time point.
#' @param vol.normal If TRUE it will normalize the volume. Default \code{FALSE}.
#' @param log.volume If TRUE log of the volume will be used for response calculation. Default \code{FALSE}
#' @param verbose Default \code{TRUE} will print information.
#'
#' @return  Returns model or batch drug response object.
#'
#' @details At present the following response measures are implemented
#' * mRECIST Computes mRECIST for individual PDX models
#' * slope Computes slope of the fitted individual PDX curves
#' * AUC  Computes area under a PDX curve for individual PDX models
#' * angle Computes angle between treatment and control PDX curves
#' * abc Computes area between the treatment and control PDX curves
#' * TGI Computes tumor growth inhibition using treatment and control PDX curves
#' * lmm Computes linear mixed model (lmm) statistics for a PDX batch
#' @md
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
response <- function(object, model.id=NULL, batch=NULL,
                     res.measure=c("mRECIST", "slope", "AUC", "angle", "abc", "TGI", "lmm"),
                     treatment.only=FALSE, max.time=NULL, impute.value=TRUE,
                     min.time=10, concurrent.time =TRUE, vol.normal=FALSE,
                     log.volume=FALSE, verbose=TRUE)
{
  if(is.null(model.id) & is.null(batch)) #Name) & is.null(expDig))
  { stop("'model.id', 'batch' all NULL") }

  ##------------- for model ----------------------------------------------------
  if(!is.null(model.id))
  {
    dl <- getExperiment(object, model.id=model.id[1], treatment.only=treatment.only,
                        max.time=max.time, vol.normal=vol.normal,
                        log.volume=log.volume, impute.value=impute.value)

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

    ###--------compute lmm for batch ---------------------------------------------
    if(res.measure=="lmm")
    {
      rtx <- lmm(dl$model)
      return(rtx)
    }
  }
}
