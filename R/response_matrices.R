#' compute PDX volume
#' \code{volume} computes PDXs tumor volume from tumor width and length
#'
#' @param w width or first measurement vector. See \code{Details}
#' @param l length or second measurement vector. See \code{Details}
#'
#' @return  Returns a vector of volume
#'
#' @details This function applies the (l x w^2)/2 formula to compute volume. However
#' before computing volume it sorts data so that the smaller value is always taken as
#' width and  the bigger value as length.
#'
#' @examples
#' volume(w=1:5, l=6:10)
#' ## smaller value is always considered as width.
#' volume(w=6:10, l=1:5)
#'
#' @export
volume <- function(w, l)
{
  v = c()
  for(i in 1:max(c(length(w), length(l))))
  {
    mr <- sort(c(w[i], l[i]))
    v <- c(v, ((mr[1]**2) * mr[2]/2 ))
  }
  return(v)
}



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
  if(x$name=="mRECIST")
  {
    z <- sprintf("mRECIST = %s\nBest average response = %3.4f\nBest response = %3.4f\n",
                 x$fit$mRECIST, x$fit$best.average.response, x$fit$best.response)
    cat(z)
  } else
  {
  z <- sprintf("%s = %f\n", x$name, x$value)
  cat(z)
  }
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
  if(x$name=="bmRECIST")
  {
    cat(sprintf("%s\n", x$name))
    if(!is.null(x$control))
    { cat(sprintf("control = %s\n", x$control$value)) }
    if(!is.null(x$treatment))
    { cat(sprintf("treatment = %s\n", x$treatment$value)) }

  } else
  {
    cat(sprintf("%s = %f\n", x$name, x$value))

    if(!is.null(x$control))
    { cat(sprintf("control = %f\n", x$control$value)) }

    if(!is.null(x$treatment))
    { cat(sprintf("treatment = %f\n", x$treatment$value)) }
  }
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
#' * bmRECIST Computes mRECIST for control and treatment arms of a PDX batch
#' @md
#'
#' @examples
#' data(brca)
#' brca  <- setResponse(brca, res.measure = c("mRECIST"), verbose=FALSE)
#' @export
setResponse <- function(object,
                        res.measure=c("mRECIST", "slope", "AUC",
                                      "angle", "abc", "TGI", "lmm", "bmRECIST"),
                        min.time=10, treatment.only=TRUE, max.time=NULL,
                        vol.normal=FALSE, impute.value=TRUE, concurrent.time =FALSE,
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
      { sen$model[mid, si] <- mr$fit[[si]] }
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
      #cat(sprintf("TGI = %f", sl$value))
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
                     log.volume=log.volume,
                     verbose=verbose)
      sen$batch[bid, c("lmm")] <- sl$value
    }
  }
  ##--------------code for batch level mR --------------------------------------
  ##----------------------------------------------------------------------------
  ###--------compute bmRECIST for batch ----------------------------------------
  if("bmRECIST" %in% res.measure)
  {
    mrvars <- c("mRECIST", "best.response", #"best.response.time",
                "best.average.response" #, "best.average.response.time"
                )
    var.control  <- paste0(mrvars, ".control")
    var.treatment<- paste0(mrvars, ".treatment")

    sen$batch[, c(var.control, var.treatment)] <- NA
    for(bid in batchInfo(object))
    {
      sl <- response(object, batch = bid, res.measure="bmRECIST",
                     treatment.only=treatment.only, max.time=max.time,
                     impute.value=impute.value, min.time=min.time,
                     concurrent.time=concurrent.time,
                     vol.normal=vol.normal, log.volume=log.volume,
                     verbose=verbose)

      for(va in mrvars)
      {
        vac <- paste0(va, ".control")
        sen$batch[bid, vac] <- sl$control$fit[[va]]

        vat <- paste0(va, ".treatment")
        sen$batch[bid, vat] <- sl$control$fit[[va]]
      }
    }
  }

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
#' @param treatment.only Default \code{TRUE}. If \code{TRUE}, give data for time>=0 and if \code{FALSE}, uses all avaliable data.
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
#' * bmRECIST Computes mRECIST for control and treatment arms of a PDX batch
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
                     res.measure=c("mRECIST", "slope", "AUC",
                                   "angle", "abc", "TGI", "lmm", "bmRECIST"),
                     treatment.only=FALSE, max.time=NULL, impute.value=TRUE,
                     min.time=10, concurrent.time =FALSE, vol.normal=FALSE,
                     log.volume=FALSE, verbose=TRUE)
{
  if(is.null(model.id) & is.null(batch)) 
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
      mr <- mRECIST(dl$time, dl$volume, min.time=min.time)#, return.detail=TRUE)
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
    if(verbose==TRUE){
      if(is.character(batch))
      {bName <- batch} else {bName <- batch$batch.name}
      cat(sprintf("computing %s for batch %s\n",res.measure, bName))
    }
    
    dl <- getExperiment(object, batch=batch,
                        treatment.only=treatment.only, max.time=max.time,
                        vol.normal=vol.normal, impute.value=impute.value,
                        concurrent.time=concurrent.time)

    cInd <- dl$batch$exp.type == "control"
    tInd <- dl$batch$exp.type == "treatment"

    contr.time <- contr.volume <- treat.time <- treat.volume <- NULL
    
    rtx <- batch_response_class(name=res.measure, value=NA)
    
    if(sum(cInd)>1)
    { 
      contr.time <- dl$batch$time[cInd]
      contr.volume <- dl$batch$mean[cInd] 
    }

    if(sum(tInd)>1)
    { 
      treat.time <- dl$batch$time[tInd]; treat.volume <- dl$batch$mean[tInd]
    }
    
    .checkInputLen <- function(contr.time, contr.volume, treat.time,treat.volume)
    {
      rt <- TRUE
      if(length(contr.time) < 2 | length(contr.volume) < 2 )
      {
        msg <- sprintf("In batch %s\ncontrol not present or only 1 data point for given time, setting response to NA", batch)
        warning(msg)
        rt <- FALSE
      }
      
      if(length(treat.time) < 2 | length(treat.volume) < 2 )
      {
        msg <- sprintf("In batch %s\ntreatment not present or only 1 data point for given time, setting response to NA", batch)
        warning(msg)
        rt <- FALSE
      }
      return(rt)
    }

    ###--------compute angle for batch -----------------------------------------
    if(res.measure =="angle")
    {
      if(.checkInputLen(contr.time, contr.volume, treat.time,treat.volume)==TRUE)
      {
        rtx <- angle(contr.time, contr.volume, treat.time,treat.volume, degree=TRUE)
      }
      return(rtx)
    }

    ###--------compute abc for batch ---------------------------------------------
    if(res.measure=="abc")
    {
      if(.checkInputLen(contr.time, contr.volume, treat.time,treat.volume)==TRUE)
      { rtx <- ABC(contr.time, contr.volume, treat.time, treat.volume) }
      return(rtx)
    }

    ###--------compute TGI for batch ---------------------------------------------
    if(res.measure=="TGI")
    {
      if(.checkInputLen(contr.time, contr.volume, treat.time,treat.volume)==TRUE)
      { rtx <- TGI(contr.volume, treat.volume)}
      return(rtx)
    }

    ###--------compute lmm for batch ---------------------------------------------
    if(res.measure=="lmm")
    {
      if(.checkInputLen(contr.time, contr.volume, treat.time,treat.volume)==TRUE)
      {
        rtx <- list(value=NA)
        exp.type <- unique(dl$model$exp.type)
        if(length(exp.type)==2 &
           "control" %in% exp.type &
           "treatment" %in% exp.type)
        {
          rtx <- lmm(dl$model, log.volume)
        } else
        {
          warning("'control' or 'treatment' missing from exp.type or batch")
        }
      }
      return(rtx)
    }

    ###--------compute bmRECIST for batch ---------------------------------------------
    if(res.measure=="bmRECIST")
    {
      rtx <- bmRECIST(contr.time, contr.volume, treat.time,treat.volume, min.time=10)
      return(rtx)
    }

  }
}
