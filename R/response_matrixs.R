
model_response_class <- function(name, value=NA, fit=NA)
{
  mr <- structure(list(name=name, value=value, fit=fit),
                  class = "modelResponse")
  return(mr)
}

print.modelResponse <- function(mr)
{
  z <- sprintf("%s = %f\n", mr$name, mr$value)
  cat(z)
}


batch_response_class <- function(name, value=NA, control=NA, treatment=NA)
{
  br <- structure(list(name=name, value=value, control=control, treatment=treatment),
                 class = "batchResponse")
  return(br)
}

print.batchResponse <- function(br)
{
  z <- sprintf("%s = %f\ncontrol = %f\ntreatment = %f\n", br$name, br$value,
               br$control$value, br$treatment$value)
  cat(z)
}



#' \code{setResponse} sets response of an Xeva object
#'
#' @param object Xeva object
#' @param res.measure response measure, multipal measure allowed
#' @param min.time default \strong{10} days. Used for \emph{mRECIST} computation
#' @param treatment.only Default \code{FALSE}. If TRUE give data only for non-zero dose periode (if dose data avalible)
#' @param max.time maximum time for data
#' @param vol.normal default \code{TRUE} will use
#' @param impute.value default \code{FALSE}. If TRUE will impute the values
#' @param concurrent.time default \code{FALSE}. If TRUE will cut the batch data such that control and treatment will end at same time point
#' @param verbose default \code{TRUE} will print infromation
#'
#' @return  returns updated Xeva object
#' @examples
#'
#' @export
setResponse <- function(object, res.measure=c("angle", "mRECIST", "AUC", "angle", "abc"),
                        min.time=10, treatment.only=TRUE, max.time=NULL,
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
      sen$model[mid, "slope"] <- sl$slope
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
      sen$model[mid, "AUC"] <- auc
    }
  }

  ###--------compute doubling time ---------------------------------------------


  ##----------------------------------------------------------------------------
  ##-----------------for batch -------------------------------------------------

  ###--------compute angle for batch -------------------------------------------
  if("angle" %in% res.measure)
  {
    sen$batch[, c("slope.control", "slope.treatment", "angle")] <- NA
    for(bid in batchNames(object))
    {
      sl <- response(object, batchName = bid, res.measure="angle",
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
    ##get abc from response

  }

  ##--------------code for batch level mR --------------------------------------

  slot(object, "sensitivity") <- sen
  return(object)
}

#'
#' @examples
#' data(pdxe)
#' response(pdxe, model.id="X.1270.LK36", res.measure="mRECIST")
#' response(pdxe, model.id="X.1270.LK36", res.measure="AUC")
#'
#' response(pdxe, batchName="X-1270.HDM201", res.measure="angle")
#' response(pdxe, batchName="X-1270.HDM201", res.measure="abc")
#' response(pdxe, batchName="X-007.BGJ398", res.measure="angle")
#'
#' object <- pdxe
#' batchName=c("X-007.BGJ398", "X-007.binimetinib", "X-007.BKM120", "X-007.BYL719")
#' model.id=expDig=max.time=NULL
#' treatment.only=vol.normal=impute.value=verbose=TRUE
#'
#'
#' batchName=expDig=max.time=NULL
#' model.id=c("X.1228.LC61.pael","X.1228.pael")
#' treatment.only=vol.normal=impute.value=verbose=TRUE
#'
#' batchName=model.id=max.time=NULL
#' treatment.only=vol.normal=impute.value=verbose=TRUE
#' b1 <- list(batch.name="b1", treatment=c("X.1228.LC61.pael","X.1228.pael"), control=c("X.1228.uned"))
#' b2 <- list(batch.name="b2", treatment=c("X.1270.LK36"), control=c("X.1270.uned"))
#' expDig <- list(b1, b2)

response <- function(object, model.id=NULL, batchName=NULL, expDig=NULL,
                     res.measure=c("angle", "mRECIST", "AUC", "angle", "abc"),
                     treatment.only=TRUE, max.time=NULL, impute.value=TRUE,
                     min.time=10, concurrent.time =TRUE, verbose=TRUE)
{
  if(is.null(model.id) & is.null(batchName) & is.null(expDig))
  { stop("'model.id', 'batchName' and 'expDig' all NULL") }

  model.measure <- c("mRECIST", "best.response", "best.response.time",
                     "best.average.response", "best.average.response.time",
                     "slope", "AUC")

  batch.measure <- c("angle")

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
    dl <- getExperiment(object, batchName=batchName, expDig=expDig,
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

    printMessage <- function(meas, batchName, expDig)
    {
      if(!is.null(batchName))
      { cat(sprintf("computing %s for batch %s\n",meas, bid)) } else
      {
        if(!is.null(expDig$batch.name))
        {cat(sprintf("computing %s for batch %s\n",meas,expDig$batch.name)) } else
        {cat(sprintf("computing %s for experiemnt\n", meas))}
      }
    }

    if(verbose==TRUE){ printMessage(res.measure, batchName, expDig) }
    ###--------compute angle for batch -----------------------------------------
    if(res.measure =="angle")
    {
      # vl.con <- vl.tre <- vl.batch <- NA
      # if(sum(cInd)>1)
      # { vl.con <- slope(dl$batch$time[cInd], dl$batch$mean[cInd], degree=TRUE) }
      #
      # if(sum(tInd)>1)
      # { vl.tre <- slope(dl$batch$time[tInd], dl$batch$mean[tInd], degree=TRUE) }
      #
      # if(all(c(!is.na(vl.con), !is.na(vl.tre))))
      # { vl.batch <- vl.con$slope - vl.tre$slope}
      #
      # rtx <- batch_response_class(name = "angle", value = vl.batch,
      #                              control=vl.con,  treatment=vl.tre)

      rtx <- angle(contr.time, contr.volume, treat.time,treat.volume, degree=TRUE)
      return(rtx)
    }

    ###--------compute abc for batch ---------------------------------------------
    if(res.measure=="abc")
    {
      # vl.con <- vl.tre <- vl.batch <- NA
      # if(sum(cInd)>1)
      # { vl.con <- AUC(dl$batch$time[cInd], dl$batch$mean[cInd]) }
      #
      # if(sum(tInd)>1)
      # { vl.tre <- AUC(dl$batch$time[tInd], dl$batch$mean[tInd]) }
      #
      # if(all(c(!is.na(vl.con), !is.na(vl.tre))))
      # { vl.batch <- vl.con - vl.tre}
      #
      # rtx <- batch_response_class(name = "abc", value = vl.batch,
      #                              control=vl.con,  treatment=vl.tre)
      rtx <- ABC(contr.time, contr.volume, treat.time, treat.volume)
      return(rtx)
    }
  }

}
