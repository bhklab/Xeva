.normalizeVolume <- function(X)
{
  if(is.na(X[1])==TRUE)
  {
    warning("First value not numeric.")
    return(rep(NA, length(X)))
  }
  if(X[1]==0)
  {
    warning("start volume zero, adding 1 to compute volume.normal")
    X <- X+1
  }
  rtx <- (X- X[1])/X[1]
  return(rtx)
}


.smoothCurve <- function(trDF, tsDF, x, y, approx.method="approx")
{
  lmdf <- data.frame(X=trDF[,x], Y=trDF[,y])
  prdf <- data.frame(X=tsDF[,x])
  if(approx.method=="lm")
  {
    lmod <- lm(Y~X, data = lmdf)
    tsDF[,y] <- predict(lmod, prdf, se.fit = F)
  }

  if(approx.method=="approx")
  {
    prdValue <- approx(lmdf$X, lmdf$Y, xout = prdf$X)
    tsDF[,y] <- prdValue$y
  }
  return(tsDF)
}

.smoothModel <- function(dw, timeVec, var = "volume")
{
  dw$impute.value <- "NO"
  t2imp <- setdiff(timeVec, dw$time)
  t2imp <- t2imp[t2imp <= max(dw$time)]
  if(length(t2imp)>0)
  {
      nwRws <- data.frame(matrix(data = NA, nrow = length(t2imp), ncol = ncol(dw)),
                          stringsAsFactors = F)
      colnames(nwRws) <- colnames(dw)
      nwRws$time <- t2imp
      nwRws$model.id <- dw$model.id[1]
      nwRws$drug.join.name <- dw$drug.join.name[1]
      nwRws$impute.value <- "YES"
      nwRws <- .smoothCurve(dw, nwRws, x="time", y=var)
      dt <- rbind(dw, nwRws)
      dt <- BBmisc::sortByCol(dt, c("time"))
      return(dt)
  }
  return(dw)
}

.getExperimentDataFromAExpID <- function(object, model.id, treatment.only)
{
  mod <- slot(object, "experiment")[[model.id]]
  if(is.null(mod))
  {
    msg <- sprintf("model.id '%s' not present in object", model.id)
    stop(msg)
  }

  mod.data <- mod$data
  mod.data$model.id <- mod$model.id
  mod.data$drug.join.name <- mod$drug$join.name

  mod.data = .removeNAcol(mod.data)
  mod.data = .reorderCol(mod.data, "model.id", 1)
  mod.data = .reorderCol(mod.data, "drug.join.name", 2)

  if(treatment.only==TRUE & !is.null(mod.data$dose))
  {
    tretIndx = extractBetweenTags(mod.data$dose, start.tag=0, end.tag=0)
    mod.data = mod.data[tretIndx, ]
  }

  mod.data$volume.normal <- NA
  mod.data$volume.normal <- .normalizeVolume(mod.data$volume)
  return(mod.data)
}

.getExperimentMultipalIDs <- function(object, mids, treatment.only=TRUE,
                                      max.time=NULL, return.list=TRUE,
                                      impute.value=FALSE, var="volume")
{
  rtx <- list()
  for(i in mids)
  {
    miD <- .getExperimentDataFromAExpID(object, model.id=i, treatment.only=treatment.only)
    if(!is.null(max.time)) { miD <- miD[miD$time<=max.time, ]}
    rtx[[i]] <- miD
  }

  if(length(rtx)>0)
  {
    if(impute.value==TRUE)
    {
      inLst2 <- list()
      timeVec <- sort(unique(unlist(lapply(rtx, "[[", "time" ))))
      for(mid in names(rtx))
      {
        inLst2[[mid]] <- .smoothModel(rtx[[mid]], timeVec=timeVec, var =var)
        inLst2[[mid]]$volume.normal <- .normalizeVolume(inLst2[[mid]]$volume)
      }
      rtx <- inLst2
    }

    if(return.list==FALSE)
    { rtx <- .rbindListOfDataframs(rtx) }
  }
  return(rtx)
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# To tolrate both types of expDig list
# expDig <- list(batch.name="myBatch",treatment=c("a","b"), control=c("c","d"))
# .reformatExpDesig(expDig)
# .reformatExpDesig(list(myBatch=expDig))
#
.reformatExpDesig <- function(expDig)
{
  if(any(c("batch.name", "treatment", "control") %in% names(expDig)))
  {
    expDig <- list(expDig)
  }
  names(expDig) <- sapply(expDig, "[[", "batch.name")
  return(expDig)
}
.getBatchData <- function(object, batchName=NULL, expDig =NULL, treatment.only=FALSE,
                          max.time=NULL, return.list=TRUE, impute.value=FALSE)
{
  if(is.null(batchName) & is.null(expDig))
  { stop("please provide 'batchName' or 'expDig'") }

  if(!is.null(batchName))
  { expDig <- expDesign(object, batchName) }
  expDig <- .reformatExpDesig(expDig)[[1]]

  dfp <- list()
  if(length(expDig$control)>0)
  {
    dfp$control <- .getExperimentMultipalIDs(object, mids=expDig$control,
                                             treatment.only=treatment.only,
                                             max.time=max.time,
                                             return.list = return.list,
                                             impute.value=impute.value)
    if(return.list==TRUE)
    {
      for(ci in 1:length(dfp$control))
      { dfp$control[[ci]]$exp.type <- "control"}
    } else{ dfp$control$exp.type <- "control"}
  }

  if(length(expDig$treatment)>0)
  {
    dfp$treatment <- .getExperimentMultipalIDs(object, mids=expDig$treatment,
                                               treatment.only=treatment.only,
                                               max.time=max.time,
                                               return.list = return.list,
                                               impute.value=impute.value)

    if(return.list==TRUE)
    {
      for(ti in 1:length(dfp$treatment))
      { dfp$treatment[[ti]]$exp.type <- "treatment"}
    } else{ dfp$treatment$exp.type <- "treatment"}
  }

  if(return.list==TRUE)
  {return(dfp)}

  dfp <- .rbindListOfDataframs(dfp)
  return(dfp)
}

##--------------------------------------------------------------------------------------------------
##----- get experiment data in flat data.fram ------------------------------------------------------
#' For a given  model.id, it will return a data.fram
#' containing all data stored in experiment slot
#'
#' @examples
#' data(pdxe)
#' getExperiment(pdxe, model.id="X.6047.uned", treatment.only=TRUE)
#' getExperiment(pdxe, model.id=c("X.6047.uned", "X.6047.pael"), treatment.only=TRUE)
#' getExperiment(pdxe, batchName="X-6047.paclitaxel", treatment.only=TRUE)
#'
#' @param object The \code{XevaSet}
#' @param model.id The \code{model.id} for which data is required, multipal allowed
#' @param batchName batch name from the Xeva set
#' @param expDig Experiment design
#' @param treatment.only Default \code{FALSE}. If TRUE give data only for non-zero dose periode (if dose data avalible)
#' @param max.time maximum time for data
#' @param return.list default \code{FALSE} will return a datafram
#'
#' @return a \code{data.fram} will all the the values stored in experiment slot
setGeneric(name = "getExperiment",
           def = function(object, model.id=NULL, batchName=NULL, expDig=NULL,
                          treatment.only=FALSE, max.time=NULL,
                          return.list = FALSE, impute.value=FALSE)
           {standardGeneric("getExperiment")} )

#' @export
setMethod( f=getExperiment,
           signature="XevaSet",
           definition=function(object, model.id=NULL, batchName=NULL, expDig=NULL,
                               treatment.only=FALSE, max.time=NULL,
                               return.list = FALSE, impute.value=FALSE)
           {
             if(is.null(model.id) & is.null(batchName) & is.null(expDig))
             {
               msg = sprintf("'model.id' 'batchName' and 'expDig' all NULL")
               stop(msg)
             }

             if(!is.null(model.id))
             {
               mids <- unique(c(model.id))
               rtz <- .getExperimentMultipalIDs(object, mids=mids,
                                                treatment.only=treatment.only,
                                                max.time=max.time,
                                                return.list = return.list,
                                                impute.value=impute.value)

               if(!is.null(max.time))
               { rtz <- rtz[rtz$time<=max.time, ] }
             }

             if(is.null(model.id))
             {
               rtz <- .getBatchData(object, batchName=batchName, expDig=expDig,
                                    treatment.only=treatment.only, max.time=max.time,
                                    return.list = return.list,
                                    impute.value=impute.value)
             }
             return(rtz)
           })


###-----------------------------------------------------------------------------
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
## collapse time-vol data based on expDesign
.collapseRplicate <- function(inLst, var = "volume")
{
  if(is.null(names(inLst))){names(inLst) <- sapply(inLst, function(x){ x$model.id[1]})}
  # if(impute.value==TRUE)
  # {
  #   inLst2 <- list()
  #   timeVec <- sort(unique(unlist(lapply(inLst, "[[", "time" ))))
  #   for(mid in names(inLst))
  #   {
  #     dw <- inLst[[mid]]
  #     inLst2[[mid]] <- .smoothModel(dw, timeVec=timeVec, var =var)
  #   }
  #   inLst <- inLst2
  # }

  timeAll <- sort(unique(unlist(lapply(inLst, "[[", "time" ))))
  rd <- data.frame()
  for(t in timeAll)
  {
    vx <- unlist(sapply(inLst, function(x){ x[x$time==t, var]}))
    vx <- vx[!is.na(vx)]
    vz <- as.list(Rmisc::STDERR(vx))
    rd <- rbind(rd, data.frame(time=t, mean=vz$mean, upper= vz$upper, lower= vz$lower))
  }
  return(rd)
}



#' Get time vs volume data with standard error
#'
#' Given a batch (treatment and control model ids)
#' it will return a data.fram with time vs volume (or any other variable)
#' with standard error calculated. Note that this function do not check if
#' model.id in given batch belongs to same patient
#'
#' @examples
#' data(pdxe)
#' # creat a experiment desing
#' ExpDesign = list(batch.name="myBatch", treatment=c("X.010.BG98", "X.010.BG98"), control=c("X.010.uned"))
#' df = getTimeVarData(object=pdxe, ExpDesign, var = "volume")
#'
#' @param object The \code{Xeva} dataset
#' @param ExpDesign A list with batch.name, treatment and control
#' @param var Name of the variable, default \code{volume}
#' @param drug.name \code{FALSE}. If \code{TRUE} will return drug name also
#' @return a \code{data.fram} with treatment, control and batch.name
#' @export
getTimeVarData <- function(object, batchName=NULL, expDig=NULL, var = "volume",
                           treatment.only=FALSE, drug.name=FALSE,
                           vol.normal=FALSE, impute.value=TRUE, max.time=NULL)
{

  dl <- .getBatchData(object, batchName=batchName, expDig=expDig,
                       treatment.only=treatment.only, max.time=max.time,
                       return.list=TRUE, impute.value=impute.value)

  df <- data.frame()
  if(!is.null(dl$control))
  {
    dfc <- .collapseRplicate(dl$control, var = var)
    dfc$exp.type <- "control"

    if(drug.name==TRUE)
    {
      drugAll <- sort(unique(unlist(lapply(dl$control, "[[", "drug.join.name"))))
      if(length(drugAll)>1)
      { warning("multipal drugs for batch, will colleps by ;") }
      dfc$drug.name <- paste(drugAll, collapse = ";")
    }

    df <- rbind(df, dfc)
  }

  if(!is.null(dl$treatment))
  {
    dft <- .collapseRplicate(dl$treatment, var = var)
    dft$exp.type <- "treatment"

    if(drug.name==TRUE)
    {
      drugAll <- sort(unique(unlist(lapply(dl$treatment, "[[", "drug.join.name"))))
      if(length(drugAll)>1)
      { warning("multipal drugs for batch, will colleps by ;") }
      dft$drug.name <- paste(drugAll, collapse = ";")
    }

    df <- rbind(df, dft)
  }
  return(df)
}




