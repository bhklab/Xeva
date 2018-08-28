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

  mod.data <- slot(mod, "data")
  mod.data$model.id <- slot(mod, "model.id")
  mod.data$drug.join.name <- slot(mod, "drug")$join.name


  mod.data <- .removeNAcol(mod.data)
  mod.data <- .reorderCol(mod.data, "model.id", 1)
  mod.data <- .reorderCol(mod.data, "drug.join.name", 2)

  if(treatment.only==TRUE & !is.null(mod.data$dose))
  {
    tretIndx <- extractBetweenTags(mod.data$dose, start.tag=0, end.tag=0)
    mod.data <- mod.data[tretIndx, ]
  }

  mod.data$volume.normal <- NA
  mod.data$volume.normal <- .normalizeVolume(mod.data$volume)
  return(mod.data)
}

.getExperimentMultipalIDs <- function(object, mids, treatment.only=TRUE,
                                      max.time=NULL, return.list=TRUE,
                                      impute.value=FALSE, vol.normal=FALSE,
                                      var="volume")
{
  rtx <- list()
  for(i in mids)
  {
    miD <- .getExperimentDataFromAExpID(object, model.id=i, treatment.only=treatment.only)
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

    for(i in names(rtx))
    {
      miD <- rtx[[i]]
      if(vol.normal==TRUE)
      {
        miD$volume.raw <- miD$volume
        miD$volume <- miD$volume.normal
        miD$volume.normal <- NULL
      }

      if(!is.null(max.time))
      { miD <- miD[miD$time<=max.time, ]}

      rtx[[i]]<- miD
    }

    if(return.list==FALSE)
    { rtx <- .rbindListOfDataframs(rtx) }
  }
  return(rtx)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

.reformatExpDesig <- function(expDig)
{
  if(any(c("batch.name", "treatment", "control") %in% names(expDig)))
  {
    expDig <- list(expDig)
  }
  names(expDig) <- sapply(expDig, "[[", "batch.name")
  return(expDig)
}

## collapse time-vol data based on expDesign
.collapseRplicate <- function(inLst, var = "volume")
{
  if(is.null(names(inLst))){names(inLst) <- sapply(inLst, function(x){ x$model.id[1]})}
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

.getTimeVarData <- function(dfp, drug.name=TRUE, var="volume")
{
  df <- data.frame()
  if(!is.null(dfp$control))
  {
    dfc <- .collapseRplicate(dfp$control, var = var)
    dfc$exp.type <- "control"

    if(drug.name==TRUE)
    {
      drugAll <- sort(unique(unlist(lapply(dfp$control, "[[", "drug.join.name"))))
      if(length(drugAll)>1)
      { warning("multipal drugs for batch, will colleps by ;") }
      dfc$drug.name <- paste(drugAll, collapse = ";")
    }
    df <- rbind(df, dfc)
  }

  if(!is.null(dfp$treatment))
  {
    dft <- .collapseRplicate(dfp$treatment, var = var)
    dft$exp.type <- "treatment"

    if(drug.name==TRUE)
    {
      drugAll <- sort(unique(unlist(lapply(dfp$treatment, "[[", "drug.join.name"))))
      if(length(drugAll)>1)
      { warning("multipal drugs for batch, will colleps by ;") }
      dft$drug.name <- paste(drugAll, collapse = ";")
    }

    df <- rbind(df, dft)
  }
  return(df)
}

.getBatchData <- function(object, batchName=NULL, expDig =NULL, treatment.only=FALSE,
                          max.time=NULL, return.list = TRUE, impute.value=FALSE,
                          vol.normal=FALSE, drug.name=TRUE, concurrent.time=FALSE)
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
                                             return.list = TRUE, #return.list,
                                             impute.value=impute.value,
                                             vol.normal=vol.normal)
    for(ci in 1:length(dfp$control))
    { dfp$control[[ci]]$exp.type <- "control"}
  }

  if(length(expDig$treatment)>0)
  {
    dfp$treatment <- .getExperimentMultipalIDs(object, mids=expDig$treatment,
                                               treatment.only=treatment.only,
                                               max.time=max.time,
                                               return.list = TRUE, #return.list,
                                               impute.value=impute.value,
                                               vol.normal=vol.normal)
    for(ti in 1:length(dfp$treatment))
    { dfp$treatment[[ti]]$exp.type <- "treatment"}
  }

  dfp$batch <- .getTimeVarData(dfp, drug.name = drug.name)

  if(concurrent.time==TRUE)
  {
    mc <- mt <- NA
    if("control" %in% dfp$batch$exp.type)
    { mc <- max(dfp$batch[dfp$batch$exp.type=="control", "time"])}
    if("treatment" %in% dfp$batch$exp.type)
    { mt <- max(dfp$batch[dfp$batch$exp.type=="treatment", "time"])}

    tx <- c(mc, mt); tx <- tx[!is.na(tx)]
    if(length(tx)>1)
    {
      dfp$batch <- dfp$batch[dfp$batch$time<=min(tx), ]
      if(!is.null(dfp$control))
      {
        for(ci in 1:length(dfp$control))
        { dfp$control[[ci]] <- dfp$control[[ci]][dfp$control[[ci]]$time<=min(tx),]}
      }

      if(!is.null(dfp$treatment))
      {
        for(ci in 1:length(dfp$treatment))
        { dfp$treatment[[ci]] <- dfp$treatment[[ci]][dfp$treatment[[ci]]$time<=min(tx),]}
      }
    }
  }

  if(return.list==TRUE)
  { rtx <- dfp } else
  {
    rtx <- list(batch= dfp$batch)
    dfcnt <- data.frame()
    dftre <- data.frame()

    if(length(dfp$control)>0)
    { dfcnt <- .rbindListOfDataframs(dfp$control) }

    if(length(dfp$treatment)>0)
    { dftre <- .rbindListOfDataframs(dfp$treatment) }

    if(nrow(dfcnt)>0 & nrow(dftre)>0)
    {
      allClNames <- unique(c(colnames(dfcnt), colnames(dftre)))
      dfcnt[, setdiff(allClNames, colnames(dfcnt))] <- NA
      dftre[, setdiff(allClNames, colnames(dftre))] <- NA

      dfcnt <- dfcnt[, allClNames]
      dftre <- dftre[, allClNames]
    }
    rtx$model <- rbind(dfcnt, dftre)
  }

  return(rtx)
}

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##----- get experiment data in flat data.fram ----------------------------------
#' For a given  model.id, it will return a data.fram
#' containing all data stored in experiment slot
#'
#' @param object The \code{XevaSet}
#' @param model.id The \code{model.id} for which data is required, multipal allowed
#' @param batchName batch name from the Xeva set
#' @param expDig Experiment design
#' @param treatment.only Default \code{FALSE}. If TRUE give data only for non-zero dose periode (if dose data avalible)
#' @param max.time maximum time for data
#' @param vol.normal default \code{TRUE} will use
#' @param return.list default \code{FALSE} will return a datafram
#' @param impute.value default \code{FALSE}. If TRUE will impute the values
#' @param concurrent.time default \code{FALSE}. If TRUE will cut the batch data such that control and treatment will end at same time point
#'
#' @examples
#' data(brca)
#' getExperiment(brca, model.id="X.6047.uned", treatment.only=TRUE)
#' getExperiment(brca, model.id=c("X.6047.uned", "X.6047.pael"), treatment.only=TRUE)
#' getExperiment(brca, batchName="X-6047.paclitaxel", treatment.only=TRUE)
#' expDesign <- list(batch.name="myBatch", treatment=c("X.6047.LJ16","X.6047.LJ16.trab"),
#'              control=c("X.6047.uned"))
#' getExperiment(brca, expDig=expDesign)
#'
#' @return a \code{data.fram} will all the the values stored in experiment slot
setGeneric(name = "getExperiment",
           def = function(object, model.id=NULL, batchName=NULL, expDig=NULL,
                          treatment.only=FALSE, max.time=NULL, vol.normal=FALSE,
                          return.list = FALSE, impute.value=FALSE,
                          concurrent.time=FALSE)
           {standardGeneric("getExperiment")} )

#' @export
setMethod( f=getExperiment,
           signature="XevaSet",
           definition=function(object, model.id=NULL, batchName=NULL, expDig=NULL,
                               treatment.only=FALSE, max.time=NULL, vol.normal=FALSE,
                               return.list = FALSE, impute.value=FALSE,
                               concurrent.time=FALSE)
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
                                                impute.value=impute.value,
                                                vol.normal=vol.normal)

               if(!is.null(max.time))
               { rtz <- rtz[rtz$time<=max.time, ] }
             }

             if(is.null(model.id))
             {
               rtz <- .getBatchData(object, batchName=batchName, expDig=expDig,
                                    treatment.only=treatment.only, max.time=max.time,
                                    return.list = return.list,
                                    impute.value=impute.value,
                                    vol.normal=vol.normal,
                                    concurrent.time=concurrent.time)
             }
             return(rtz)
           })
