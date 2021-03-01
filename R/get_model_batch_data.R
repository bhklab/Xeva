getBatchFormatted <- function(object, batch=NULL, patient.id=NULL, drug=NULL,
                              control.name=NULL)
{
  if(all(c(is.null(batch), is.null(patient.id), is.null(drug), is.null(control.name))))
  {  stop("all variables NULL") }

  if(!is.null(batch))
  {
    if(is.character(batch))
    {
      bt <- slot(object, "expDesign")[[batch]]
      if(is.null(bt))
      {
        msg <- sprintf("batch name %s not present in object", batch)
        stop(msg)
      }
      return(bt)
    }

    if(is.list(batch))
    {
      if(!"batch.name" %in% names(batch))
      { stop(sprintf("'batch.name' is required")) }

      if(is.null(batch$treatment) & is.null(batch$control))
      { stop(sprintf("'treatment' and 'control' both can't be NULL")) }

      allMod <- c(batch$treatment, batch$control)
      modNotPr <- setdiff(allMod, rownames(modelInfo(object)))
      if(length(modNotPr)>0)
      {
        stop(sprintf("model.id not present in dataset: %s", paste(modNotPr, collapse = ", ")))
      }

      return(batch)
    }
  } else {
    mid <- modelInfo(object)
    mid <- mid[mid$patient.id==patient.id, ]
    rtx <- list(name=patient.id)
    rtx$treatment <- mid[mid$drug==drug, "model.id"]
    if(!is.null(control.name))
    { rtx$control <- mid[mid$drug==control.name, "model.id"]}
    return(rtx)
  }
}


.normalizeByElement1 <- function(x) 
{
  if(x[1]==0)
  {x <- x + 1}
  return((x - x[1])/x[1])
}

.normalizeVolume <- function(time, volume, atT0)
{
  if(atT0==TRUE)
  { min.posTime <- min(time[time>=0]) }

  if(atT0==FALSE)
  { min.posTime <- min(time) }

  nrv <- volume[time==min.posTime]
  if(nrv == 0)
  {
    warning("treatment start volume zero, adding 1 to compute volume.normal")
    nrv <- nrv + 1
    volume <- volume+1
  }
  vn <- (volume - nrv) / nrv
  return(vn)
}

#' @import grDevices
#' @import stats
#' @import utils
.smoothCurve <- function(trDF, tsDF, x, y, approx.method = "approx")
{
  lmdf <- data.frame(X = trDF[, x], Y = trDF[, y])
  prdf <- data.frame(X = tsDF[, x])
  if (approx.method == "lm")
  {
    lmod <- lm(Y ~ X, data = lmdf)
    tsDF[, y] <- predict(lmod, prdf, se.fit = FALSE)
  }

  if (approx.method == "approx")
  {
    prdValue <- approx(lmdf$X, lmdf$Y, xout = prdf$X)
    tsDF[, y] <- prdValue$y
  }
  return(tsDF)
}

.smoothModel <- function(dw, timeVec, var = "volume")
{
  dw$impute.value <- "NO"
  t2imp <- setdiff(timeVec, dw$time)
  t2imp <- t2imp[t2imp <= max(dw$time)]
  if (length(t2imp) > 0)
  {
    nwRws <-
      data.frame(matrix(
        data = NA,
        nrow = length(t2imp),
        ncol = ncol(dw)
      ),
      stringsAsFactors = FALSE)
    colnames(nwRws) <- colnames(dw)
    nwRws$time <- t2imp
    nwRws$model.id <- dw$model.id[1]
    nwRws$drug.join.name <- dw$drug.join.name[1]
    nwRws$impute.value <- "YES"
    nwRws <- .smoothCurve(dw, nwRws, x = "time", y = var)
    dt <- rbind(dw, nwRws)
    dt <- BBmisc::sortByCol(dt, c("time"))
    return(dt)
  }
  return(dw)
}

.getExperimentDataFromAExpID <-
  function(object, model.id, treatment.only)
  {
    mod <- slot(object, "experiment")[[model.id]]
    if (is.null(mod))
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

    if (treatment.only == TRUE & !is.null(mod.data$dose))
    {
      #tretIndx <-
      #  extractBetweenTags(mod.data$dose, start.tag = 0, end.tag = 0)
      #mod.data <- mod.data[tretIndx,]
      mod.data <- mod.data[mod.data$time>=0, ]
    }

    mod.data$volume.normal <- .normalizeVolume(mod.data$time, mod.data$volume,
                                               atT0=TRUE)

    return(mod.data)
  }

.getExperimentMultipalIDs <-
  function(object,
           mids,
           treatment.only = TRUE,
           max.time = NULL,
           return.list = TRUE,
           impute.value = FALSE,
           vol.normal = FALSE,
           log.volume=FALSE,
           var = "volume")
  {
    rtx <- list()
    for (i in mids)
    {
      miD <-
        .getExperimentDataFromAExpID(object, model.id = i, treatment.only = treatment.only)
      rtx[[i]] <- miD
    }

    if (length(rtx) > 0)
    {
      if (impute.value == TRUE)
      {
        inLst2 <- list()
        timeVec <- sort(unique(unlist(lapply(rtx, "[[", "time"))))
        for (mid in names(rtx))
        {
          inLst2[[mid]] <- .smoothModel(rtx[[mid]], timeVec = timeVec, var = var)
          inLst2[[mid]]$volume.normal <- .normalizeVolume( inLst2[[mid]]$time,
                                                           inLst2[[mid]]$volume,
                                                           atT0 = TRUE)
        }
        rtx <- inLst2
      }

      for (i in names(rtx))
      {
        miD <- rtx[[i]]

        if(log.volume==TRUE)
        {
          miD$volume <- log(miD$volume+1)
          miD$volume.normal <- .normalizeVolume(miD$time, miD$volume,
                                                atT0 = TRUE)
        }

        if (vol.normal == TRUE)
        {
          miD$volume.raw <- miD$volume
          miD$volume <- .normalizeVolume(miD$time, miD$volume, atT0 = TRUE)
        }

        if (vol.normal == "all")
        {
          miD$volume.raw <- miD$volume
          miD$volume <- .normalizeVolume(miD$time, miD$volume, atT0 = FALSE)
        }

        if (!is.null(max.time))
        {
          miD <- miD[miD$time <= max.time,]
        }

        rtx[[i]] <- miD
      }

      if (return.list == FALSE)
      {
        rtx <- .rbindListOfDataframs(rtx)
      }
    }
    return(rtx)
  }


.reformatExpDesig <- function(expDig)
{
  if (any(c("batch.name", "treatment", "control") %in% names(expDig)))
  {
    expDig <- list(expDig)
  }
  names(expDig) <- unlist(lapply(expDig, "[[", "batch.name"))

  return(expDig)
}

.collapseRplicate <- function(inLst, var = "volume")
{
  if (is.null(names(inLst)))
  {
    names(inLst) <- vapply(inLst, function(x)
      {x$model.id[1]}, FUN.VALUE = character(1))
  }
  timeAll <- sort(unique(unlist(lapply(
    inLst, "[[", "time"
  ))))
  rd <- data.frame()
  for (t in timeAll)
  {
    vx <- unlist(lapply(inLst, function(x) {x[x$time == t, var] }))
    vx <- vx[!is.na(vx)]
    vz <- as.list(Rmisc::STDERR(vx))
    rd <-
      rbind(rd,
            data.frame(
              time = t,
              mean = vz$mean,
              upper = vz$upper,
              lower = vz$lower
            ))
  }
  return(rd)
}

.getTimeVarData <- function(dfp,
                            drug.name = TRUE,
                            var = "volume")
{
  df <- data.frame()
  if (!is.null(dfp$control))
  {
    dfc <- .collapseRplicate(dfp$control, var = var)
    dfc$exp.type <- "control"

    if (drug.name == TRUE)
    {
      drugAll <-
        sort(unique(unlist(
          lapply(dfp$control, "[[", "drug.join.name")
        )))
      if (length(drugAll) > 1)
      {
        txt <-
          sprintf(
            "multiple drugs for batch (in control arm), will collapse by ;\nDrugs are %s\n",
            paste0(drugAll, collapse = ",")
          )
        warning(txt)
      }
      dfc$drug.name <- paste(drugAll, collapse = ";")
    }
    df <- rbind(df, dfc)
  }

  if (!is.null(dfp$treatment))
  {
    dft <- .collapseRplicate(dfp$treatment, var = var)
    dft$exp.type <- "treatment"

    if (drug.name == TRUE)
    {
      drugAll <-
        sort(unique(unlist(
          lapply(dfp$treatment, "[[", "drug.join.name")
        )))
      if (length(drugAll) > 1)
      {
        txt <-
          sprintf(
            "multiple drugs for batch (in treatment arm), will collapse by ;\nDrugs are %s\n",
            paste0(drugAll, collapse = ",")
          )
        warning(txt)
      }
      dft$drug.name <- paste(drugAll, collapse = ";")
    }

    df <- rbind(df, dft)
  }
  return(df)
}

.getBatchData <-
  function(object,
           batch = NULL,
           patient.id = NULL,
           drug = NULL,
           control.name = NULL,
           treatment.only = FALSE,
           max.time = NULL,
           return.list = TRUE,
           impute.value = FALSE,
           vol.normal = FALSE,
           log.volume =FALSE,
           drug.name = TRUE,
           concurrent.time = FALSE)
  {
    expDig <- getBatchFormatted(object, batch, patient.id, drug, control.name)
    expDig <- .reformatExpDesig(expDig)[[1]]

    dfp <- list()
    if (length(expDig$control) > 0)
    {
      dfp$control <-
        .getExperimentMultipalIDs(
          object,
          mids = expDig$control,
          treatment.only = treatment.only,
          max.time = max.time,
          return.list = TRUE,
          impute.value = impute.value,
          vol.normal = vol.normal,
          log.volume =log.volume
        )
      for (ci in seq_along(dfp$control))
      {
        dfp$control[[ci]]$exp.type <- "control"
      }
    }

    if (length(expDig$treatment) > 0)
    {
      dfp$treatment <-
        .getExperimentMultipalIDs(
          object,
          mids = expDig$treatment,
          treatment.only = treatment.only,
          max.time = max.time,
          return.list = TRUE,
          impute.value = impute.value,
          vol.normal = vol.normal,
          log.volume =log.volume
        )
      for (ti in seq_along(dfp$treatment))
      {
        dfp$treatment[[ti]]$exp.type <- "treatment"
      }
    }

    dfp$batch <- .getTimeVarData(dfp, drug.name = drug.name)

    if (concurrent.time == TRUE)
    {
      mc <- mt <- NA
      if ("control" %in% dfp$batch$exp.type)
      {
        mc <- max(dfp$batch[dfp$batch$exp.type == "control", "time"])
      }
      if ("treatment" %in% dfp$batch$exp.type)
      {
        mt <- max(dfp$batch[dfp$batch$exp.type == "treatment", "time"])
      }

      tx <- c(mc, mt)
      tx <- tx[!is.na(tx)]
      if (length(tx) > 1)
      {
        dfp$batch <- dfp$batch[dfp$batch$time <= min(tx),]
        if (!is.null(dfp$control))
        {
          for (ci in seq_along(dfp$control))
          {
            dfp$control[[ci]] <-
              dfp$control[[ci]][dfp$control[[ci]]$time <= min(tx), ]
          }
        }

        if (!is.null(dfp$treatment))
        {
          for (ci in seq_along(dfp$treatment))
          {
            dfp$treatment[[ci]] <-
              dfp$treatment[[ci]][dfp$treatment[[ci]]$time <= min(tx), ]
          }
        }
      }
    }

    if (return.list == TRUE)
    {
      rtx <- dfp
    } else
    {
      rtx <- list(batch = dfp$batch)
      dfcnt <- data.frame()
      dftre <- data.frame()

      if (length(dfp$control) > 0)
      {
        dfcnt <- .rbindListOfDataframs(dfp$control)
      }

      if (length(dfp$treatment) > 0)
      {
        dftre <- .rbindListOfDataframs(dfp$treatment)
      }

      if (nrow(dfcnt) > 0 & nrow(dftre) > 0)
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

##----- get experiment data in flat data.fram ----------------------------------
#' Get PDX experiment data
#'
#' For a given \code{model.id}, \code{getExperiment} will
#'
#' @param object The \code{XevaSet} object.
#' @param model.id The \code{model.id} for which data is required, multiple IDs are allowed.
#' @param batch Batch name from the \code{XevaSet} or experiment design.
#' @param patient.id Patient id from the \code{XevaSet}. Default \code{NULL}.
#' @param drug Name of the drug.
#' @param control.name Name of drug used as control. Default \code{NULL}.
#' @param treatment.only Default \code{FALSE}. If \code{TRUE}, give data for time>=0 only.
#' @param max.time Maximum time for data.
#' @param vol.normal Default \code{FALSE} will return raw volume. If TRUE it will normalize the volume by treatment start volume (volume at time>=0). If set to \code{'all'} will use full data to volume normalization.
#' @param log.volume If TRUE log of the volume will be used. Default \code{FALSE}.
#' @param return.list Default \code{FALSE} will return a \code{data.frame}.
#' @param impute.value Default \code{FALSE}. If \code{TRUE}, impute the missing values.
#' @param concurrent.time Default \code{FALSE}. If \code{TRUE}, cut the batch data such that control and treatment will end at same time point.
#'
#' @examples
#' data(brca)
#'
#' getExperiment(brca, model.id="X.6047.uned", treatment.only=TRUE)
#'
#' getExperiment(brca, model.id=c("X.6047.uned", "X.6047.pael"), treatment.only=TRUE)
#'
#' getExperiment(brca, batch="X-6047.paclitaxel", treatment.only=TRUE)
#'
#' ed <- list(batch.name="myBatch", treatment=c("X.6047.LJ16","X.6047.LJ16.trab"),
#'              control=c("X.6047.uned"))
#'
#' getExperiment(brca, batch=ed)
#'
#' @return a \code{data.fram} will all the the values stored in experiment slot
setGeneric(
  name = "getExperiment",
  def = function(object,
                 model.id = NULL,
                 batch = NULL,
                 patient.id = NULL,
                 drug = NULL,
                 control.name = NULL,
                 treatment.only = FALSE,
                 max.time = NULL,
                 vol.normal = FALSE,
                 log.volume = FALSE,
                 return.list = FALSE,
                 impute.value = FALSE,
                 concurrent.time = FALSE)
  {
    standardGeneric("getExperiment")
  }
)

#' @rdname getExperiment
#' @export
setMethod(
  f = getExperiment,
  signature = "XevaSet",
  definition = function(object,
                        model.id = NULL,
                        batch = NULL,
                        patient.id = NULL,
                        drug = NULL,
                        control.name = NULL,
                        treatment.only = FALSE,
                        max.time = NULL,
                        vol.normal = FALSE,
                        log.volume = FALSE,
                        return.list = FALSE,
                        impute.value = FALSE,
                        concurrent.time = FALSE)
  {
    if (is.null(model.id) & is.null(batch) & is.null(patient.id))
    {
      msg <- sprintf("'model.id' 'batch' and 'patient.id' all NULL")
      stop(msg)
    }

    if(!as.character(vol.normal)%in% c("TRUE","FALSE","all"))
    {
      msg <- sprintf("allowed 'vol.normal' values are: TRUE, FALSE,'all' ")
      stop(msg)
    }

    if (!is.null(model.id))
    {
      mids <- unique(c(model.id))
      rtz <- .getExperimentMultipalIDs(
        object,
        mids = mids,
        treatment.only = treatment.only,
        max.time = max.time,
        return.list = return.list,
        impute.value = impute.value,
        vol.normal = vol.normal,
        log.volume =log.volume
      )

      if (!is.null(max.time))
      {
        rtz <- rtz[rtz$time <= max.time,]
      }
    }

    if (!is.null(batch) | !is.null(patient.id))
    {
      rtz <- .getBatchData(
        object,
        batch = batch,
        patient.id = patient.id,
        drug = drug,
        control.name = control.name,
        treatment.only = treatment.only,
        max.time = max.time,
        return.list = return.list,
        impute.value = impute.value,
        vol.normal = vol.normal,
        log.volume =log.volume,
        concurrent.time = concurrent.time
      )
    }
    return(rtz)
  }
)
