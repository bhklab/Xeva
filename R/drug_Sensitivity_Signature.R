

.getSensitivityVal <- function(object, sensitivity.measure, mdf, drug, collapse.by="mean")
{
  senType <- "model"
  if(is.element(sensitivity.measure, colnames(object@sensitivity$model))==FALSE)
  {
    msg1 <- sprintf("sensitivity.measure '%s' not present", sensitivity.measure)
    stop(msg1)
  }

  mdfI <- mdf[mdf$drug==drug,]
  modNotPresent <- setdiff(mdfI$model.id, rownames(sensitivity(object, senType)))
  if(length(modNotPresent)>0)
  {
    msg1 <- sprintf("models not present in sensitivity slot:\n%s\n", paste(modNotPresent, collapse = "\n"))
    warning(msg1)
    mdfI <- mdfI[!(mdfI$model.id %in% modNotPresent), ]
  }

  mdfI[,sensitivity.measure] <- sensitivity(object, senType)[mdfI$model.id, sensitivity.measure]

  dupBID <- mdfI$biobase.id[duplicated(mdfI$biobase.id)]
  if(length(dupBID)>0)
  {
    dupDF <- mdfI[mdfI$biobase.id %in%dupBID,]
    msg1 <- sprintf("model.ids have same 'biobase.id'\n%s", printAndCapture(dupDF))
    warning(msg1)

    msg2 <- sprintf("collapsing same 'biobase.id' using %s", collapse.by)
    warning(msg2)
    stop("code not done")
  }

  return(mdfI)
}



.getBioIdSensitivityDF <- function(object, molData, drug, sensitivity.measure,
                                   collapse.by="mean", model.ids, mDataType,
                                   model2bidMap)
{
  mdf <- modelInfo(object)
  if(!is.null(model.ids))
  {
    mdf <- mdf[mdf$model.id %in% model.ids, ]
    if(nrow(mdf)==0)
    {
      msg <- sprintf("'model.ids' are not present in Xeva object")
      stop(msg)
    }
  }
  mdf[,"biobase.id"] <- NA

  for(I in 1:nrow(mdf))
  {
    bid <- model2bidMap[model2bidMap$model.id==mdf[I, "model.id"], "biobase.id"]
    if(length(bid)==0){ bid <- NA }
    mdf[I,"biobase.id"] <- bid[1]
  }
  mdf <- mdf[ !is.na(mdf[,"biobase.id"]), ]
  mdf <- mdf[ as.character(mdf[,"biobase.id"]) %in% colnames(molData),]
  if(nrow(mdf)==0)
  {
    msg <- sprintf("No '%s' ids are comman in molecular data and experimental data", mDataType)
    stop(msg)
  }

  mdfI <- .getSensitivityVal(object, sensitivity.measure, mdf, drug=drug, collapse.by=collapse.by)
  return(mdfI)
}

##====== drugSensitivitySig for one drug ==========================
#' Shows drug sensitivity values
#'
#' @description
#' Given a Xeva object and drug name, this function will return sensitivity values for all the genes/features.
#'
#' @param object The \code{Xeva} dataset.
#' @param drug Name of the drug.
#' @param mDataType Molecular data type.
#' @param molData External data matrix. Rows as features and columns as samples.
#' @param features Set which molecular data features to use. Default \code{NULL} will use all features.
#' @param model.ids Set which \code{model.id} to use from the dataset. Default \code{NULL} will use all \code{model.id}s.
#' @param model2bidMap A \code{data.frame} with \code{model.id} and \code{biobase.id}. Default \code{NULL} will use internal mapping.
#' @param sensitivity.measure Name of the sensitivity measure.
#' @param fit Default \code{lm}.
#' @param standardize Default \code{SD}. Name of the method to use for data standardization befor fitting.
#' @param nthread number of threads
#' @param tissue tissue type. Default \code{NULL} uses \code{'tissue'} from \code{object}.
#' @param verbose Default \code{TRUE} will show information
#'
#' @return A \code{data.frame} with features and values.
#'
#' @examples
#' data(brca)
#' senSig <- drugSensitivitySig(object=brca, drug="tamoxifen",
#'                              mDataType="RNASeq", features=1:5,
#'                              sensitivity.measure="slope", fit = "lm")
#' @details A matrix of values can be directly passed to molData. \code{fit} can be \code{lm}, \code{maxCor}, or \code{gam}.
#' In case where a \code{model.id} maps to multiple \code{biobase.id}s, the first \code{biobase.id} in the \code{data.frame} will be used.
#'
setGeneric(name = "drugSensitivitySig",
           def = function(object, drug,
                          mDataType=NULL, molData=NULL, features=NULL,
                          model.ids=NULL, model2bidMap = NULL,
                          sensitivity.measure="slope",
                          fit = c("lm"),
                          standardize=c("SD", "rescale", "none"),
                          nthread=1, tissue=NULL, verbose=TRUE)
            {standardGeneric("drugSensitivitySig")}
          )

#' @rdname drugSensitivitySig
#' @export
setMethod(f= "drugSensitivitySig",
          signature=c("XevaSet"),
          definition=function(object, drug,
                               mDataType=NULL, molData=NULL, features=NULL,
                               model.ids=NULL, model2bidMap = NULL,
                               sensitivity.measure="slope",
                               fit = c("lm"),
                               standardize=c("SD", "rescale", "none"),
                               nthread=1, tissue=NULL, verbose=TRUE)
  {
  if(is.null(mDataType)& is.null(molData))
  {
    stop("'mDataType' and 'molData' both can't be NULL ")
  }

  if(is.null(molData))
  {
    molData <- Biobase::exprs(getMolecularProfiles(object, mDataType))
  }

  molData <- as.matrix(molData)

  if(is.null(model2bidMap))
  {model2bidMap <- model2BiobaseIdMap(object, mDataType)}

  rtLx <- list()
  drugIx <- c(drug)[1]

  if(verbose==TRUE){cat(sprintf("Running for drug %s\n\n", drugIx))}
  mdfI <- .getBioIdSensitivityDF(object, molData, drugIx, sensitivity.measure,
                                   collapse.by="mean", model.ids, mDataType,
                                   model2bidMap)

  if(nrow(mdfI)<2)
  {
      msg <- sprintf("Too few samples for drug %s\nNumber of samples %d", drugIx, nrow(mdfI))
      stop(msg)
  }

  if(is.null(features))
  { features <- rownames(molData)}

  ##---------------------------------------------------------------
  if(!is.null(tissue))
  {
    if(length(tissue) == 1)
    {
      cat(sprintf("setting 'tissue' = %s for all models", tissue[1]))
      tt <- rep(tissue[1], nrow(mdfI))
      names(tt) <- mdfI$model.id
    }

  if(length(tissue) > 1)
  {
      if(length(tissue)!= nrow(mdfI))
      {stop("length of type should be equeal to length of models")}

      if(is.null(names(tissue)))
      {
        msg <- sprintf("'tissue' have no names. Plese provide a named list")
        stop(msg)
      }

    tt <- tissue
    }
  } else
  {
    if("tissue" %in% colnames(modelInfo(object)))
    {
      typeDF <- mapModelSlotIds(object, id=mdfI$model.id, id.name = "model.id",
                                map.to = "tissue", unique = FALSE)
      tt <- typeDF[, "tissue"]
      names(tt) <- typeDF$model.id
    } else
    {
      warning("'tissue' not present in modelInfo, setting tissue = 'tumor' for all models")
      tt <- rep("tumor", nrow(mdfI))
      names(tt) <- mdfI$model.id
    }
  }

  mdfI[, "tissue"] <- tt[mdfI$model.id]
  ## -------------------------------------------------------------------------

  x <- t(molData[features, mdfI$biobase.id])
  x <- removeZeroVar(x, varCutoff=0, sort=FALSE)

  fetDiff <- ncol(t(molData[features, mdfI$biobase.id])) - ncol(x)
  if(fetDiff>0)
  {
    msg1 <- sprintf("%d features removed because of 0 variance", fetDiff)
    warning(msg1)
  }

  rtx <- .runFit(x = x,
                 y = mdfI[,sensitivity.measure],
                 fit = fit[1],
                 nthread= nthread, type=mdfI[, "tissue"],
                 standardize=standardize[1], verbose=verbose)

  rtx$drug <- drugIx
  rownames(rtx) <- NULL
  rtx <- .reorderCol(rtx, "drug", 2)
  rtLx<- .appendToList(rtLx, rtx)
  rtDF <- do.call("rbind", rtLx)
  rownames(rtDF) <- NULL
  return(rtDF)
  })


####-------------------------------------------------------------------------

#' @import parallel
#' @import doSNOW
#' @import foreach
#' @importFrom stats setNames
.runFit <- function(x, y, fit = c("lm"), nthread=1,
                    type=NULL, standardize='SD', verbose=TRUE)
{
  fit = fit[1]
  #if(class(x)!= "matrix")
  if(is(x, "matrix")==FALSE)
  { stop("x must be a matrix") }

  if(is.null(colnames(x)))
  { colnames(x) <- 1:ncol(x) }

  if(standardize=="SD"){ x <- scale(x)[,]}
  if(standardize=="rescale"){ x <- as.matrix(apply(x,2, .normalize01))}

  ##----------------------------------------------------------------------
  if(fit=="lm")
  {
    #rr <- PharmacoGx:::rankGeneDrugSensitivity(
    rr <- rankGeneDrugSensitivity(data= x,
                                               drugpheno= y,
                                               type= type, #batch=batch,
                                               single.type=FALSE,
                                               standardize=standardize,
                                               nthread=nthread,
                                               verbose=verbose)
    rr <- data.frame(rr[[1]], stringsAsFactors = FALSE)
    rr$feature <- colnames(x)
    rr <- .reorderCol(rr, "feature", 1)
    return(rr)
  }
}


