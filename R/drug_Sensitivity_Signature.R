

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
#' drugSensitivitySig
#'
#' @description
#' Given a Xeva object, and drug name it will return sensitivity value for all the genes/fetures
#'
#' @param object The \code{Xeva} dataset
#' @param drug Name of the drug
#' @param mDataType molecular data type
#' @param molData External data matrix. Rows as features and columns as samples
#' @param features which  molecular data fetures to use. Default \code{NULL} will use all fetures
#' @param model.ids which model.id to use from the dataset. Default \code{NULL} will use all model.id
#' @param model2bidMap a datafram with model.id and biobase.id. Default \code{NULL} will use internal mapping
#' @param sensitivity.measure Name of the sensitivity measure
#' @param fit Default \code{lm}. Name of the model to be fitted. Options are "lm", "maxCor", "gam"
#' @param type Tissue type. Default is NULL which will use \code{'tumor.type'} from \code{object}
#' @return A datafram with fetures and values
#'
#' @examples
#' data(pdxe)
#' ## select BRCA samples
#' mid <- modelInfo(pdxe)[modelInfo(pdxe)$tumor.type=="BRCA", ]
#' drugSensitivitySig(object=pdxe, drug=c("paclitaxel","tamoxifen"),
#'                    mDataType="RNASeq", features=1:5,
#'                    model.ids = mid$model.id,
#'                    sensitivity.measure="slope", fit = "lm")
#' @details A matrix of values can be directly passed to molData. \code{fit} can be "lm", "maxCor" or "gam".
#' In case where a model.id map to multipal biobase.id the first biobase.id in the datafram will be used.
#'
setGeneric(name = "drugSensitivitySig",
           def = function(object, drug,
                          mDataType=NULL, molData=NULL, features=NULL,
                          model.ids=NULL, model2bidMap = NULL,
                          sensitivity.measure="slope",
                          fit = c("lm", "maxCor", "gam"),
                          standardize=c("SD", "rescale", "none"),
                          nthread=1, tumor.type=NULL, verbose=TRUE)
            {standardGeneric("drugSensitivitySig")}
          )

#' @export
setMethod(f= "drugSensitivitySig",
          signature=c("XevaSet"),
          definition=function(object, drug,
                               mDataType=NULL, molData=NULL, features=NULL,
                               model.ids=NULL, model2bidMap = NULL,
                               sensitivity.measure="slope",
                               fit = c("lm", "maxCor", "gam"),
                               standardize=c("SD", "rescale", "none"),
                               nthread=1, tumor.type=NULL, verbose=TRUE)
  {
  #molecular.summary.stat=c("mean", "median", "first", "last", "or", "and"),
  #sensitivity.summary.stat=c("mean", "median", "first", "last"),
  #returnValues=c("estimate", "pvalue", "fdr"),
  #sensitivity.cutoff,

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
  for(drugIx in c(drug))
  {
    if(verbose==TRUE){printf("Running for drug %s\n\n", drugIx)}
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
    if(!is.null(tumor.type))
    {
      #
      if(length(tumor.type) == 1)
      {
        printf("setting 'tumor.type' = %s for all models", tumor.type[1])
        tt <- rep(tumor.type[1], nrow(mdfI))
        names(tt) <- mdfI$model.id
      }

      if(length(tumor.type) > 1)
      {
        if(length(tumor.type)!= nrow(mdfI))
        {stop("length of type should be equeal to length of models")}

        if(is.null(names(tumor.type)))
        {
          msg <- sprintf("'tumor.type' have no names. Plese provide a named list")
          stop(msg)
        }

        tt <- tumor.type
      }

    } else
    {
      if("tumor.type" %in% colnames(modelInfo(object)))
      {
        typeDF <- mapModelSlotIds(object, id=mdfI$model.id, id.name = "model.id",
                                  map.to = "tumor.type", unique = FALSE)
        tt <- typeDF[, "tumor.type"]
        names(tt) <- typeDF$model.id

      } else
      {
        warning("'tumor.type' not present in modelInfo, setting tumor.type = 'tumor' for all models")
        tt <- rep("tumor", nrow(mdfI))
        names(tt) <- mdfI$model.id
      }
    }
    mdfI[, "tumor.type"] <- tt[mdfI$model.id]
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
                   nthread= nthread, type=mdfI[, "tumor.type"],
                   standardize=standardize[1], verbose=verbose)

    rtx$drug <- drugIx
    rownames(rtx) <- NULL
    rtx <- .reorderCol(rtx, "drug", 2)
    rtLx<- .appendToList(rtLx, rtx)
  }

  rtDF <- do.call("rbind", rtLx)
  rownames(rtDF) <- NULL
  return(rtDF)
  })


####-------------------------------------------------------------------------

#' @import parallel
#' @import doSNOW
#' @import foreach
.runFit <- function(x, y, fit = c("lm", "maxCor", "gam"), nthread=1,
                    type=NULL, standardize='SD', verbose=TRUE)
{
  fit = fit[1]
  if(class(x)!= "matrix")
  {
    stop("x must be a matrix")
  }

  if(is.null(colnames(x)))
  { colnames(x) <- 1:ncol(x) }

  if(standardize=="SD"){ x <- scale(x)[,]}
  if(standardize=="rescale"){ x <- as.matrix(apply(x,2, .normalize01))}

  ##----------------------------------------------------------------------
  if(fit=="lm")
  {
    rr <- PharmacoGx:::rankGeneDrugSensitivity(data= x,
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

  ##---------------------------------------------------------------------
    #library(doSNOW)
    cl <- makeCluster(nthread)
    registerDoSNOW(cl)

    #result <- foreach (i=colnames(x), .final = function(i) setNames(i, colnames(x))) %dopar%
    #{ .nonLinerFits(x[,i], y, fit = fit )}

    result <- foreach (i=colnames(x),
                       .final = function(i) {setNames(i, colnames(x))},
                       .export=c(".nonLinerFits")) %dopar%
    { .nonLinerFits(x[,i], y, fit = fit )}


    stopCluster(cl)

    if(fit == "maxCor")
    {
      rtx <- data.frame(feature= colnames(x),
                        maxCor = sapply(result, function(i) i[1,1]))
    }

    if(fit == "gam")
    {
      rtx <- .convertListToDataFram(result)
      rtx$feature <- rownames(rtx)
      rtx <- .reorderCol(rtx, "feature", 1)
    }
    return(rtx)
}


.nonLinerFits <- function(x, y, fit)
{
  switch(fit,
         ##--------- Maximal correlation ---------
         "maxCor" = { argmax <- acepack::ace(x, y)
                      value <- stats::cor(argmax$tx, argmax$ty)},
         ##------- generalized additive model  -------------
         "gam" = { g <- mgcv::gam(y ~ s(x))
                   val <- mgcv::summary.gam(g)
                   value <- sapply(c("r.sq", "dev.expl"), function(i) val[[i]]) },
         value <- NA
         )

  return(value)
}



nonLinerFitExample <- function()
{
  x = -100 : 100
  y = (100^2 - x^2)^0.5    # let's make a circle

  cor(x, y) # a circular line has zero linear correlation
  cor(x, y, method='spearman') ## and zero rank correlation

  ##--------- Maximal correlation ---------
  library(acepack)
  argmax = ace(x, y)
  cor(argmax$tx, argmax$ty)

  ##------- generalized additive model  -------------
  library(mgcv)
  g <- gam(y ~ s(x))
  summary(g)
  plot(g,scheme=2)

}






########-----------------------------------------------------
########-----------------------------------------------------
##====== geneSensitivityPlot for one drug ==========================
#' geneSensitivityPlot
#'
#' Plot a gene expression against sensitivity signatures for a drug
#' @description
#' Given a Xeva object, feture name and drug name it will plot feture values against sensitivity value
#'
#' @examples
#' data(cm.pdxe)
#' geneSensitivityPlot(object=cm.pdxe, mDataType="RNASeq", feature="A1BG", drug="binimetinib",
#' sensitivity.measure="slope", standardize="log")
#'
#' @export
#' @import ggplot2
geneSensitivityPlot <- function(object, mDataType, feature, drug,
                                sensitivity.measure="slope",
                                standardize=c("SD", "rescale", "log", "none"))
{
  molData <- Biobase::exprs(getMolecularProfiles(object, mDataType))
  sdf <- .getBioIdSensitivityDF(object, molData, drug, sensitivity.measure, collapse.by="mean" )
  ##---plot gene vs  sensitivity.measure---------------------
  #feature = rownames(exprs(molData))[1]#:5]

  sdf[, feature] <- molData[feature, sdf$biobase.id]
  if(standardize[1]=="SD"){sdf[, feature] <- as.vector(scale(sdf[, feature]))}
  if(standardize[1]=="rescale"){sdf[, feature] <- .normalize01(sdf[, feature])}
  if(standardize[1]=="log")
  {
    if(min(sdf[, feature])<= 0)
    {
      sdf[, feature] <- sdf[, feature] + abs(min(sdf[, feature]))+0.0001
    }
    sdf[, feature] <- log(sdf[, feature])
  }

  plt <- ggplot(sdf, aes_string(x=sensitivity.measure, y=feature))#, color= "type", group="model.id"))
  plt <- plt + geom_line(linetype = 1)+ geom_point()
  .ggplotEmptyTheme(plt)
}









