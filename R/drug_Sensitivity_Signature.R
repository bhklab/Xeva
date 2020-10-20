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
    msg1 <- sprintf("models not present in sensitivity slot:\n%s\n",
                    paste(modNotPresent, collapse = "\n"))
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

  for(I in seq_len(nrow(mdf)))
  {
    bid <- model2bidMap[model2bidMap$model.id==mdf[I, "model.id"], "biobase.id"]
    if(length(bid)==0){ bid <- NA }
    mdf[I,"biobase.id"] <- bid[1]
  }
  mdf <- mdf[ !is.na(mdf[,"biobase.id"]), ]
  mdf <- mdf[ as.character(mdf[,"biobase.id"]) %in% colnames(molData),]
  if(nrow(mdf)==0)
  {
    msg <- sprintf("No '%s' ids are common in molecular data and experimental data",
                   mDataType)
    stop(msg)
  }

  mdfI <- .getSensitivityVal(object, sensitivity.measure, mdf, drug=drug,
                             collapse.by=collapse.by)
  return(mdfI)
}

#' @import methods
#' @import Biobase
.getExpressionSet <- function(tx, y, sensitivity.measure, tissue=NULL)
{
  pd <- data.frame(name=colnames(tx), stringsAsFactors = FALSE)
  rownames(pd) <- as.character(pd$name)
  pd[, sensitivity.measure] <- y
  if(!is.null(tissue))
  { pd$tissue <- tissue }

  eSet <- Biobase::ExpressionSet(assayData=as.matrix(tx),
                        phenoData=new("AnnotatedDataFrame",
                                      data=data.frame(pd)))
  return(eSet)
}

##====== drugSensitivitySig for one drug ==========================
#' get drug sensitivity values
#'
#' @description
#' Given a Xeva object and drug name, this function will return sensitivity values for all the genes/features.
#'
#' @param object The \code{Xeva} dataset.
#' @param drug Name of the drug.
#' @param sensitivity.measure Name of the sensitivity measure.
#' @param mDataType Molecular data type.
#' @param standardize Default \code{SD}. Name of the method to use for data standardization before fitting.
#' @param features Set which molecular data features to use. Default \code{NULL} will use all features.
#' @param fit Association method to use, can be 'lm', 'CI', 'pearson' or 'spearman'. Default \code{lm}.
#' @param nthread number of threads
#' @param verbose Default \code{TRUE} will show information
#'
#' @return A \code{data.frame} with features and values.
#'
#' @examples
#' data(brca)
#' senSig <- drugSensitivitySig(object=brca, drug="tamoxifen",
#'                              mDataType="RNASeq",
#'                              sensitivity.measure="slope", fit = "lm")
#'
#' ## example to compute the Pearson correlation between gene expression and PDX response
#' senSig <- drugSensitivitySig(object=brca, drug="tamoxifen",
#'                              mDataType="RNASeq",
#'                              sensitivity.measure="slope", fit = "pearson",
#'                              features=c(1,2,3,4,5,6))
#'
#' @details Method to compute association can be specified by \code{fit}. It can be one of the:
#' \itemize{
#' \item{"lm" for linear models}
#' \item{"CI" for concordance index}
#' \item{"pearson" for Pearson correlation}
#' \item{"spearman" for Spearman correlation}
#' }
#'
#' If fit is set to NA, processed data (an ExpressionSet) will be returned.
#'
#' A matrix of values can be directly passed to molData.
#' In case where a \code{model.id} maps to multiple \code{biobase.id}s, the first \code{biobase.id} in the \code{data.frame} will be used.
#'
#' @export
#' @import Biobase
drugSensitivitySig <- function(object, drug, sensitivity.measure,
                               mDataType,
                               standardize=c("SD", "rescale", "none"),
                               features=NULL,
                               fit = c("lm", "CI", "pearson", "spearman"),
                               nthread=1, verbose=TRUE)
{
  df <- summarizeData(object, drug=drug, mDataType=mDataType,
                      sensitivity.measure=sensitivity.measure)

  if(!is.null(features))
  {
    df <- df[features, ]
  }

  sm <- pData(df)[, sensitivity.measure]
  nonNAsample <- sampleNames(df)[!is.na(sm)]

  if(length(nonNAsample)==0)
  { stop("no sample with non NA sensitivity.measure") }

  if(verbose==TRUE)
  {
    cat(sprintf("Running for %d samples with non NA sensitivity.measure\n",
                length(nonNAsample)))
  }


  df <- df[, nonNAsample]
  x <- t(exprs(df))
  allxcol <- ncol(x)
  x <- removeZeroVar(x, varCutoff=0, sort=FALSE)
  
  y <- pData(df)[, sensitivity.measure]
  
  if(ncol(x)==0)
  { stop("no non-zero variance feature left in the data") }
  
  if(verbose==TRUE)
  {
    if(allxcol-ncol(x) >0)
    {
      txt <- sprintf("Removing %d features because of zero variance\n", allxcol-ncol(x))
      cat(txt)
    }
    
    txt <- sprintf("Running for %d samples with non NA sensitivity.measure and %d features\n",
                   nrow(x), ncol(x))
    cat(txt)
  }
  
  rtx <-compute_association(x, y, fit = fit[1], nthread= nthread,
                            standardize=standardize[1], verbose=verbose)

  rtx$drug <- drug
  rtx <- .reorderCol(rtx, "drug", 2)
  rownames(rtx) <- NULL
  return(rtx)
}


