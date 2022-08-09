#' Get molecular profiles from a XevaSet object
#'
#' This function serves to get molecular profiles from a \code{XevaSet} object.
#'
#' @param object The \code{XevaSet}.
#' @param data.type \code{character}, where one of the molecular data types is needed.
#' @return An \code{ExpressionSet} where sample names are the \code{biobase.id} of the model.
#' @examples
#' data(brca)
#' brca.RNA <- getMolecularProfiles(brca, data.type="RNASeq")
#' @export
getMolecularProfiles <- function(object, data.type)
{
  if(is.element(data.type, names(slot(object, "molecularProfiles")))==FALSE)
  {
    msg = sprintf("available molecular data are\n%s\n",
                  paste(names(object@molecularProfiles), collapse ="\n"))
    stop(msg)
  }
  expset <- slot(object, "molecularProfiles")[[data.type]]
  return(expset)
}

.modelID2biobaseID <- function(object, mDataType, drug=NULL, tissue=NULL, unique.model=TRUE)
{
  modIn <- modelInfo(object, mDataType = mDataType)
  if(unique.model==TRUE)
  {
    modIn <- modIn[unique(modIn$model.id), ]
  }

  if(!is.null(drug))
  {
    modIn <- modIn[modIn$drug %in% c(drug), ]
    if(nrow(modIn)==0)
    {
      msg <- sprintf("No model present with drug %s ", drug)
      stop(msg)
    }
  }

  if(!is.null(tissue))
  {
    modIn <- modIn[modIn$tissue %in% c(tissue), ]
    if(nrow(modIn)==0)
    {
      msg <- sprintf("No model present with tissue %s ", tissue)
      stop(msg)
    }
  }

  bioName <- sprintf("biobase.id.%s", mDataType)
  modIn <- modIn[ !is.na(modIn[, bioName]), ]
  if(nrow(modIn)==0)
  {
    msg <- sprintf("No model present for drug %s with molecular data type %s", drug, mDataType)
    if(!is.null(tissue))
    {
      msg <- sprintf("No model present for drug %s and tissue %s with molecular data type %s",
                     drug, tissue, mDataType)
    }
    stop(msg)
  }
  return(list(data=modIn, bioName=bioName))
}

#####================= Summarize Molecular Profiles ==================
.summarizeMolecularANDSens <- function(object, drug, mDataType, tissue=NULL,
          sensitivity.measure="all", sen.type=c("batch", "model"),
          unique.model=TRUE, batch=NULL)
{
  sen.type <- match.arg(sen.type)

  if(sensitivity.measure=="all")
  {
    if(sen.type=="batch")
    {
      sm <- colnames(slot(object, "sensitivity")$batch)
      sensitivity.measure <- sm[sm!="batch.name"]
      if(length(sensitivity.measure) < 1)
      { stop("batch sensitivity not avaliable") }
    } else {
      sm <- colnames(slot(object, "sensitivity")$model)
      sensitivity.measure <- sm[sm!="model.id"]
      if(length(sensitivity.measure) < 1)
      { stop("model sensitivity not avaliable") }
    }
  } else
    {
    ms <- sensitivity.measure %in% colnames(slot(object, "sensitivity")$model)
    bs <- sensitivity.measure %in% colnames(slot(object, "sensitivity")$batch)

    sen.type <- NA
    if(sum(ms)==length(sensitivity.measure)){sen.type <- "model" }
    if(sum(bs)==length(sensitivity.measure)){sen.type <- "batch" }

    if(is.na(sen.type))
    {
      if(sum(c(ms, bs))==0)
      {
        msg<-sprintf("sensitivity measure '%s' not present in model or batch sensitivity",
                     sensitivity.measure)
        stop(msg)
      } else {
        msg<-sprintf("All sensitivity measure should be either model or batch type")
        stop(msg)
      }
    }
  }

  senType <- sen.type

  if(senType=="model")
  {
    sm <- sensitivity(object, type = senType, sensitivity.measure = sensitivity.measure)
    if(is.null(sensitivity.measure))
    { sensitivity.measure <- colnames(sm)[!(colnames(sm) %in% c("model.id", "batch.name"))] }

    modInX <- .modelID2biobaseID(object, mDataType, drug=drug, tissue=tissue,
                                 unique.model=unique.model)
    modIn <- modInX$data; bioName <- modInX$bioName

    modIn[,c(sensitivity.measure)] <- sm[modIn$model.id, c(sensitivity.measure)]
    molP <- getMolecularProfiles(object, mDataType)
    molP <- molP[, modIn[, bioName]]
    colnames(molP) <- rownames(modIn)
    rownames(Biobase::pData(molP)) <- rownames(modIn)
    pd <- Biobase::pData(molP)
    for(i in colnames(modIn))
    { pd[,i] <- modIn[,i] }
    Biobase::pData(molP) <- pd
    return(molP)
  }

  if(senType=="batch")
  {
    molP <- getBatchAndMolData(object, mDataType, drug, tissue, unique.model,
                               senType, sensitivity.measure)
    return(molP)
  }
}


#' Summarize molecular profiles and drug response data
#'
#' This function serves to get molecular profiles from a \code{XevaSet} object.
#'
#' @param object The \code{XevaSet}.
#' @param drug Name of the drug.
#' @param mDataType \code{character}, where one of the molecular data types is needed.
#' @param tissue Default \code{NULL} will return all tissue types.
#' @param sensitivity.measure Default \code{'all'} will return all sensitivity measures.
#' @param sen.type Type of sensitivity measure. Options are 'batch' (default) or 'model'
#' @param unique.model Default \code{TRUE} will return only one sequncing ID, in the case where one model ID maps to several sequencing IDs.
#' @param batch Name of the batch. Default \code{NULL}.
#' @param summarizedExp If \code{TRUE} will return SummarizedExperiment class object. Default \code{TRUE} will return ExpressionSet.
#'
#' @return An \code{ExpressionSet} where sample names are \code{model.id} and sensitivity measures will be presented in \code{pData}.
#'
#' @examples
#' data(brca)
#'
#' pacRNA <- summarizeData(brca, drug="paclitaxel", mDataType="RNASeq",
#'                     tissue= "BRCA", sensitivity.measure="mRECIST")
#' print(pacRNA)
#'
#' #to get all batch level response
#' pacRNA <- summarizeData(brca, drug="paclitaxel", mDataType="RNASeq",
#'                     tissue= "BRCA", sensitivity.measure="all")
#' print(pacRNA)
#'
#' #to get all model level response
#' pacRNA <- summarizeData(brca, drug="paclitaxel", mDataType="RNASeq",
#'                     tissue= "BRCA", sensitivity.measure="all",
#'                     sen.type="model")
#' print(pacRNA)
#' @details
#' \itemize{
#' \item {If a sequencing sample belongs to multiple models, this function
#' will create a separate column for each model.}
#' \item {All models without molecular data will be removed from the output object.}
#' }
#' @import SummarizedExperiment
#' @export
setGeneric("summarizeData", function(object, drug, mDataType, tissue=NULL,
                                     sensitivity.measure="all", sen.type=c("batch", "model"),
                                     unique.model=TRUE, batch=NULL,
                                     summarizedExp=FALSE) standardGeneric("summarizeData"))
setMethod('summarizeData',
          signature=signature(object = "XevaSet"),
          definition=function(object, drug, mDataType, tissue=NULL,
                              sensitivity.measure="all", sen.type=c("batch", "model"),
                              unique.model=TRUE, batch=NULL, summarizedExp=FALSE)
          {
            df <- .summarizeMolecularANDSens(object, drug, mDataType, tissue,
                                      sensitivity.measure, sen.type,
                                      unique.model, batch)
            if(summarizedExp==TRUE & is(df, "ExpressionSet"))
            {
              df <- SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(df)
            }
            return(df)
          })


#' @rdname summarizeData
#' @export
summarizeMolecularProfiles <- function(...)
{
  warning("summarizeMolecularProfiles has been deprecated, please use summarizeData ")
  summarizeData(...)
}


.batch2DataFram <- function(object, batchName=NULL, expDig=NULL)
{
  if(is.null(batchName) & is.null(expDig))
  { stop("please provide 'batchName' or 'expDig'") }

  if(!is.null(batchName))
  {
    expDig <- batchInfo(object, batch = batchName)
  }

  bat2mods <- data.frame()
  for(i in seq_along(expDig))
  {
    expdi <- expDig[[i]]

    if(!is.null(expdi$control))
    {
      bat2mods <- rbind(bat2mods, data.frame(batch.name=expdi$batch.name,
                                             model.id=expdi$control,
                                             exp.type="control",
                                             stringsAsFactors = FALSE))
    }
    if(!is.null(expdi$treatment))
    {
      bat2mods <- rbind(bat2mods, data.frame(batch.name=expdi$batch.name,
                                             model.id=expdi$treatment,
                                             exp.type="treatment",
                                             stringsAsFactors = FALSE))
    }
  }
  return(bat2mods)
}
