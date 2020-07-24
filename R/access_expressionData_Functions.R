#####================= getMolecularProfiles ==================
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
#' Summarize molecular profiles
#'
#' This function serves to get molecular profiles from a \code{XevaSet} object.
#'
#' @param object The \code{XevaSet}.
#' @param drug Name of the drug.
#' @param mDataType \code{character}, where one of the molecular data types is needed.
#' @param tissue Default \code{NULL} will return all tissue types.
#' @param sensitivity.measure Default \code{NULL} will return all sensitivity measures.
#' @param unique.model Default \code{TRUE} will return only one sequncing ID, in the case where one model ID maps to several sequencing IDs.
#' @param batch Name of the batch. Default \code{NULL}.
#' @return An \code{ExpressionSet} where sample names are \code{model.id} and sensitivity measures will be presented in \code{pData}.
#'
#' @examples
#' data(brca)
#' pacRNA <- summarizeMolecularProfiles(brca, drug="paclitaxel", mDataType="RNASeq",
#'                                      tissue= "BRCA", sensitivity.measure="mRECIST")
#' print(pacRNA)
#' @details
#' \itemize{
#' \item {If a sequencing sample belongs to multiple models, \code{summarizeMolecularProfiles}
#' will create a separate column for each model.}
#' \item {All models without molecular data will be removed from the output \code{ExpressionSet}.}
#' }
#' @export
summarizeMolecularProfiles <- function(object, drug, mDataType, tissue=NULL,
                                       sensitivity.measure=NULL, unique.model=TRUE,
                                       batch=NULL)
{
  senType <- NULL
  if(is.null(sensitivity.measure))
  { senType <- "model" } else
  {
    if(sensitivity.measure %in% colnames(slot(object, "sensitivity")$model))
    {senType <- "model" }

    if(is.null(senType))
    {
      if(sensitivity.measure %in% colnames(slot(object, "sensitivity")$batch))
      {senType <- "batch" }
    }
  }

  if(is.null(senType))
  {
    msg <- sprintf("sensitivity measure '%s' not present in model or batch sensitivity", sensitivity.measure)
    stop(msg)
  }
  ##----------------------------------------------------------------------------
  ##----------------------------------------------------------------------------
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
