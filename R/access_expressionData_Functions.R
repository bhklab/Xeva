#####================= getMolecularProfiles ==================
#' Get molecular profiles from a XevaSet object
#'
#' This function serves to get molecular profiles from a \code{XevaSet} object.
#'
#' @param object The \code{XevaSet}.
#' @param data.type \code{character}, where one of the molecular data types is needed.
#' @return A \code{SummarizedExperiment} object (converted from ExpressionSet if needed), where sample names are the \code{biobase.id} of the model.
#' @examples
#' data(brca)
#' brca.RNA <- getMolecularProfiles(brca, data.type="RNASeq")
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @importFrom MultiAssayExperiment experiments
getMolecularProfiles <- function(object, data.type)
{
  if(is.element(data.type, names(slot(object, "molecularProfiles")))==FALSE)
  {
    msg = sprintf("available molecular data are\n%s\n",
                  paste(names(object@molecularProfiles), collapse ="\n"))
    stop(msg)
  }
  molData <- .getMolecularData(object, data.type)

  # Convert ExpressionSet to SummarizedExperiment if needed
  if (inherits(molData, "ExpressionSet")) {
    molData <- SummarizedExperiment::SummarizedExperiment(
      assays = list(exprs = Biobase::exprs(molData)),
      colData = S4Vectors::DataFrame(Biobase::pData(molData)),
      rowData = S4Vectors::DataFrame(Biobase::fData(molData))
    )
    message(sprintf("Note: '%s' molecular data was stored as ExpressionSet and has been converted to SummarizedExperiment.", data.type))
  }

  return(molData)
}

#' Get a specific assay (SummarizedExperiment) from a XevaSet object's molecularProfiles
#'
#' This function provides direct access to a specific assay stored in the
#' \code{MultiAssayExperiment} slot of a \code{XevaSet} object.
#'
#' @param object The \code{XevaSet} object.
#' @param data.type \code{character} name of the assay (e.g. "RNASeq").
#' @return A \code{SummarizedExperiment} corresponding to the given \code{data.type}.
#' @examples
#' data(brca)
#' rna <- getMolecularProfileAssay(brca, "RNASeq")
#' @export
getMolecularProfileAssay <- function(object, data.type) {
  mp <- object@molecularProfiles
  if (inherits(mp, "MultiAssayExperiment")) {
    if (!(data.type %in% names(MultiAssayExperiment::experiments(mp)))) {
      stop(sprintf("Data type '%s' not found in molecularProfiles", data.type))
    }
    return(MultiAssayExperiment::experiments(mp)[[data.type]])
  } else if (is.list(mp)) {
    if (!(data.type %in% names(mp))) {
      stop(sprintf("Data type '%s' not found in molecularProfiles", data.type))
    }
    se <- mp[[data.type]]
    if (inherits(se, "ExpressionSet")) {
      se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(exprs = Biobase::exprs(se)),
        colData = S4Vectors::DataFrame(Biobase::pData(se)),
        rowData = S4Vectors::DataFrame(Biobase::fData(se))
      )
    }
    return(se)
  } else {
    stop("Unsupported format in 'molecularProfiles' slot.")
  }
}

#' Internal helper to retrieve molecular data (SE) from a XevaSet
#'
#' Returns a SummarizedExperiment from the XevaSet@molecularProfiles slot,
#' regardless of whether it's a list or MultiAssayExperiment.
#'
#' @param object The \code{XevaSet} object.
#' @param data.type \code{character}, the name of the molecular data type (e.g. "RNASeq").
#' @return A \code{SummarizedExperiment} object.
#' @keywords internal
.getMolecularData <- function(object, data.type) {
  molSlot <- object@molecularProfiles

  if (inherits(molSlot, "MultiAssayExperiment")) {
    if (!(data.type %in% names(MultiAssayExperiment::experiments(molSlot)))) {
      stop(sprintf("Data type '%s' not found in molecularProfiles", data.type))
    }
    molData <- MultiAssayExperiment::experiments(molSlot)[[data.type]]
  } else if (is.list(molSlot)) {
    if (!(data.type %in% names(molSlot))) {
      stop(sprintf("Data type '%s' not found in molecularProfiles", data.type))
    }
    molData <- molSlot[[data.type]]
  } else {
    stop("Unsupported format in 'molecularProfiles' slot.")
  }

  # Convert ExpressionSet to SummarizedExperiment if needed
  if (inherits(molData, "ExpressionSet")) {
    molData <- SummarizedExperiment::SummarizedExperiment(
      assays = list(exprs = Biobase::exprs(molData)),
      colData = S4Vectors::DataFrame(Biobase::pData(molData)),
      rowData = S4Vectors::DataFrame(Biobase::fData(molData))
    )
  }

  return(molData)
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
#' @return A \code{SummarizedExperiment} object where sample names are \code{model.id}, and sensitivity measures are stored in \code{colData}.
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
#' \item {All models without molecular data will be removed from the output \code{SummarizedExperiment}.}
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
    molP <- .getMolecularData(object, mDataType) # good for backward compatibility
    # molP <- getMolecularProfileAssay(object, mDataType) # good for clarity and performance
    molP <- molP[, modIn[, bioName]]
    # Reassign sample-level metadata using colData
    colnames(molP) <- rownames(modIn)
    rownames(modIn) <- rownames(modIn)  # Ensures proper alignment

    SummarizedExperiment::colData(molP) <- S4Vectors::DataFrame(modIn)

    return(molP)
  }

  if(senType=="batch")
  {
    stop("not implemented yet")
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
