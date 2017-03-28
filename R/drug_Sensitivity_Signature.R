.getModDrugBid <- function(object, drugs=NULL)
{
  midDr <- sapply(object@experiment, "[[", c("drug", "join.name"))
  midDr <- data.frame(model.id= names(midDr), drug= midDr, stringsAsFactors = FALSE)
  if(!is.null(drugs))
  {
    midDr <- midDr[midDr$drug %in% drugs, ]
  }

  if(nrow(midDr)==0)
  {
    msg1 <- sprintf("drug %s not present in Xeva dataset", paste(drugs, collapse = "\n"))
    stop(msg1)
  }

  bid <- mapModelSlotIds(object, id=midDr$model.id, id.name = "model.id",
                         map.to = "biobase.id", unique = FALSE)
  midDr[, "biobase.id"] <- bid[, "biobase.id"]
  return(midDr)
}


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



##====== drugSensitivitySig for one drug ==========================
#' drugSensitivitySig
#'
#' Get sensitivity signatures for a drug
#'
#' @description
#' Given a Xeva object, and drug name it will return sensitivity value for all the genes/fetures
#'
#' @param object The \code{Xeva} dataset
#' @param mDataType
#' @param drug Name of the drug
#' @param features Which fetures to use from Biobase object. Default \code{NULL} will use all fetures.
#' @return A datafram with fetures and values
#'
#' @examples
#' data(cm.pdxe)
#' drugSensitivitySig(object=cm.pdxe, mDataType="RNASeq", drug="binimetinib", features=1:5,
#' sensitivity.measure="slop")
#'
#' @export
drugSensitivitySig <- function(object, mDataType, drug, features=NULL,
                               sensitivity.measure="slop",
                               #molecular.summary.stat=c("mean", "median", "first", "last", "or", "and"),
                               #sensitivity.summary.stat=c("mean", "median", "first", "last"),
                               #returnValues=c("estimate", "pvalue", "fdr"),
                               #sensitivity.cutoff,
                               standardize=c("SD", "rescale", "none"),
                               nthread=1, verbose=TRUE, ...)
{
  #mDataType = "RNASeq"
  mdf <- .getModDrugBid(object, drug)
  molD <- getMolecularProfiles(object, mDataType)

  mdf <- mdf[ as.character(mdf$biobase.id) %in% colnames(Biobase::exprs(molD)),]

  if(is.null(features))
  { features = rownames(Biobase::exprs(molD))}

  if(verbose==TRUE){printf("Running for drug %s\n\n", drug)}
  mdfI <- .getSensitivityVal(object, sensitivity.measure, mdf, drug=drug, collapse.by="mean")

  rr <- PharmacoGx:::rankGeneDrugSensitivity(data= t(Biobase::exprs(molD)[features, mdfI$biobase.id]),
                                             drugpheno= mdfI[,sensitivity.measure],
                                             #type=type, #batch=batch,
                                             single.type=FALSE,
                                             standardize=standardize,
                                             nthread=nthread,
                                             verbose=verbose)
  return(rr[[1]])
}




