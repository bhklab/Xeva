#####================= getMolecularProfiles ==================
#' Get Molecular Profiles
#'
#' Get Molecular Profiles
#'
#' @param object The \code{XevaSet}
#' @param data.type \code{character}, which one of the molecular data types is needed
#' @return a \code{ExpressionSet} where sample names are \code{biobase.id} of model
#' @examples
#' data(pdxe)
#' pdxe_RNA <- getMolecularProfiles(pdxe, data.type="RNASeq")
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


.batch2DataFram <- function(object, batchName=NULL, expDig=NULL)
{
  if(is.null(batchName) & is.null(expDig))
  { stop("please provide 'batchName' or 'expDig'") }

  if(!is.null(batchName))
  { expDig <- expDesign(object, batchName) }

  bat2mods <- data.frame()
  for(i in 1:length(expDig))
  {
    expdi <- expDig[[i]]

    if(!is.null(expdi$control))
    {
      bat2mods <- rbind(bat2mods, data.frame(batch.name=expdi$batch.name,
                                             model.id=expdi$control,
                                             exp.type="control",
                                             stringsAsFactors = F))
    }
    if(!is.null(expdi$treatment))
    {
      bat2mods <- rbind(bat2mods, data.frame(batch.name=expdi$batch.name,
                                             model.id=expdi$treatment,
                                             exp.type="treatment",
                                             stringsAsFactors = F))
    }
  }
  return(bat2mods)
}

.modelID2biobaseID <- function(object, mDataType, drug=NULL, tumor.type=NULL, unique.model=TRUE)
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

  if(!is.null(tumor.type))
  {
    modIn <- modIn[modIn$tumor.type %in% c(tumor.type), ]
    if(nrow(modIn)==0)
    {
      msg <- sprintf("No model present with tumor.type %s ", tumor.type)
      stop(msg)
    }
  }

  bioName <- sprintf("biobase.id.%s", mDataType)
  modIn <- modIn[ !is.na(modIn[, bioName]), ]
  if(nrow(modIn)==0)
  {
    msg <- sprintf("No model present for drug %s with molecular data type %s", drug, mDataType)
    if(!is.null(tumor.type))
    {
      msg <- sprintf("No model present for drug %s and tumor type %s with molecular data type %s",
                     drug, tumor.type, mDataType)
    }
    stop(msg)
  }
  #return(modIn)
  return(list(data=modIn, bioName=bioName))
}
#####================= Summarize Molecular Profiles ==================
#' summarizeMolecularProfiles
#'
#' summarizeMolecularProfiles
#'
#' @param object The \code{XevaSet}
#' @param drug Name of the drug
#' @param mDataType \code{character}, which one of the molecular data types is needed
#' @param tumor.type default \code{NULL} will return all across all tumor.type
#' @param sensitivity.measure default \code{NULL} will return all sensitivity measure
#' @param unique.model default TRUE will return only one sequncing id, in case where one model id mapes to several sequencing ids
#' @return A \code{ExpressionSet} where sample names are model.id and sensitivity measure will be present in pData
#' @examples
#' data(pdxe)
#' pacRNA <- summarizeMolecularProfiles(pdxe, drug="paclitaxel", mDataType="RNASeq",
#'                                      tumor.type= "BRCA", sensitivity.measure="mRECIST")
#' print(pacRNA)
#' @details
#' \itemize{
#' \item {If a sequencing sample belong to multipal models, summarizeMolecularProfiles
#' will creat saperate column for each model. }
#' \item {All the models without the moleculer data will be removed from the output expression set.}
#' }
#' @export
summarizeMolecularProfiles <- function(object, drug, mDataType, tumor.type=NULL,
                                       sensitivity.measure=NULL, unique.model=TRUE,
                                       batchName=NULL, expDig=NULL)
{
  senType <- NULL
  if(is.null(sensitivity.measure))
  {
    senType <- "model"
  } else
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
    msg <- sprintf("sensitivity measure '%s' not present in model or batch sensitivity. First use setSensitivity function to set sensitivity", sensitivity.measure)
    stop(msg)
  }
  ##----------------------------------------------------------------------------
  ##----------------------------------------------------------------------------
  #sm <- sensitivity(object, type = "model", sensitivity.measure = sensitivity.measure)

  #sm <- sensitivity(object, type = senType, sensitivity.measure = sensitivity.measure)
  #if(is.null(sensitivity.measure))
  #{ sensitivity.measure <- colnames(sm)[!(colnames(sm) %in% c("model.id", "batch.name"))] }

  if(senType=="model")
  {
    sm <- sensitivity(object, type = senType, sensitivity.measure = sensitivity.measure)
    if(is.null(sensitivity.measure))
    { sensitivity.measure <- colnames(sm)[!(colnames(sm) %in% c("model.id", "batch.name"))] }

    modInX <- .modelID2biobaseID(object, mDataType, drug=drug, tumor.type=tumor.type,
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
    if(is.null(batchName) & is.null(expDig))
    {
      sm <- sensitivity(object, type = senType, sensitivity.measure = sensitivity.measure)
    }

    if(!is.null(batchName))
    {
      sm <- sensitivity(object, type = senType, sensitivity.measure = sensitivity.measure)
      sm <- sm[sm$batch.name %in% batchName, ]
    }

    if(is.null(batchName) & !is.null(expDig))
    {
      #Get time var data and compute angle etc.
      stop("not implemented yet")
    }

    if(is.null(sensitivity.measure))
    { sensitivity.measure <- colnames(sm)[!(colnames(sm) %in% c("model.id", "batch.name"))] }

    btDF <- .batch2DataFram(object, sm$batch.name)
    btDF[, sensitivity.measure] <- NA
    unqBtach <- intersect(sm$batch.name, unique(btDF$batch.name))
    for(bn in unqBtach)
    {
      btDF[btDF$batch.name==bn, sensitivity.measure] <- sm[sm$batch.name==bn, sensitivity.measure]
    }

    modInX <- .modelID2biobaseID(object, mDataType, drug=drug, tumor.type=tumor.type,
                                 unique.model=unique.model)
    modIn <- modInX$data; bioName <- modInX$bioName

    commanMod <- intersect(btDF$model.id, modIn$model.id)
    if(length(commanMod)==0)
    {
      msg <- sprintf("No model present for given batch, drug and tissue type")
      stop(msg)
    }

    bdf <- merge(btDF, modIn, by.x = "model.id", by.y="model.id")
    bdf <- bdf[!is.na(bdf[,bioName]), ]
    rownames(bdf) <- paste0("R", 1:nrow(bdf))

    bt2bid <- data.frame()
    for(bn in unique(bdf$batch.name))
    {
      v <- bdf[bdf$batch.name==bn, ]
      for(I in 1:nrow(v))
      {
        if(nrow(bt2bid)==0)
        { bt2bid <- rbind(bt2bid, v[I,])}
        if(!(v[I, bioName] %in% bt2bid[, bioName]))
        { bt2bid <- rbind(bt2bid, v[I,]) }
      }
    }

    bt2bid[, bioName] <- as.character(bt2bid[, bioName])
    rownames(bt2bid) <- make.names(bt2bid[, bioName], unique = T)

    molP <- getMolecularProfiles(object, mDataType)
    molP <- molP[, bt2bid[, bioName]]
    colnames(molP) <- rownames(bt2bid)
    rownames(Biobase::pData(molP)) <- rownames(bt2bid)
    pd <- Biobase::pData(molP)
    for(i in colnames(bt2bid))
    { pd[,i] <- bt2bid[,i] }
    Biobase::pData(molP) <- pd
    return(molP)
  }

}



# old_summarizeMolecularProfiles <- function(object, drug, mDataType, tumor.type=NULL,
#                                            sensitivity.measure=NULL, unique.model=TRUE)
# {
#   modIn <- modelInfo(object, mDataType = mDataType)
#   if(unique.model==TRUE)
#   {
#     modIn <- modIn[unique(modIn$model.id), ]
#   }
#   modIn <- modIn[modIn$drug %in% c(drug), ]
#   if(nrow(modIn)==0)
#   {
#     msg <- sprintf("No model present with drug %s ", drug)
#     stop(msg)
#   }
#
#   if(!is.null(tumor.type))
#   {
#     modIn <- modIn[modIn$tumor.type %in% c(tumor.type), ]
#     if(nrow(modIn)==0)
#     {
#       msg <- sprintf("No model present with tumor.type %s ", tumor.type)
#       stop(msg)
#     }
#   }
#
#   bioName <- sprintf("biobase.id.%s", mDataType)
#   modIn <- modIn[ !is.na(modIn[, bioName]), ]
#   if(nrow(modIn)==0)
#   {
#     msg <- sprintf("No model present for drug %s with molecular data type %s", drug, mDataType)
#     if(!is.null(tumor.type))
#     {
#       msg <- sprintf("No model present for drug %s and tumor type %s with molecular data type %s",
#                      drug, tumor.type, mDataType)
#     }
#     stop(msg)
#   }
#
#   sm <- sensitivity(object, type = "model", sensitivity.measure = sensitivity.measure)
#
#   if(is.null(sensitivity.measure))
#   {
#     sensitivity.measure <- colnames(sm)[colnames(sm)!="model.id"]
#   }
#
#   modIn[,c(sensitivity.measure)] <- sm[modIn$model.id, c(sensitivity.measure)]
#
#   molP <- getMolecularProfiles(object, mDataType)
#   molP <- molP[, modIn[, bioName]]
#   colnames(molP) <- rownames(modIn)
#   rownames(Biobase::pData(molP)) <- rownames(modIn)
#
#   pd <- Biobase::pData(molP)
#   for(i in colnames(modIn))
#   { pd[,i] <- modIn[,i] }
#   Biobase::pData(molP) <- pd
#   return(molP)
# }


