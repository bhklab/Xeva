
.batchID2biobaseID <- function(object, bid, modIn)
{
  bi <- batchInfo(object, batch = bid)
  allMid <- c(bi[[1]]$treatment)
  allMid <- allMid[allMid %in% as.character(modIn$model.id)]

  if(length(allMid)>0)
  {
    rt = modIn[allMid, ]
    rt$batch.name = bid
    return(list(TRUE, rt))
  }
  return(list(FALSE, FALSE))
}



##--------------------------
#' data(brca)
#' object=brca; mDataType="RNASeq"; drug="BGJ398"; tissue=NULL; unique.model=TRUE
#' senType="batch"; sensitivity.measure="abc"
#'
#' object=readRDS("~/CXP/XG/TNBC_analysis/data/TNBC_XevaObj_04Sept20.rds")
#' mDataType="RNASeq"; drug="TAXOL"; tissue=NULL; unique.model=TRUE
#' senType="batch"; sensitivity.measure="abc"

getBatchAndMolData <- function(object, mDataType, drug, tissue, unique.model,
                               senType, sensitivity.measure)
{
  modInX <- .modelID2biobaseID(object, mDataType, drug=drug, tissue=tissue,
                               unique.model=unique.model)
  modIn <- modInX$data; bioName <- modInX$bioName

  sm <- sensitivity(object, type = senType, sensitivity.measure = sensitivity.measure)
  if(is.null(sensitivity.measure))
  { sensitivity.measure <- colnames(sm)[!(colnames(sm) %in% c("model.id", "batch.name"))] }

  bsmat <- data.frame()
  for(i in 1:nrow(sm))
  {
    rt <- .batchID2biobaseID(object, sm$batch.name[i], modIn)
    if(rt[[1]]==TRUE)
    {
      for(q in 1:nrow(rt[[2]]))
      { rt[[2]][q, c(sensitivity.measure)] <- sm[i, c(sensitivity.measure)] }
      bsmat <- rbind(bsmat, rt[[2]])
    }
  }

  if(nrow(bsmat)==0)
  {
    msg <- sprintf("Drug %s not present in treatment arm\n", drug)
    stop(msg)
  }
  biobID2sen <- unique(bsmat[,c(bioName, "batch.name")])
  otherCol <- setdiff(colnames(bsmat), colnames(biobID2sen))
  biobID2sen[, otherCol] <- NA

  for(i in 1:nrow(biobID2sen))
  {
    v=bsmat[bsmat[,bioName]== biobID2sen[i, bioName] &
            bsmat[, "batch.name"]== biobID2sen[i, "batch.name"], ]

    for(ci in otherCol)
    {
      z <- unique(v[,ci])
      if(length(z)>1) { z <- paste0(c(z), collapse = ";")}
      biobID2sen[i, ci] <- z
    }
  }

  rownames(biobID2sen) <- NULL

  molP <- getMolecularProfiles(object, mDataType)
  molP <- molP[, biobID2sen[, bioName]]
  pd <- Biobase::pData(molP)
  for(i in colnames(biobID2sen))
  { pd[,i] <- biobID2sen[,i] }
  Biobase::pData(molP) <- pd
  return(molP)
}



