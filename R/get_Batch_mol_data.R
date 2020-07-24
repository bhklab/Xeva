
.batchID2biobaseID <- function(object, bid, modIn)
{
  bi = batchInfo(object, batch = bid)
  allMid = c(bi[[1]]$treatment, bi[[1]]$control)
  allMid = allMid[allMid %in% as.character(modIn$model.id)]

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

getBatchAndMolData <- function(object, mDataType, drug, tissue, unique.model,
                               senType, sensitivity.measure)
{
  modInX <- .modelID2biobaseID(object, mDataType, drug=drug, tissue=tissue,
                               unique.model=unique.model)
  modIn <- modInX$data; bioName <- modInX$bioName

  sm <- sensitivity(object, type = senType, sensitivity.measure = sensitivity.measure)
  if(is.null(sensitivity.measure))
  { sensitivity.measure <- colnames(sm)[!(colnames(sm) %in% c("model.id", "batch.name"))] }

  bsmat=data.frame()
  for(i in 1:nrow(sm))
  {
    bid = sm$batch.name[i]
    rt = .batchID2biobaseID(object, bid, modIn)
    if(rt[[1]]==TRUE)
    {
      for(q in 1:nrow(rt[[2]]))
      { rt[[2]][q, c(sensitivity.measure)] <- sm[i, c(sensitivity.measure)] }
      bsmat <- rbind(bsmat, rt[[2]])
    }
  }

  biobID2sen <- unique(bsmat[,c(bioName, "batch.name")])
  rmCol= setdiff(colnames(bsmat), colnames(biobID2sen))
  biobID2sen[, rmCol] <- NA
  #biobID2sen <- data.frame(matrix(NA, nrow = length(unique(bsmat[,bioName])),
  #                                ncol = ncol(bsmat)), stringsAsFactors = F)
  #colnames(biobID2sen) <- colnames(bsmat)
  #biobID2sen[,bioName] <- unique(bsmat[,bioName])
  otherCol <- setdiff(colnames(biobID2sen), c(bioName, "batch.name"))

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

  rownames(biobID2sen) <- biobID2sen[, bioName]

  molP <- getMolecularProfiles(object, mDataType)
  molP <- molP[, biobID2sen[, bioName]]
  pd <- Biobase::pData(molP)
  for(i in colnames(biobID2sen))
  { pd[,i] <- biobID2sen[,i] }
  Biobase::pData(molP) <- pd
  return(molP)
}



