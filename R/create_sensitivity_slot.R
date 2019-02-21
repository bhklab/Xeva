##---- sensitivity ----------
#@model
#@batch

.checkUnqLength <- function(inVec)
{ length(inVec)== length(unique(inVec)) }

.creatSensitivitySlot <- function(modelSensitivity, batchSensitivity, expSlot, expDesign)
{
  ##------------ for modelSensitivity ------------------------------------------------
  if(nrow(modelSensitivity)==0)
  {
    modelSensitivity <-data.frame(model.id= names(expSlot), stringsAsFactors = FALSE)
  }

  for(mid in names(expSlot))
  {
    if(is.element(mid, modelSensitivity$model.id) ==FALSE)
    {
      msg <- sprintf("provide modelSensitivity for all models\nmodelSensitivity missing for %s", mid)
      stop(msg)
    }
  }

  if( .checkUnqLength(modelSensitivity$model.id)==FALSE)
  {stop("model.ids are not unique")}
  rownames(modelSensitivity) <- as.character(modelSensitivity$model.id)
  modelSensitivity <- modelSensitivity[names(expSlot), ,drop=FALSE]
  ##-----------------------------------------------------------------------------
  ##------------ for Batch Sensitivity -------------------------------------------

  if(nrow(batchSensitivity)>0)
  {
    if(is.element("batch.name", colnames(batchSensitivity)) ==FALSE)
    {
      stop("in 'batchSensitivity' datafram one column must be 'batch.name'")
    }

    if(nrow(batchSensitivity)!= length(names(expDesign)))
    {
      batchSensitivity <- .reorderCol(batchSensitivity, "batch.name", 1)
      missingId <- unique(setdiff(names(expDesign), batchSensitivity$batch.name))
      bsN <- data.frame(matrix(NA, nrow = length(missingId), ncol = ncol(batchSensitivity)))
      colnames(bsN) <- colnames(batchSensitivity)
      bsN$batch.name <- missingId
      batchSensitivity <- rbind(batchSensitivity, bsN)
    }
  }

  if(nrow(batchSensitivity)==0)
  {
    batchSensitivity <- data.frame(batch.name= names(expDesign), stringsAsFactors = FALSE)
  }

  if( .checkUnqLength(batchSensitivity$batch.name)==FALSE)
  {stop("batch names are not unique")}
  rownames(batchSensitivity) <- as.character(batchSensitivity$batch.name)
  batchSensitivity <- batchSensitivity[names(expDesign), ,drop=FALSE]

  ##--------------------------------------------------------------------------------
  if(is(modelSensitivity, "data.frame")==FALSE | is(batchSensitivity, "data.frame")==FALSE)
  {stop("slot class error")}

  rtx <- list(model = modelSensitivity,
              batch = batchSensitivity)
  return(rtx)
}

