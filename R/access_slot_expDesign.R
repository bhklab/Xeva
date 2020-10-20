addExtraBatchInfo <- function(moinf, btRTX, infoCol)
{
  #moinf <- modelInfo(object)
  if(!infoCol %in% colnames(moinf))
  { stop(sprintf("%s is not present in modelInfo", infoCol)) }
  for(b in names(btRTX))
  {
    mids <- c()
    if( infoCol == "drug")
    { mids <- btRTX[[b]]$treatment } else 
    { mids <- c(btRTX[[b]]$control, btRTX[[b]]$treatment) }
    
    idx <- c() 
    if( length(mids) > 0 ) { idx <- moinf[mids, infoCol] }
    
    if(length(idx)>0)
    {
      btRTX[[b]][[infoCol]] = names(sort(table(idx), decreasing = T))
    }
  }
  return(btRTX)
}

.convertBatchToDf <- function(btRTX)
{
  allRows <- vapply(btRTX, function(i)i$batch.name, FUN.VALUE = character(1))
  bndf=data.frame(batch.name=allRows, stringsAsFactors = FALSE)
  rownames(bndf) <- NULL
  
  allCols <- unique(unlist(lapply(btRTX, names)))
  col2take <- setdiff(allCols, c("batch.name"))
  bndf[, col2take] <- NA
  
  for(i in 1:nrow(bndf))
  {
    b <- btRTX[[bndf$batch.name[i]]]
    #if(!is.null(b$treatment)){bndf$treatment[i] = paste0(b$treatment, collapse = ";")}
    #if(!is.null(b$control)){bndf$control[i] = paste0(b$control, collapse = ";")}
    for(ci in col2take)
    {
      v <- b[[ci]]
      if(length(v)>0)
      { bndf[i, ci] <- paste0(v, collapse = ";") }
    }
  }
  
  
  # bndf=data.frame(batch.name=vapply(btRTX, function(i)i$batch.name, FUN.VALUE = character(1)),
  #                 treatment=NA, control=NA, patient.id=NA, drug=NA, stringsAsFactors = FALSE)
  # for(i in 1:nrow(bndf))
  # {
  #   b = btRTX[[bndf$batch.name[i]]]
  #   if(!is.null(b$treatment)){bndf$treatment[i] = paste0(b$treatment, collapse = ";")}
  #   if(!is.null(b$control)){bndf$control[i] = paste0(b$control, collapse = ";")}
  #   if(!is.null(b$patient.id)){bndf$patient.id[i] = paste0(b$patient.id, collapse = ";")}
  #   if(!is.null(b$drug)){bndf$drug[i] = paste0(b$drug, collapse = ";")}
  # }
  return(bndf)
}

#' Get batch information
#'
#' Get batch information from a Xeva dataset.
#'
#' @param object The Xeva object from which batch information is obtained.
#' @param batch Name of the batch. Default \code{NULL}.
#' @param model.id Model ID for which need to be searched in the batches. Default \code{NULL}.
#' @param model.id.type Type of the model ID in a batch. See the Details section below.
#'
#' @param patient.id If TRUE, will return patient.id . Default \code{FALSE}.
#' @param drug If TRUE, will return drug . Default \code{FALSE}.
#' @param return.df If TRUE, will return a data.frame (see the Details). Default \code{FALSE}.
#'
#' @details By default this function will return the names of all the batches present in the
#' dataset. If a batch specified, it will return the experiment design (control
#' and treatment model IDs) of that particular batch. If \code{model.id} is specified,
#' it will return the names of all the batches where this particuler \code{model.id} is present.
#' If both \code{batch} and \code{model.id} are specified, \code{batch} will take precedent.
#'
#' For \code{model.id.type}, the default value \code{'any'} will return all the batch IDs
#' where the given model ID is present in any arm (ie. control or treatment) of the
#' batch. It can also be restricted to look only for treatment (or control) arm by
#' specifying the type.
#'
#' If \code{return.df=TRUE}, it will return a data.frame where each batch will be a row.
#' Multiple model.ids will be merged by ; . It will also return patient.id and drug.
#'
#' @examples
#' data(brca)
#' ##to get all the batch names
#' batch.name <- batchInfo(brca)
#'
#' ##to get a specific batch
#' batch.design <- batchInfo(brca, batch=batch.name[1])
#'
#' ##to get all the batches where a model.id is present in control arm
#' bn <- batchInfo(brca, model.id="X.6047.uned", model.id.type="control")
#' 
#' ##get batch level information with tissue column added
#' batch.name <- batchInfo(brca)
#' batch.df <- batchInfo(brca, batch=batch.name[1:10], other.col="tissue",return.df = TRUE)
#'
#' @return A \code{Vector} with batch names.
#' @export
batchInfo <- function(object,
                      batch = NULL,
                      model.id = NULL,
                      model.id.type = c("any", "control", "treatment"),
                      other.col = NULL,
                      return.df = FALSE)
{
  if (is.null(batch) & is.null(model.id))
  {
    rtx <- names(slot(object, "expDesign"))
    return(rtx)
  }
  
  if (is.null(batch) & !is.null(model.id))
  {
    model.id.type <- match.arg(model.id.type)
    rtx <- list()
    for (ed in slot(object, "expDesign"))
    {
      if (model.id.type == "any")
      {
        if (is.element(model.id, ed$treatment) |
            is.element(model.id, ed$control))
        {
          rtx <- .appendToList(rtx, ed$batch.name)
        }
      }
      
      if (model.id.type == "control" &
          is.element(model.id, ed$control))
      {
        rtx <- .appendToList(rtx, ed$batch.name)
      }
      
      if (model.id.type == "treatment" &
          is.element(model.id, ed$treatment))
      {
        rtx <- .appendToList(rtx, ed$batch.name)
      }
    }
    
    return(unique(unlist(rtx)))
  }
  
  if (!is.null(batch) & is.null(model.id))
  {
    btRTX <- list()
    for (bn in c(batch))
    {
      bt <- slot(object, "expDesign")[[bn]]
      if (is.null(bt))
      {
        msg <-
          sprintf("batch name %s not present\nuse batchInfo(object) to see all batch names",
                  bn)
        stop(msg)
      }
      btRTX[[bn]] <- bt
    }
    
    
    extraInfoCol <- c("drug", "patient.id")
    if (!is.null(other.col))
    {
      extraInfoCol <- unique(c(extraInfoCol, other.col))
    }
    
    moinf <- modelInfo(object)
    wrongCol <-
      extraInfoCol[!extraInfoCol %in% colnames(moinf)]
    if (length(wrongCol) != 0)
    {
      stop(sprintf(
        "These columns not present in modelInfo:\n%s\n",
        paste0(wrongCol, collapse = "\n")
      ))
    }
    
    for (infoCol in extraInfoCol)
    {
      btRTX <- addExtraBatchInfo(moinf, btRTX, infoCol)
    }
    
    # if(patient.id==TRUE)
    # {
    #   btRTX <- addExtraBatchInfo(object, btRTX, infoCol="patient.id")
    # }
    #
    # if(drug==TRUE)
    # {
    #   btRTX <- addExtraBatchInfo(object, btRTX, infoCol="drug")
    # }
    
    
    if (return.df == TRUE)
    {
      btRTX <- .convertBatchToDf(btRTX)
    }
    return(btRTX)
  }
}
