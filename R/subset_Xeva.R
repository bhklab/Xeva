.resetBatchDesign <- function(batch, vid)
{
  rtx <- list()
  for(i in names(batch))
  {
    rtB <- list()
    rtB$batch.name <- batch[[i]]$batch.name

    if(!is.null(batch[[i]]$treatment))
    {
      rtB$treatment <- batch[[i]]$treatment[batch[[i]]$treatment %in% vid]
      if(length(rtB$treatment)==0){ rtB$treatment<- NULL}
    }

    if(!is.null(batch[[i]]$control))
    {
      rtB$control <- batch[[i]]$control[batch[[i]]$control %in% vid]
      if(length(rtB$control)==0){ rtB$control<- NULL}
    }

    rtx[[i]] <- rtB
    if(is.null(rtx[[i]]$treatment) & is.null(rtx[[i]]$control))
    {rtx[[i]] <- NULL}
  }
  return(rtx)
}



#' Subset Xeva object.
#'
#' @examples
#' data(brca)
#' print(brca)
#' df <- subsetXeva(brca, ids = c("X-1004", "X-1008", "X-1286"), id.name="patient.id", keep.batch=TRUE)
#' print(df)
#' @param object The \code{XevaSet} object.
#' @param ids IDs to be selected for.
#' @param id.name Names of the IDs.
#' @param keep.batch Default \code{TRUE}. If \code{FALSE}, remove all other \code{model.ids} from the experiemt design that do not belong to selection.
#' @return New Xeva object.
#' @export
subsetXeva <- function(object, ids, id.name, keep.batch=TRUE)
{
  md <- modelInfo(object)
  if(is.element(id.name, colnames(md))==FALSE)
  {
    msg <- sprintf("'id.name = %s' not present in modelInfo\nValid 'id.name' are\n%s\n",
                   id.name, paste(colnames(md), collapse = "\n"))
    stop(msg)
  }
  ids <- unique(c(ids))
  mdn <- md[md[, id.name] %in% ids, ]
  if(nrow(mdn)==0)
  {
    warning("No model for input ids present, returning NULL")
    return(NULL)
  }

  expDesign <- slot(object, "expDesign")
  expDeNew <- list()
  for(i in seq_along(expDesign))
  {
    bn <- expDesign[[i]]$batch.name
    tr <- expDesign[[i]]$treatment
    ct <- expDesign[[i]]$control
    if(!is.null(tr) & any(tr %in% mdn$model.id))
    {
      expDeNew[[bn]] <- expDesign[[i]]
      next()
    }

    if(!is.null(ct) & any(ct %in% mdn$model.id))
    {
      expDeNew[[bn]] <- expDesign[[i]]
      next()
    }
  }

  if(keep.batch==FALSE)
  { expDeNew <- .resetBatchDesign(expDeNew, mdn$model.id) }

  nwModId <- unique( unlist(lapply(expDeNew,
                                   function(x){c(x[["treatment"]],x[["control"]])})))
  newModId <- unique(c(mdn$model.id, nwModId))
  newModId <- newModId[newModId %in% rownames(md)]
  mdn <- md[newModId, ]

  slot(object, "expDesign") <- expDeNew
  slot(object, "experiment") <- slot(object, "experiment")[mdn$model.id]
  slot(object, "model") <- slot(object, "model")[mdn$model.id, ]
  slot(object, "drug")  <- slot(object, "drug")[slot(object, "drug")$drug.id
                                                %in% mdn$drug,]
  sn <- slot(object, "sensitivity")
  sn$model <- sn$model[sn$model$model.id %in% mdn$model.id,  , drop=F]
  sn$batch <- sn$batch[sn$batch$batch.name %in% names(expDeNew), , drop=F]
  slot(object, "sensitivity") <- sn

  m2b <- slot(object, "modToBiobaseMap")
  slot(object, "modToBiobaseMap") <- m2b[m2b$model.id %in% mdn$model.id, ]
  m2b <- slot(object, "modToBiobaseMap")
  for(mold in names(slot(object, "molecularProfiles")) )
  {
    ids2take <- unique(m2b[m2b$mDataType== mold, "biobase.id"])
    ids2take <- ids2take[!is.na(ids2take)]
    if(length(ids2take)>0)
    {
      mol <- slot(object, "molecularProfiles")[[mold]]
      slot(object, "molecularProfiles")[[mold]] <- mol[, ids2take]
    }
  }

  nonNAExp <- vapply(slot(object, name = "experiment"), function(mod){!is.null(mod)},
                     FUN.VALUE = logical(1))
  slot(object, "experiment") <- slot(object, "experiment")[nonNAExp]

  return(object)
}
