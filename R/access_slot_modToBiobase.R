# .modelId2BiobaseId_old <- function(object, mid, mDataType)
# {
#   idMap <- slot(object, "modToBiobaseMap")
#   idMap <- idMap[idMap$mDataType == mDataType, ]
#
#   if(nrow(idMap)>0)
#   {
#     rtx <- list() #data.frame()
#     unqMid <- unique(idMap$model.id)
#     for(mx in c(mid))
#     {
#       if(is.element(mx, unqMid))
#       {
#         bid <- idMap[idMap$model.id==mx, "biobase.id"]
#       } else
#       { bid <- NA }
#       rtx[[length(rtx)+1]] <- setNames(c(mx, bid), c("model.id", "biobase.id"))
#     }
#   } else
#   {
#     msg <- sprintf("mDataType %s not present in moleculer profile\n", mDataType)
#     stop(msg)
#   }
#
#   rtz <- data.frame(model.id = sapply(rtx, `[[`, "model.id"),
#                     biobase.id=sapply(rtx, `[[`, "biobase.id"),
#                     stringsAsFactors = FALSE)
#   return(rtz)
# }


#' model2BiobaseIdMap
#' Gives model.id to biobase.id mapping datafram
#'
#' @examples
#' data(pdxe)
#' idMap <- model2BiobaseIdMap(pdxe, mDataType="RNASeq")
#' head(idMap)
#' @param object The \code{XevaSet}
#' @param mDataType Data type for which ids to be retrive. Default \code{NULL} will return full datafram
#' @return a \code{data.frame} with the model.id and biobase.id
#' @export
model2BiobaseIdMap <- function(object, mDataType=NULL)
{
  idMap <- slot(object, "modToBiobaseMap")
  if(!is.null(mDataType))
  {
    idMap <- idMap[idMap$mDataType == mDataType, ]
    if(nrow(idMap)==0)
    {
      msg <- sprintf("mDataType %s not present in moleculer profile\n", mDataType)
      stop(msg)
    }
  }

  # if(!is.null(model.id))
  # {
  #   model.id <- c(model.id)
  #   idx <- unlist(sapply(model.id, function(x) {which(idMap$model.id==x)}))
  #   idMap <- idMap[idx, ]
  #   npx <- setdiff(model.id, idMap$model.id)
  #   if(length(npx)>0)
  #   {
  #     for(x in npx)
  #     {idMap[nrow(idMap)+1, "model.id" ] <- x}
  #   }
  # }
  rownames(idMap) <- NULL
  return(idMap)
}





#
# .modelId2BiobaseId <- function(object, model.id=NULL, mDataType=NULL)
# {
#   idMap <- slot(object, "modToBiobaseMap")
#   if(!is.null(mDataType))
#   {
#     idMap <- idMap[idMap$mDataType == mDataType, ]
#     if(nrow(idMap)==0)
#     {
#       msg <- sprintf("mDataType %s not present in moleculer profile\n", mDataType)
#       stop(msg)
#     }
#   }
#
#   if(!is.null(model.id))
#   {
#     model.id <- c(model.id)
#     idx <- unlist(sapply(model.id, function(x) {which(idMap$model.id==x)}))
#     idMap <- idMap[idx, ]
#     npx <- setdiff(model.id, idMap$model.id)
#     if(length(npx)>0)
#     {
#       for(x in npx)
#       {idMap[nrow(idMap)+1, "model.id" ] <- x}
#     }
#   }
#   rownames(idMap) <- NULL
#   return(idMap)
# }
#


#library(profvis)
#profvis(
#   mx <- .modelId2BiobaseId(object, modI$model.id, mdtyp)
# )
#
#
# #require(profr)
# #require(ggplot2)
# x = profr(.modelId2BiobaseId(object, modI$model.id, mdtyp))
# ggplot(x)

