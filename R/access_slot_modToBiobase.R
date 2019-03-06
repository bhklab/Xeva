#' model2BiobaseIdMap
#' Gives model.id to biobase.id mapping datafram
#'
#' @examples
#' data(brca)
#' idMap <- model2BiobaseIdMap(brca, mDataType="RNASeq")
#' head(idMap)
#' @param object The \code{XevaSet}.
#' @param mDataType Data type for which ids to be retrive. Default \code{NULL} will return a full data frame.
#' @return a \code{data.frame} with the model.id and biobase.id
#' @keywords internal
#' @noRd
##### @export
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
  rownames(idMap) <- NULL
  return(idMap)
}
