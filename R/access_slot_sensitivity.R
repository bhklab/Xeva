## get sensitivity for an Xeva object
#'
#' Get sensitivity for an Xeva object
#' @description
#' Given a Xeva object, it will return sensitivity datafram
#'
#' @examples
#' data(cm.pdxe)
#' head(sensitivity(cm.pdxe, type="batch"))
#' head(sensitivity(cm.pdxe, type="model"))
#' @param object The \code{Xeva} dataset
#' @param type sensitivity type (either model or batch)
#'
#' @return a \code{data.fram} with model or batch id and sensitivity values
#' @export
sensitivity <- function(object, type)
{
  if(is.element(type, names(object@sensitivity))==TRUE )
  {
    return(object@sensitivity[[type]])
  } else
  {
    msg <- sprintf("sensitivity 'type' can be:\n%s", paste(names(object@sensitivity), collapse = "\n"))
    stop(msg)
  }
}
