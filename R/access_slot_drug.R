#' Get drug information
#' Get the drug information slot from a {XevaSet} object.
#'
#' @examples
#' data(brca)
#' head(drugInfo(brca))
#' @param object The \code{XevaSet} to retrieve drug information from.
#' @return A \code{data.frame} with the drug annotations.
setGeneric(name = "drugInfo", def = function(object) {standardGeneric("drugInfo")} )

#' @rdname drugInfo
#' @export
setMethod( f=drugInfo, signature="XevaSet",
           definition=function(object)
             { slot(object, "drug") } )
