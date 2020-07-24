#' Get drug information
#' Get the drug information slot from a {XevaSet} object.
#'
#' @examples
#' data(brca)
#' head(drugInform(brca))
#' @param object The \code{XevaSet} to retrieve drug information from.
#' @return A \code{data.frame} with the drug annotations.
setGeneric(name = "drugInform", def = function(object) {standardGeneric("drugInform")} )

#' @rdname drugInform
#' @export
setMethod( f=drugInform, signature="XevaSet",
           definition=function(object)
             { slot(object, "drug") } )
