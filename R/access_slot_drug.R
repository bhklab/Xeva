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

#' Set new drug information in a XevaSet
#' A method to set the drug information slot in a \code{XevaSet} object.
#'
#' @examples
#' data(brca)
#' drugInfo(brca)<- drugInfo(brca)
#'
#' @param object The \code{XevaSet} to replace drug info in.
#' @param value A \code{data.frame} with the new drug annotations.
#' @return updated \code{XevaSet}
#'
#' @docType methods
#' @rdname SetDrugInfo-methods
setGeneric(name= "drugInfo<-", def = function(object, value) {standardGeneric("drugInfo<-")} )

#' @rdname SetDrugInfo-methods
#' @export
setMethod( f="drugInfo<-",
           signature=c(object = "XevaSet", value="data.frame"),
           definition=function(object, value)
           {
             slot(object, "drug") <- value
             return(object)
           } )

