#' get drug information
#' get drug information slot
#'
#' @examples
#' data(brca)
#' drugInfo(brca)
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{data.frame} with the drug annotations
setGeneric(name = "drugInfo", def = function(object) {standardGeneric("drugInfo")} )
#' @export
setMethod( f=drugInfo, signature="XevaSet",
           definition=function(object)
             { slot(object, "drug") } )

#' set drug information
#' set drug information slot
#' @examples
#' data(brca)
#' drugInfo(brca)<- drugInfo(brca)
#' @param object The \code{XevaSet} to replace drug info in
#' @param value A \code{data.frame} with the new drug annotations
#' @return updated \code{XevaSet}
setGeneric(name= "drugInfo<-", def = function(object, value) {standardGeneric("drugInfo<-")} )
#' @export
setMethod( f="drugInfo<-",
           signature=c(object = "XevaSet", value="data.frame"),
           definition=function(object, value)
           {
             slot(object, "drug") <- value
             return(object)
           } )

