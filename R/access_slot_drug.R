##----- get drugInfo -------------
#' drugInfo Generic
#' Generic for drugInfo method
#'
#' @examples
#' data(pdxe)
#' drugInfo(pdxe)
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{data.frame} with the drug annotations
setGeneric(name = "drugInfo", def = function(object) {standardGeneric("drugInfo")} )

#' @export
setMethod( f=drugInfo, signature="XevaSet",
           definition=function(object)
             { object@drug } )


#' drugInfo<- Generic
#' Generic for drugInfo replace method
#' @examples
#' data(pdxe)
#' drugInfo(pdxe) <- drugInfo(pdxe)
#' @param object The \code{XevaSet} to replace drug info in
#' @param value A \code{data.frame} with the new drug annotations
#' @return Updated \code{XevaSet}
setGeneric(name= "drugInfo<-", def = function(object, value) {standardGeneric("drugInfo<-")} )

###### @describeIn PharmacoSet Update the drug annotations
#' @export
setMethod( f="drugInfo<-",
           signature=c(object = "XevaSet", value="data.frame"),
           definition=function(object, value)
           {
             object@drug <- value
             return(object)
           } )






