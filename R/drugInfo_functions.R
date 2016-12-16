
##----- get drugInfo -------------
#' drugInfo Generic
#' Generic for drugInfo method
#'
#' @examples
#' data(pdxe)
#' drugInfo(pdxe)
#' @param object The \code{XenoSet} to retrieve drug info from
#' @return a \code{data.frame} with the drug annotations
setGeneric(name = "drugInfo", def = function(object) {standardGeneric("drugInfo")} )

#### @describeIn PharmacoSet Returns the annotations for all the drugs tested in the PharmacoSet
#' @export
setMethod( f=drugInfo, signature="XenoSet", definition=function(object){ dim(object@drug) } )



#' drugInfo<- Generic
#' Generic for drugInfo replace method
#' @examples
#' data(pdxe)
#' drugInfo(pdxe) <- drugInfo(pdxe)
#' @param object The \code{XenoSet} to replace drug info in
#' @param value A \code{data.frame} with the new drug annotations
#' @return Updated \code{XenoSet}
setGeneric(name= "drugInfo<-", def = function(object, value) {standardGeneric("drugInfo<-")} )

###### @describeIn PharmacoSet Update the drug annotations
#' @export
setMethod( f="drugInfo<-",
           signature="XenoSet",
           definition=function(object, value)
           {
             object@annotation$drugInfo = value
             #object@drugInfo = value  ##This will not work as slot drugInfo already have to be present
             object
           } )

