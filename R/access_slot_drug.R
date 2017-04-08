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




##' modelDrugInfo
##' Returns model and drug name
##' @examples
##' data(pdxe)
##' md <- modelDrugInfo(pdxe)
##' md <- modelDrugInfo(pdxe, model.id=c("X.007.biib", "X.6047.LJ16"))
##' md <- modelDrugInfo(pdxe, drug=c("binimetinib", "encorafenib"))
##' @param object The \code{XevaSet} to replace drug info in
##' @param model.id a vector with model.id
##' @param drug a vector with drug names
##' @return Updated \code{XevaSet}
##' #' @export
# modelDrugInfo <- function(object, model.id=NULL, drug=NULL)
# {
#   modDrg <- sapply(slot(object,name="experiment"), "[[", c("drug", "join.name"))
#   modDrgDF <- data.frame(model.id= names(modDrg),
#                          drug = modDrg,
#                          stringsAsFactors = FALSE)
#
#   rownames(modDrgDF) <- modDrgDF$model.id
#
#   if(!is.null(model.id))
#   {
#
#     modDrgDF <- modDrgDF[model.id, ]
#     modDrgDF <- modDrgDF[!is.na(modDrgDF$model.id),]
#   }
#
#   if(!is.null(drug))
#   {
#     modDrgDF <- modDrgDF[modDrgDF$drug %in% c(drug), ]
#   }
#
#   #if(!is.null(other.info))
#   #{
#   #  bid <- mapModelSlotIds(object, id=modDrgDF$model.id, id.name = "model.id",
#   #                         map.to = other.info, unique = FALSE)
#   #  modDrgDF[,other.info] <- bid
#   #}
#
#   return(modDrgDF)
# }







