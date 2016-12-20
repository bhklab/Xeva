##----- get modelInfo -------------
#' modelInfo Generic
#' Generic for ModelInfo method
#'
#' @examples
#' data(pdxe)
#' ModelInfo(pdxe)
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{data.frame} with the model annotations
setGeneric(name = "ModelInfo", def = function(object) {standardGeneric("ModelInfo")} )

#' @export
setMethod( f=ModelInfo, signature="XevaSet",
           definition=function(object)
           { object@model } )


#' ModelInfo<- Generic
#' Generic for ModelInfo replace method
#' @examples
#' data(pdxe)
#' ModelInfo(pdxe) <- ModelInfo(pdxe)
#' @param object The \code{XevaSet} to replace drug info in
#' @param value A \code{data.frame} with the new model annotations
#' @return Updated \code{XevaSet}
setGeneric(name= "ModelInfo<-", def = function(object, value) {standardGeneric("ModelInfo<-")} )

#' @export
setMethod( f="ModelInfo<-",
           signature=c(object = "XevaSet", value="data.frame"),
           definition=function(object, value)
           {
             object@model = value
             return(object)
           } )

##-----------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------

.checkIfColPresentinModel <- function(object, nameCol)
{
  if(is.element(nameCol, colnames(object@model)) ==FALSE)
  {
    msg = sprintf("%s is not a valid id.name\nValid id.names are:\n\n%s",
                  nameCol, paste(colnames(object@model), collapse = "\n"))
    stop(msg)
  }
}

##---------------------------------------------------
##---------------------------------------------------
#' map one id type to another in model slot
#' model.id given a "biobase.id" or "patient id" or any other id that is a column in model slot
#'
#' @examples
#' data(pdxe)
#' mapModelSlotIds(object=pdxe, id="X-007", id.name="biobase.id", map.to="model.id")
#' @param object The \code{Xeva} dataset
#' @param id The \code{id}
#' @param id.name The \code{id} name
#' @param map.to The name of the mapped id
#' @return a \code{data.fram} with id and mapped id
setGeneric(name = "mapModelSlotIds", def = function(object, id, id.name, map.to) {standardGeneric("mapModelSlotIds")})

#' @export
setMethod( f=mapModelSlotIds,
           signature=c(object="XevaSet"),
           definition= function(object, id, id.name, map.to="all")
           {
             id = c(id)
             .checkIfColPresentinModel(object, id.name)

             rtd = object@model[object@model[,id.name] %in% id, ]
             if(map.to!="all")
             {
               .checkIfColPresentinModel(object, map.to)
               rtd = rtd[, c(id.name, map.to)]
               if(id.name==map.to){rtd = rtd[,id.name, drop=FALSE]}
              }
             return(rtd)
           })






