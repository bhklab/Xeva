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





