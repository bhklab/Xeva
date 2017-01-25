##----- get modelInfo -------------
#' modelInfo Generic
#' Generic for modelInfo method
#'
#' @examples
#' data(pdxe)
#' modelInfo(pdxe)
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{data.frame} with the model annotations
setGeneric(name = "modelInfo", def = function(object) {standardGeneric("modelInfo")} )

#' @export
setMethod( f=modelInfo, signature="XevaSet",
           definition=function(object)
           { slot(object,name="model") } )


#' modelInfo<- Generic
#' Generic for modelInfo replace method
#' @examples
#' data(pdxe)
#' modelInfo(pdxe) <- modelInfo(pdxe)
#' @param object The \code{XevaSet} to replace drug info in
#' @param value A \code{data.frame} with the new model annotations
#' @return Updated \code{XevaSet}
setGeneric(name= "modelInfo<-", def = function(object, value) {standardGeneric("modelInfo<-")} )

#' @export
setMethod( f="modelInfo<-",
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
##
#' Map ids of model slot
#'
#'
#' Map one id type to another in model slot.
#' For example map a model.id to biobase.id
#'
#' @examples
#' data(pdxe)
#' mapModelSlotIds(object=pdxe, id="X-007", id.name="biobase.id", map.to="model.id")
#' @param object The \code{Xeva} dataset
#' @param id The \code{id}
#' @param id.name The \code{id} name
#' @param map.to The name of the mapped id. Default \code{all}
#' @param unique Default \code{TRUE}. If unique=FALSE output will be mapped to input
#' @return a \code{data.fram} with id and mapped id
setGeneric(name = "mapModelSlotIds", def = function(object, id, id.name, map.to="all",unique=TRUE) {standardGeneric("mapModelSlotIds")})

#' @export
setMethod( f=mapModelSlotIds,
           signature=c(object="XevaSet"),
           definition= function(object, id, id.name, map.to="all", unique=TRUE)
           {
             id = c(as.character(id))

             if(id.name=="batch.name")
             {
               rtd = .mapBatchName2Id(object, id, map.to)
             } else{

               .checkIfColPresentinModel(object, id.name)
               rtd = object@model[object@model[,id.name] %in% id, ]
               if(map.to!="all")
               {
                 .checkIfColPresentinModel(object, map.to)
                 rtd = rtd[, c(id.name, map.to)]
                 if(id.name==map.to){rtd = rtd[,id.name, drop=FALSE]}
                 rtd = unique(rtd)
               }

               if(unique==FALSE)
               { rtd = rtd[match(id, rtd[,id.name]),] }
             }

             return(rtd)
           })

##--------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
##-----map batch to patient.id --------------------------------------------------------
##
.mapBatchName2Id <- function(object, id, map.to)
{
  btMapRet = data.frame()
  for(batch.name in id)
  {
    bt <- expDesign(object, batch.name = batch.name)
    bt.Mod <- unique(c(bt$treatment, bt$control))
    btMap = mapModelSlotIds(object, id= bt.Mod, id.name="model.id", map.to=map.to, unique=TRUE)
    btMap[, "batch.name"] = batch.name
    btMapRet = rbind(btMapRet, btMap)
  }
  btMapRet <- btMapRet[, c("batch.name", map.to)]
  btMapRet <- unique(btMapRet); rownames(btMapRet)=NULL
  return(btMapRet)
}




