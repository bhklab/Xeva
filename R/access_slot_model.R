<<<<<<< HEAD
.addBioBaseIDCol <- function(object, bid, modI)
{
  orgColNames <- colnames(modI)
  mapX <- model2BiobaseIdMap(object, bid)
  rtx <- merge(modI, mapX, by.x=c("model.id"), by.y=c("model.id"))
  newCol <- setdiff(colnames(rtx), colnames(modI))
  modI[, newCol] <- NA
  modNotPres <- setdiff(modI$model.id, rtx$model.id)
  rtx <- rbind(rtx, modI[modNotPres, ])

  colnames(rtx) <- sapply(colnames(rtx), function(x)
                          {
                            if(is.element(x, orgColNames)==FALSE)
                            { x <- sprintf("%s.%s", x,bid)}
                            return(x)
                          })
  return(rtx)
}


=======
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
##----- get modelInfo -------------
#' modelInfo Generic
#' Generic for modelInfo method
#'
#' @examples
#' data(pdxe)
<<<<<<< HEAD
#' mid <- modelInfo(pdxe)
#' head(mid)
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{data.frame} with the model annotations
setGeneric(name = "modelInfo", def = function(object, mDataType=NULL) {standardGeneric("modelInfo")} )

#' @export
setMethod( f=modelInfo, signature="XevaSet",
           definition=function(object, mDataType=NULL)
=======
#' modelInfo(pdxe)
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{data.frame} with the model annotations
setGeneric(name = "modelInfo", def = function(object) {standardGeneric("modelInfo")} )

#' @export
setMethod( f=modelInfo, signature="XevaSet",
           definition=function(object)
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
           {
             modI <- slot(object,name="model")
             drgMod <- sapply(slot(object,name="experiment"), "[[", c("drug", "join.name"))
             modI$drug <- drgMod[ modI$model.id ]
<<<<<<< HEAD

             if(!is.null(mDataType))
             {
               for(bid in c(mDataType))
               {
                 modI <- .addBioBaseIDCol(object, bid, modI)
               }
             }
=======
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
             return(modI)
           } )


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
             warning("This will not update drug information in experiment slot\n write the code to do so...")
             return(object)
           } )

##-----------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------

.checkIfColPresentinModel <- function(object, nameCol)
{
  if(is.element(nameCol, colnames(modelInfo(object))) ==FALSE)
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
<<<<<<< HEAD
#' For example map a model.id to patient.id
#'
#' @examples
#' data(pdxe)
#' mapModelSlotIds(object=pdxe, id="X-007", id.name="patient.id", map.to="model.id")
=======
#' For example map a model.id to biobase.id
#'
#' @examples
#' data(pdxe)
#' mapModelSlotIds(object=pdxe, id="X-007", id.name="biobase.id", map.to="model.id")
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
#' ##map batch ids
#' mapModelSlotIds(pdxe, id= "X-011.INC280", id.name = "batch.name", map.to = "tumor.type")
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

             ##------------------------------------
             if(id.name=="batch.name")
             {
               rtd = .mapBatchName2Id(object, id, map.to)
             } else{

               .checkIfColPresentinModel(object, id.name)
<<<<<<< HEAD
               rtd <- object@model[object@model[,id.name] %in% id, ]
=======
               rtd = object@model[object@model[,id.name] %in% id, ]
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
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




