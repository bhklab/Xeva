.addBioBaseIDCol <- function(object, bid, modI)
{
  orgColNames <- colnames(modI)
  mapX <- model2BiobaseIdMap(object, bid)
  rtx <- merge(modI, mapX, by.x=c("model.id"), by.y=c("model.id"))
  newCol <- setdiff(colnames(rtx), colnames(modI))
  modI[, newCol] <- NA
  modNotPres <- setdiff(modI$model.id, rtx$model.id)
  rtx <- rbind(rtx, modI[modNotPres, ])

  colnames(rtx) <- vapply(colnames(rtx), function(x)
                          {
                            if(is.element(x, orgColNames)==FALSE)
                            { x <- sprintf("%s.%s", x,bid)}
                            return(x)
                          }, FUN.VALUE = character(1))
  return(rtx)
}


##----- get modelInfo -------------
#' modelInfo Generic
#' Generic for modelInfo method
#'
#' @param object Xeva object
#' @param mDataType Molecular data type.
#'
#' @examples
#' data(brca)
#' mid <- modelInfo(brca)
#' head(mid)
#'
#' @return A \code{data.frame} with the model annotations.
setGeneric(name = "modelInfo", def = function(object, mDataType=NULL) {standardGeneric("modelInfo")} )

#' @rdname modelInfo
#' @export
setMethod( f=modelInfo, signature="XevaSet",
           definition=function(object, mDataType=NULL)
           {
             modI <- slot(object,name="model")
             drgMod <- vapply(slot(object,name="experiment"), function(mod)
                              { slot(mod,name="drug")[["join.name"]] },
                              FUN.VALUE = character(1))

             modI$drug <- drgMod[ modI$model.id ]

             if(!is.null(mDataType))
             {
               for(bid in c(mDataType))
               {
                 modI <- .addBioBaseIDCol(object, bid, modI)
               }
             }
             rownames(modI) <- make.names(modI$model.id, unique = TRUE)
             #rownames(modI) <- as.character(modI$model.id)
             return(modI)
           } )

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

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
#' For example map a model.id to patient.id
#'
#' @examples
#' data(brca)
#' mapModelSlotIds(brca, id="X-1004", id.name="patient.id", map.to="model.id")
#' ##map batch ids
#' mapModelSlotIds(brca, id="X-1004.BGJ398", id.name="batch.name", map.to="tissue")
#' @param object The \code{Xeva} dataset
#' @param id The \code{id}
#' @param id.name The \code{id} name
#' @param map.to The name of the mapped id. Default \code{all}
#' @param unique Default \code{TRUE}. If unique=FALSE output will be mapped to input
#' @return a \code{data.fram} with id and mapped id
#' @keywords internal
#' @noRd
### @export
mapModelSlotIds <- function(object, id, id.name, map.to="all", unique=TRUE)
{
  id <- c(as.character(id))

  if(id.name=="batch.name")
  {
   rtd <- .mapBatchName2Id(object, id, map.to)
  } else{
   .checkIfColPresentinModel(object, id.name)
   rtd <- object@model[object@model[,id.name] %in% id, ]
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
}

##--------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
##-----map batch to patient.id --------------------------------------------------------
##
.mapBatchName2Id <- function(object, id, map.to)
{
  btMapRet <- data.frame()
  for(bn in c(id))
  {
    #bt <- expDesign(object, batch.name = bn)
    bt <- batchInfo(object, batch = bn)
    bt.Mod <- unique(c(bt[[bn]]$treatment, bt[[bn]]$control))
    btMap <- mapModelSlotIds(object, id= bt.Mod, id.name="model.id",
                             map.to=map.to, unique=TRUE)

    btMap[, "batch.name"] <- bn
    btMapRet <- rbind(btMapRet, btMap)
  }
  btMapRet <- btMapRet[, c("batch.name", map.to)]
  btMapRet <- unique(btMapRet); rownames(btMapRet)=NULL
  return(btMapRet)
}

