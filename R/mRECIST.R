##--- to calculate mRECIST --------

tumorVolumeChange <- function(volume)
{
  Vini = volume[1]
  return( sapply(volume, function(vt){ 100*(vt-Vini)/Vini }) )
}

avgResponse <- function(volume.change)
{
  ar = cumsum(volume.change)/seq(along=volume.change)
  return(ar)
}


checkNumericIntCharZero <- function(u)
{
  if(length(u)==0){return(NA)}
  return(u)
}

getBestResponse <- function(exdf, ResColName, min.time=10)
{
  exdfMinAge = exdf[exdf$time >= min.time, ]
  if(dim(exdfMinAge)[1]==0)
    {exdfMinAge = exdf}
  minIndxA = which.min( exdfMinAge[, ResColName] )
  minIndx = minIndxA + dim(exdf)[1] - dim(exdfMinAge)[1]

  rtz = list(time = exdf[minIndx, "time"],
             value= exdf[minIndx, ResColName],
             index= minIndx)

  rtz$time  = checkNumericIntCharZero(rtz$time)
  rtz$value = checkNumericIntCharZero(rtz$value)
  rtz$index = checkNumericIntCharZero(rtz$index)
  return(rtz)
}

calculateResponses <- function(exdf, responseName = "volume", min.time=10)
{
  exdf$volume.change = tumorVolumeChange(exdf[,responseName])
  exdf$average.response= avgResponse(exdf$volume.change)

  best.response = getBestResponse(exdf, ResColName ="volume.change", min.time=min.time)
  best.average.response = getBestResponse(exdf, ResColName ="average.response", min.time=min.time)

  return(list(data = exdf,
              best.response=best.response,
              best.average.response=best.average.response))
}


##check PDTX curve

######################################################################
#' Computes the mRECIST
#'
#' \code{computemRECIST} returns the mRECIST for given volume response
#'
#' @param best.response Value of best response
#' @param best.average.response Value of best average response
#'
#' @return  Returns the mRECIST for given volume response
#'
#' @examples
#' computemRECIST(best.response=8.722, best.average.response=8.722)
#' @export
computemRECIST <- function(best.response, best.average.response)
{
  mRecist = NA
  if(is.na(best.response) | is.na(best.average.response) )
  {return(mRecist)}

  if(!is.na(best.response) & !is.na(best.average.response))
  {
    # if(best.response < -95 & best.average.response < -40)
    # {mRecist = "mCR"}
    #
    # if(best.response < -50 & best.average.response < -20)
    # {mRecist = "mPR"}
    #
    # if(best.response <  35 & best.average.response <  30)
    # {mRecist = "mSD"}
    #
    # if(is.na(mRecist)){mRecist = "mPD"}

    ####---- the order of aRecist assignment is really important ----------
    mRecist = "PD"

    if(best.response <  35 & best.average.response <  30)
    {mRecist = "SD"}

    if(best.response < -50 & best.average.response < -20)
    {mRecist = "PR"}

    if(best.response < -95 & best.average.response < -40)
    {mRecist = "CR"}
  }

  return(mRecist)
}



mRECISTForModel <- function(modx)
{
  if(is.null(modx$best.response$value))
  {
    modxDataMat = calculateResponses(modx$data, responseName = "volume")
    modx$data = modxDataMat$data
    modx$best.response = modxDataMat$best.response
    modx$best.average.response = modxDataMat$best.average.response
  }

  if(is.null(modx$data$body.weight.change))
  {
    modx$data$body.weight.change = tumorVolumeChange(modx$data$body.weight)
  }

  modx$mRECIST = computemRECIST(best.response = modx$best.response$value,
                                best.average.response= modx$best.average.response$value)
  return(modx)
}


###===============================================================================================
###===============================================================================================

#' setmRECIST<- Generic for setmRECIST replace method
#' @examples
#' data(pdxe)
#' #calculate mRECIST for each experiment
#' setmRECIST(pdxe)<- setmRECIST(pdxe)
#' getmRECIST(pdxe)
#' @param object The \code{XevaSet} object
#' @return Updated \code{XevaSet}
setGeneric(name= "setmRECIST", def = function(object) {standardGeneric("setmRECIST")} )

#' @export
setMethod( f="setmRECIST",
           signature = "XevaSet",
           definition= function(object)
           {
             if(is(object, "XevaSet"))
             {
               for(I in names(object@experiment))
               {
                 object@experiment[[I]] = mRECISTForModel(object@experiment[[I]])
               }
             }
             return(object)
           } )


setGeneric(name= "setmRECIST<-", def = function(object,value) {standardGeneric("setmRECIST<-")} )

## @describeIn PharmacoSet Returns the annotations for all the drugs tested in the PharmacoSet
#' @export
setMethod( f="setmRECIST<-",
           signature= signature(object="XevaSet"),
           definition=function(object, value)
           {
             object = value
             object
           } )


###===============================================================================================
###===============================================================================================
#' getmRECIST Generic
#' Generic for getmRECIST method
#'
#' @examples
#' data(pdxe)
#' # calculate mRECIST for each experiment
#' setmRECIST(pdxe)<- setmRECIST(pdxe)
#' getmRECIST(pdxe, group.by="biobase.id")
#' @param object The \code{XevaSet} to retrieve mRECIST from
#' @param group.by The name of column which will be mapped to model.id
#' @return a \code{data.frame} with the mRECIST values, rows are drugs and columns are model.id
setGeneric(name = "getmRECIST", def = function(object, group.by="biobase.id") {standardGeneric("getmRECIST")} )


#' @export
setMethod( f=getmRECIST,
           signature="XevaSet",
           definition= function(object, group.by)
           {
             rtx = data.frame(matrix(NA, nrow = length(object@experiment), ncol = 3))
             colnames(rtx) = c("drug.join.name", "model.id", "mRECIST")
             dfI = 1
             for(I in object@experiment)
             {
               if(is.null(I$mRECIST))
                 {
                 msg = sprintf("mRECIST not present for model %s\nRun setmRECIST(object)<- setmRECIST(object) first\n", I$model.id)
                 stop(msg)
                 }
              rtx[dfI, ] <- c(I$drug$join.name, I$model.id, I$mRECIST)
              dfI = dfI+1
             }
             rownames(rtx)= NULL

             ##----map to patient id -----------------
             #rtx[, group.by] = subset(object@model, object@model$model.id %in% rtx$model.id)[,group.by]
             rtx = merge(rtx, object@model[, c("model.id", group.by)], by.x = "model.id", by.y = "model.id")

             dataColName = c(group.by, "model.id", "drug.join.name", "mRECIST")
             rtx = BBmisc::sortByCol(rtx , dataColName, asc = rep(TRUE, length(dataColName)))
             return(rtx[,dataColName])
           }
           )






