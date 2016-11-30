##--- to calculate mRECIST --------

tumorVolumeChange <- function(volume)
{
  Vini = volume[1]
  return( sapply(exdf$volume, function(vt){ 100*(vt-Vini)/Vini }) )
}

avgResponse <- function(volume.change)
{
  ar = cumsum(volume.change)/seq(along=volume.change)
  return(ar)
}

getBestResponse <- function(exdf, ResColName, min.time=10)
{
  exdfMinAge = exdf[exdf$time >= min.time, ]
  minIndxA = which.min( exdfMinAge[, ResColName] )
  minIndx = minIndxA + dim(exdf)[1] - dim(exdfMinAge)[1]
  return(list(time = exdf[minIndx, "time"],
              value= exdf[minIndx, ResColName],
              index= minIndx))
}

calculateResponses <- function(exdf, responseName = "volume")
{
  exdf$volume.change = tumorVolumeChange(exdf[,responseName])
  exdf$average.response= avgResponse(exdf$volume.change)

  best.response = getBestResponse(exdf, ResColName ="volume.change", min.time=10)
  best.average.response = getBestResponse(exdf, ResColName ="average.response", min.time=10)
  return(list(data = exdf,
              best.response=best.response,
              best.average.response=best.average.response))
}


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
  if(best.response < -95 & best.average.response < -40)
  {mRecist = "mCR"}

  if(best.response < -50 & best.average.response < -20)
  {mRecist = "mPR"}

  if(best.response <  35 & best.average.response <  30)
  {mRecist = "mSD"}

  if(is.na(mRecist)){mRecist = "mPD"}

  return(mRecist)
}


getmRECIST <- function(exdf, responseName = "volume")
{
  #modx = pdxe@experiment[[4]] #$data
  if(is.null(modx$best.response$value))
  {
    modxDataMat = calculateResponses(modx$data, responseName = "volume")
    modx$data = modxDataMat$data
    modx$best.response = modxDataMat$best.response
    modx$best.average.response = modxDataMat$best.average.response
  }

  modx$mRECIST = computemRECIST(best.response = modx$best.response$value,
                               best.average.response= modx$best.average.response$value)

}








