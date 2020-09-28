.tumorVolumeChange <- function(volume)
{
  Vini <- volume[1]

  return(vapply(volume, function(vt){ 100*(vt-Vini)/Vini }, FUN.VALUE = numeric(1)) )
}

.avgResponse <- function(volume.change)
{
  ar <- cumsum(volume.change)/seq(along=volume.change)
  return(ar)
}


.checkNumericIntCharZero <- function(u)
{
  if(length(u)==0){return(NA)}
  return(u)
}

.getBestResponse <- function(time, response, min.time=NULL)
{
  rtz <- list()
  rtz[c("time", "value", "index")] <- NA

  exdf <-  data.frame(time= time, response=response)
  if(!is.null(min.time))
  {
    exdfMinAge = exdf[exdf$time >= min.time, ]
  } else { exdfMinAge <- exdf }

  if(dim(exdfMinAge)[1]==0)
  {
    return(rtz)
  }

  minIndxA <- which.min(exdfMinAge$response)
  minIndx  <- minIndxA + nrow(exdf) - nrow(exdfMinAge)

  rtz$time  <- .checkNumericIntCharZero(exdf[minIndx, "time"])
  rtz$value <- .checkNumericIntCharZero(exdf[minIndx, "response"])
  rtz$index <- .checkNumericIntCharZero(minIndx)
  return(rtz)
}


empty_model_response_class_mRECIST <- function()
{
  exdf <- list()
  exdf[c("volume.change", "average.response", "best.response", "best.response.time",
         "best.average.response", "best.average.response.time", "mRECIST")] <- NA

  model_response_class(name = "mRECIST", value = NA, fit=exdf)
}

######################################################################
#' Computes the mRECIST
#'
#' \code{mRECIST} Returns the mRECIST for given volume response.
#'
#' @param time Value of best response.
#' @param volume Value of best average response.
#' @param min.time Minimum time after which tumor volume will be considered.
#### @param return.detail Default \code{FALSE}. If \code{TRUE}, return all intermediate values.
#' @return  Returns the mRECIST.
#' @examples
#' time  <- c(0, 3, 7, 11, 18, 22, 26, 30, 32, 35)
#' volume<- c(250.8, 320.4, 402.3, 382.6, 384, 445.9, 460.2, 546.8, 554.3, 617.9)
#' mRECIST(time, volume, min.time=10, return.detail=FALSE)
#' @export
mRECIST <- function(time, volume, min.time=10 #, return.detail=FALSE
                    )
{
  if(volume[1]==0) {volume <- volume + 1}

  exdf <- list()
  exdf[c("volume.change", "average.response", "best.response", "best.response.time",
    "best.average.response", "best.average.response.time", "mRECIST")] <- NA

  exdf$volume.change <- .tumorVolumeChange(volume)
  exdf$average.response <- .avgResponse(exdf$volume.change)

  df <- data.frame(time=time, volume=volume)
  df <- df[df$time>=min.time, ]
  if(nrow(df)<2)
  {
    warning(sprintf("insufficient data after time %d", min.time))
  } else
  {
    br <- .getBestResponse(time, exdf$volume.change, min.time=min.time)
    exdf$best.response <- br$value
    exdf$best.response.time <- br$time

    bar <- .getBestResponse(time, exdf$average.response, min.time=min.time)
    exdf$best.average.response <- bar$value
    exdf$best.average.response.time <- bar$time

    best.response <- exdf$best.response
    best.average.response <- exdf$best.average.response

    mRecist <- NA; exdf$mRECIST <- mRecist

    if(!is.na(best.response) & !is.na(best.average.response))
    {
      ####---- the order of mRECIST assignment is important ----------
      mRecist <- "PD"

      if(best.response <  35 & best.average.response <  30)
      {mRecist <- "SD"}

      if(best.response < -50 & best.average.response < -20)
      {mRecist <- "PR"}

      if(best.response < -95 & best.average.response < -40)
      {mRecist <- "CR"}

      exdf$mRECIST <- mRecist
    }
  }

  #if(return.detail==FALSE)
  #{return(exdf$mRECIST)}
  #return(exdf)

  mr <- empty_model_response_class_mRECIST()
  mr$value <- exdf$mRECIST
  mr$fit <- exdf
  return(mr)
}


#' compute mRECIST for a batch
#' Computes the mRECIST for a batch.
#'
#' @param contr.time Time vector for control.
#' @param contr.volume Volume vector for control.
#' @param treat.time Time vector for treatment.
#' @param treat.volume Volume vector for treatment.
#' @param min.time Minimum time after which tumor volume will be considered.
#'
#' @return Returns batch response object.
#' @examples
#' contr.time <- treat.time  <- c(0, 3, 7, 11, 18, 22, 26, 30, 32, 35)
#' contr.volume<- contr.time * tan(60*pi/180)
#' treat.volume<- treat.time * tan(15*pi/180)
#' bmr <- bmRECIST(contr.time, contr.volume, treat.time, treat.volume)
#' print(bmr)
#' @export
bmRECIST <- function(contr.time=NULL, contr.volume=NULL,
                     treat.time=NULL, treat.volume=NULL,
                     min.time=10)
{
  con <- tre <- empty_model_response_class_mRECIST()

  if(!is.null(contr.time) & !is.null(contr.volume))
  {
    if(length(contr.volume)!=length(contr.time))
    {
      msg <- sprintf("contr.time and contr.volume should have same length")
      stop(msg)
    }
    con <- mRECIST(contr.time, contr.volume, min.time=min.time)
  }


  if(!is.null(treat.time) & !is.null(treat.volume))
  {
    if(length(treat.volume)!=length(treat.time))
    {
      msg <- sprintf("treat.time and treat.volume should have same length")
      stop(msg)
    }
    tre <- mRECIST(treat.time, treat.volume, min.time=min.time)
  }

  v <- sprintf("control = %s, treatment = %s", con$value, tre$value)
  rtx <- batch_response_class(name="bmRECIST", value=v, control=con, treatment=tre)
  return(rtx)
}
