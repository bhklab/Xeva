.trapz_AUC <- function (x, y)
{
  n <- 2:length(x)
  auc <- as.double((x[n]-x[n - 1]) %*% (y[n] + y[n - 1]))/2
  return(auc)
}


#' area under the curve
#' \code{AUC} Returns area under the curve
#'
#' @param time A \code{vector} of time points recorded for the experiment.
#' @param volume First \code{vector} of volume.
#' @param time.normalize If TRUE, AUC value will be divided by max time
#' @return  Returns \code{angle} and \code{slope} object.
#' @examples
#' time  <- c(0, 3, 7, 11, 18, 22, 26, 30, 32, 35)
#' volume1<- time * tan(30*pi/180)
#' volume2<- time * tan(45*pi/180)
#' auc1 <- AUC(time, volume1)
#' auc2 <- AUC(time, volume2)
#' par(pty="s")
#' xylimit <- range(c(time, volume1, volume2))
#' plot(time, volume1, type = "b", xlim = xylimit, ylim = xylimit, col="red")
#' abline(lm(volume1~time), col="red")
#' lines(time, volume2, type = "b", col="green")
#' abline(lm(volume2~time), col="green")
#' @export
AUC <- function(time, volume, time.normalize=TRUE, vol.normal=TRUE)
{
  if(vol.normal==TRUE)
  { volume <- .normalizeByElement1(volume) }
  
  if(time[1]!=0)
  { warning("time t0 is not zero")}
  auc <- .trapz_AUC(time, volume)
  if(time.normalize==TRUE){auc = auc/(time[length(time)] - time[1])}
  rtx <- model_response_class(name = "auc", value = auc)
  return(rtx)
}

#' area between curves
#' Computes the area between two time-volume curves.
#'
#' @param contr.time Time vector for control.
#' @param contr.volume Volume vector for control.
#' @param treat.time Time vector for treatment.
#' @param treat.volume Volume vector for treatment.
#'
#' @return Returns batch response object.
#'
#' @examples
#' contr.time <- treat.time  <- c(0, 3, 7, 11, 18, 22, 26, 30, 32, 35)
#' contr.volume<- contr.time * tan(60*pi/180)
#' treat.volume<- treat.time * tan(15*pi/180)
#' abc <- ABC(contr.time, contr.volume, treat.time, treat.volume)
#' par(pty="s")
#' xylimit <- range(c(contr.time, contr.volume, treat.time, treat.volume))
#' plot(contr.time, contr.volume, type = "b", xlim = xylimit, ylim = xylimit)
#' lines(treat.time, treat.volume, type = "b")
#' polygon(c(treat.time, rev(treat.time)), c(contr.volume, rev(treat.volume)),
#'         col = "#fa9fb5", border = NA)
#'
#' @export
ABC <- function(contr.time=NULL, contr.volume=NULL, treat.time=NULL, treat.volume=NULL)
{
  con <- tre <- model_response_class(name = "auc", value = NA)
  abc <- NA

  if(!is.null(contr.time) & !is.null(contr.volume))
  {
    if(length(contr.volume)!=length(contr.time))
    {
      msg <- sprintf("contr.time and contr.volume should have same length")
      stop(msg)
    }
    con <- AUC(contr.time, contr.volume)
  }


  if(!is.null(treat.time) & !is.null(treat.volume))
  {
    if(length(treat.volume)!=length(treat.time))
    {
      msg <- sprintf("treat.time and treat.volume should have same length")
      stop(msg)
    }
    tre <- AUC(treat.time, treat.volume)
  }

  abc <- con$value - tre$value

  rtx <- batch_response_class(name="abc", value=abc, control=con, treatment=tre)
  return(rtx)
}
