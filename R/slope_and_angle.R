#' Computes slope
#'
#' \code{slope} returns the slope for given time and volume data.
#'
#' @param time A \code{vector} of time.
#' @param volume A \code{vector} of volume.
#' @param degree Default \code{TRUE} will give angle in degrees and \code{FALSE} will return in radians.
#' @return Returns the slope and a \code{fit} object.
#' @examples
#' time  <- c(0, 3, 7, 11, 18, 22, 26, 30, 32, 35)
#' volume<- c(250.8, 320.4, 402.3, 382.6, 384, 445.9, 460.2, 546.8, 554.3, 617.9)
#' sl <- slope(time, volume)
#' par(pty="s")
#' xylimit <- range(c(time, volume))
#' plot(time, volume, type = "b", xlim = xylimit, ylim = xylimit)
#' abline(lm(volume~time))
#' @export
slope <- function(time, volume, degree=TRUE, vol.normal=TRUE)
{
  df <- data.frame(time=time, volume=volume)
  ##---- remove all non finite (Inf, NA, NaN) data --------------
  df <- df[is.finite(df$time), ]
  df <- df[is.finite(df$volume),]

  #df$time  <- df$time- df$time[1]
  #df$volume<- df$volume-df$volume[1]
  if(df$time[1]!=0)
  { warning("time t0 is not zero")}
  
  if(vol.normal==TRUE)
  { df$volume <- .normalizeByElement1(df$volume) }
  
  fit <- lm(volume~time +0, df)
  ang <- atan(coef(fit)[["time"]])
  ##----old way to compute angle ---
  #z <- sum(df$time*df$volume) / (sqrt(sum(df$time * df$time)) * sqrt(sum(df$volume * df$volume)) )
  #ang <- acos(z)
  if(degree==TRUE) { ang <- ang*180/base::pi }

  rtx <- model_response_class(name = "slope", value = ang, fit=fit)
  return(rtx)
}


#' compute angle
#' Computes the angle between two time-volume curves.
#'
#' @param contr.time Time vector for control.
#' @param contr.volume Volume vector for control.
#' @param treat.time Time vector for treatment.
#' @param treat.volume Volume vector for treatment.
#' @param degree Default \code{TRUE} will give angle in degrees and \code{FALSE} will return in radians.
#' @return Returns batch response object.
#' @examples
#' contr.time <- treat.time  <- c(0, 3, 7, 11, 18, 22, 26, 30, 32, 35)
#' contr.volume<- contr.time * tan(60*pi/180)
#' treat.volume<- treat.time * tan(15*pi/180)
#' ang <- angle(contr.time, contr.volume, treat.time, treat.volume)
#' print(ang)
#' par(pty="s")
#' xylimit <- range(c(contr.time, contr.volume, treat.time, treat.volume))
#' plot(contr.time, contr.volume, type = "b", xlim = xylimit, ylim = xylimit)
#' lines(treat.time, treat.volume, type = "b")
#' abline(lm(contr.volume~contr.time))
#' abline(lm(treat.volume~treat.time))
#' @export
angle <- function(contr.time=NULL, contr.volume=NULL, treat.time=NULL, treat.volume=NULL,
                  degree=TRUE)
{
  con <- tre <- model_response_class(name = "slope", value = NA, fit=NA)
  ang <- NA

  if(!is.null(contr.time) & !is.null(contr.volume))
  {
    if(length(contr.volume)!=length(contr.time))
    {
      msg <- sprintf("contr.time and contr.volume should have same length")
      stop(msg)
    }
    con <- slope(contr.time, contr.volume, degree= degree)
  }


  if(!is.null(treat.time) & !is.null(treat.volume))
  {
    if(length(treat.volume)!=length(treat.time))
    {
      msg <- sprintf("treat.time and treat.volume should have same length")
      stop(msg)
    }
    tre <- slope(treat.time, treat.volume, degree= degree)
  }

  #if(c#lass(con)=="modelResponse" & c#lass(tre)=="modelResponse") ##old Class command
  ang <- con$value - tre$value
  rtx <- batch_response_class(name="angle", value=ang, control=con, treatment=tre)
  return(rtx)
}
