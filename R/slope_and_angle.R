#' Computes slope
#'
#' \code{slope} returns the slope for given time and volume data
#'
#' @param time vector of time
#' @param volume vector of volume
#' @param degree default \code{TRUE} will give angle in Degree and \code{FALSE} will return Radians
#' @return  returns the slope and a fit object
#' @examples
#' time  <- c(0, 3, 7, 11, 18, 22, 26, 30, 32, 35)
#' volume<- c(250.8, 320.4, 402.3, 382.6, 384, 445.9, 460.2, 546.8, 554.3, 617.9)
#' sl <- slope(time, volume)
#' par(pty="s")
#' xylimit <- range(c(time, volume))
#' plot(time, volume, type = "b", xlim = xylimit, ylim = xylimit)
#' abline(lm(volume~time))
#' @export
slope <- function(time, volume, degree=TRUE)
{
  df <- data.frame(time=time, volume=volume)
  ##---- remove all non finite (Inf, NA, NaN) data --------------
  df <- df[is.finite(df$time), ]
  df <- df[is.finite(df$volume),]

  df$time  <- df$time- df$time[1]
  df$volume<- df$volume- df$volume[1]

  fit <- lm(volume~time +0, df)
  ang <- atan(coef(fit)[["time"]])
  #z <- sum(df$time*df$volume) / (sqrt(sum(df$time * df$time)) * sqrt(sum(df$volume * df$volume)) )
  #ang <- acos(z)
  if(degree==TRUE) { ang <- ang*180/base::pi }

  rtx <- model_response_class(name = "slope", value = ang, fit=fit)
  return(rtx)
}


#'
#' \code{angle} returns the angle between two volume data
#'
#' @param contr.time time vector for control
#' @param contr.volume volume vector for control
#' @param treat.time time vector for treatment
#' @param treat.volume volume vector for treatment
#' @param degree default \code{TRUE} will give angle in Degree and \code{FALSE} will return Radians
#' @return  returns angle and \code{slope}
#' @examples
#' contr.time <- treat.time  <- c(0, 3, 7, 11, 18, 22, 26, 30, 32, 35)
#' contr.volume<- contr.time * tan(60*pi/180)
#' treat.volume<- treat.time * tan(15*pi/180)
#' ang <- angle(contr.time, contr.volume, treat.time, treat.volume)
#' par(pty="s")
#' xylimit <- range(c(time, volume1, time, volume2))
#' plot(time, volume1, type = "b", xlim = xylimit, ylim = xylimit)
#' lines(time, volume2, type = "b")
#' abline(lm(volume1~time))
#' abline(lm(volume2~time))
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

  #if(class(con)=="modelResponse" & class(tre)=="modelResponse")
  ang <- con$value - tre$value
  rtx <- batch_response_class(name="angle", value=ang, control=con, treatment=tre)
  return(rtx)
}
