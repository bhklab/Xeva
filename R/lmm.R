#' linear mixed model
#'
#' Comput the linear mixed model (lmm) statistics for a PDX batch
#'
#' @param data a data.frame containing a batch data
#' @param log.volume FALSE if volume is raw, TRUE if volume is in log
#' @return Returns a fit object
#'
#' @details The input data.frame (data) must contain these columns: 'model.id', 'volume', 'time', 'exp.type'.
#' Lower value from lmm function indicates better response.
#'
#'
#' @examples
#' data(repdx)
#' data <- getExperiment(repdx, batch = "P1")$model
#' lmm(data, log.volume=FALSE)
#'
#' data3 <- getExperiment(repdx, batch = "P3")$model
#' lmm(data3, log.volume=FALSE)
#'
#' @export
#' @import nlme
lmm <- function(data, log.volume)
{
  if(any(!c("model.id", "volume", "time", "exp.type")%in% colnames(data)))
  {
    msg="these columns must be present, 'model.id', 'volume', 'time', 'exp.type'"
    stop(msg)
  }

  if(!log.volume %in% c(TRUE, FALSE))
  {
    msg="log.volume should be either TRUE if volume is in log or FALSE if volume is raw"
    stop(msg)
  }

  if(log.volume==FALSE)
  {
    minvol <- min(data$volume)
    if(minvol<=0){data$volume <- data$volume + abs(minvol) + 1}
    data$volume <- log(data$volume)
  }

  fit <- lme(volume~time*exp.type, data=data, random= ~1|model.id)
  fit$value <- as.numeric(fit$coefficients$fixed[4])
  return(fit)
}
