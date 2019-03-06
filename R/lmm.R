# library(nlme)
#
# repdx <- readRDS("~/CXP/Xeva_dataset/data/LungPDX_Xeva.rds")
# save(repdx)
#
#
# x <- sample(1000)
# devtools::use_data(repdx)
#
#
# d <- getExperiment(lpdx, batch = "P1")
# data=d$model


#' linear mixed model
#'
#' Comput the linear mixed model (lmm) statistics for a PDX batch
#'
#' @param data a data.frame containg a batch data
#' @return Returns a fit object
#'
#' @details The input data.frame (data) must contain these columns: model.id, volume, time, exp.type
#'
#' @examples
#' data(repdx)
#' data <- getExperiment(repdx, batch = "P1")$model
#' lmm(data)
#'
#' @export
#' @import nlme
lmm <- function(data)
{
  if(any(!c("model.id", "volume", "time", "exp.type")%in% colnames(data)))
  {
    msg="these columns must be present, 'model.id', 'volume', 'time', 'exp.type'"
    stop(msg)
  }

  fit <- lme(log(volume)~time*exp.type, data=data, random= ~1|model.id)
  fit$value <- as.numeric(fit$coefficients$fixed[4])
  return(fit)
}
