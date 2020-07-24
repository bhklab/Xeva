##----- TGI -----
#TGI = (VC - VT)/(VC0 - VT0)
#where VC and VT are the median of control and treated growth curve respectively
#at the end of the study. VC0 and VT0 indicate the initial tumor volume for
#control and treated growth curve respectively.
#
#
#' tumor growth inhibition (TGI)
#' Computes the tumor growth inhibition (TGI) between two time-volume curves
#'
#' @param contr.volume Volume vector for control
#' @param treat.volume Volume vector for treatment
#' @return Returns batch response object
#'
#' @examples
#' contr.volume <- c(1.35, 6.57, 13.94, 20.39, 32.2, 39.26, 46.9, 53.91)
#' treat.volume <- c(0.4, 1.26, 2.59, 3.62, 5.77, 6.67, 7.47, 8.98, 9.29, 9.44)
#' TGI(contr.volume, treat.volume)
#'
#' @export
TGI <- function(contr.volume, treat.volume)
{
  tgi=NA
  if(length(contr.volume)>0 & length(treat.volume)>0)
  {
    tgi = (contr.volume[length(contr.volume)] - treat.volume[length(treat.volume)])/
          (contr.volume[1] - treat.volume[1])
  }
  rtx <- batch_response_class(name="TGI", value=tgi)
  return(rtx)
}

