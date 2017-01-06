
## calculate angle between control and treatment groups
#'
#' Calculate angle between control and treatment groups
#'
#'
#' Given a batch (control and treatment model ids)
#' it will return angle (in Degree) between
#' the linear fit of treatment and control group
#'
#' @examples
#' data(pdxe)
#' # creat a experiment desing
#' ExpDesign = list(batch.name="myBatch", treatment=c("X.015.BY19"), control=c("X.015.uned"))
#' df = calculateAngle(object=pdxe, ExpDesign, var = "volume", plot=TRUE)
#'
#' @param object The \code{Xeva} dataset
#' @param ExpDesign A list with batch.name, treatment and control
#' @param var Name of the variable, default \code{volume}
#' @return a \code{data.fram} with treatment, control and batch.name
setGeneric(name = "calculateAngle", def = function(object, ExpDesign, var="volume", plot=TRUE)
{standardGeneric("calculateAngle")} )

#' @export
setMethod( f=calculateAngle,
           signature=c(object="XevaSet", ExpDesign="list"),
           definition= function(object, ExpDesign, var="volume", plot=TRUE)
           {
             DFx = getTimeVarData(object, ExpDesign, var=var)
             dfC = DFx[DFx$exp.type=="control",]
             dfT = DFx[DFx$exp.type=="treatment",]

             if(nrow(dfC)==0)
             { stop("Control have no data!") }
             if(nrow(dfT)==0)
             { stop("treatment have no data!") }

             ##-------- do lm -------------------------------------------------
             f <- paste("mean", "~", paste("time", collapse=" + "))
             fitC = lm(f, data=dfC)
             fitT = lm(f, data=dfT)

             ##-----calculate angle ----------------------------------------
             angC = atan(coef(fitC)[["time"]]) *180 / pi
             angT = atan(coef(fitT)[["time"]]) *180 / pi

             angDiff = angC - angT

             if(plot==TRUE)
             {
               rng=range(DFx$mean)
               trng = range(DFx$time)
               opar <- par()      # make a copy of current settings
               par(pty="s",xpd=FALSE) ##important make plot region square
               plot(dfT$time,dfT$mean, col="red",  pch=19, xlab = "time", ylab = "tumor volume")
               points(dfC$time,  dfC$mean, col="blue", pch=19, xlim = trng, ylim = rng)
               abline(fitT, col="red")
               abline(fitC, col="blue")
               par(xpd=TRUE)
               legend(max(trng)+10, mean(rng),
                      legend=c("Treatment", "Control"), fill=c("red", "blue"), cex=0.8)
               par(pty=opar$pty, xpd=opar$xpd) ## reset par to old setting

             }
             return(angDiff)
           })





