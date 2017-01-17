

.computSlopFun <- function(x,y)
{
  #fit = lm(y~x)

  ##-------------------------------------------------------------
  data= data.frame(x=x, y=y)
  p= x[1]; q = y[1]
  fit <- (lm(I(y-q)~I(x-p) +0, data))

  #nd = predict(fit, newdata = list(x=0))+q
  #.testPlot(x,y, nmod)
  #abline(nd, coef(nmod), col='blue')

  ##------------------------------------------------
  #ang = atan(coef(fit)[["x"]]) *180 / pi
  ang = atan(coef(fit)[["I(x - p)"]]) *180 / pi
  return(list(fit=fit, angel=ang, data=data))
  ##-------- do lm -------------------------------------------------
  #f <- paste("mean", "~", paste("time", collapse=" + "))
  #fitC = lm(f, data=dfC)
}


.plotAngelAndFit <- function(dfT, dfC, fitT, fitC)
{
  xrng = range(dfT$time, dfC$time)
  yrng = range(dfT$mean, dfC$mean)
  opar <- par()      # make a copy of current settings
  par(pty="s",  xpd=TRUE)
  plot(dfT$time, dfT$mean, col="red",  pch=19,
       xlab = "time", ylab = "tumor volume",
       xlim = xrng, ylim = yrng)
  points(dfC$time,  dfC$mean, col="blue", pch=19)

  legend("bottomright", inset=c(-0.38,0),
         legend=c("Treatment", "Control"), fill=c("red", "blue"), cex=0.8)

  par(xpd=FALSE);

  #nd = predict(fit, newdata = list(x=0))+q
  #.testPlot(x,y, nmod)
  #abline(nd, coef(nmod), col='blue')
  lxT = predict(fitT$fit, newdata = list(x=0))+dfT$mean[1]
  abline(lxT, coef(fitT$fit), col="red")

  lxC = predict(fitC$fit, newdata = list(x=0))+dfC$mean[1]
  abline(lxC, coef(fitC$fit), col="blue")

  #abline(fitT$fit, col="red"); abline(fitC$fit, col="blue")


  par(pty=opar$pty, xpd=opar$xpd) ## reset par to old setting
}


####-----------------------------------------------------------
plot_Batch_angel <- function(dt.fit, dc.fit, fitC, fitT, dfC, dfT)
{
  xrng = range(sapply(dt.fit, function(i) range(i$data$x)),
               sapply(dc.fit, function(i) range(i$data$x)))

  yrng = range(sapply(dt.fit, function(i) range(i$data$y)),
               sapply(dc.fit, function(i) range(i$data$y)))

  opar <- par()      # make a copy of current settings
  par(pty="s",  xpd=TRUE)
  plot(NA, col="red",  pch=19,
       xlab = "time", ylab = "tumor volume",
       xlim = xrng, ylim = yrng)

  .plt1Mod <- function(d.fit, col){
  for(dfX in d.fit)
  {
    #points(dfX$data$x,  dfX$data$y, col="blue", pch=19)
    #lx = predict(dfX$fit, newdata = list(x=0))+dfX$data$y[1]
    #abline(lx, coef(dfX$fit), col="blue")
    lines(dfX$data$x,  dfX$data$y, col=col,type="b", lty=3, pch=19)
  } }

  .plt1Mod(dc.fit, col = "#6baed6")
  .plt1Mod(dt.fit, col = "#fc8d59")

  ##----- add line ------
  if(1==2){
  dfX =dc.fit[[5]]
  points(dfX$data$x,  dfX$data$y, col="blue", pch=19)
  lx = predict(dfX$fit, newdata = list(x=0))+dfX$data$y[1]
  yl = coef(dfX$fit)
  #avgAng = mean(sapply(dc.fit, "[[", "angel"))
  avgAng = mean(sapply(dc.fit, function(x)coef(x$fit)))
  yv = c(avgAng*xrng[1], avgAng*xrng[2])
  lines(xrng, yv, col="red")
  #lx1 = predict(dfX$fit, newdata = list(x=xrng[1]))+dfX$data$y[1]
  #lx2 = predict(dfX$fit, newdata = list(x=xrng[2]))+dfX$data$y[1]
  #abline(lx, coef(dfX$fit), col="blue")
  #lines(c(0,110), c(4.006887*0+lx, 4.006887*110+lx), col="red")
  }


  legend("bottomright", inset=c(-0.38,0),
         legend=c("Treatment", "Control"),
         fill=c("#a50f15", "#081d58"), cex=0.8)

  par(xpd=FALSE)

  lxT = predict(fitT$fit, newdata = list(x=0))+dfT$mean[1]
  abline(lxT, coef(fitT$fit), col="#a50f15")

  lxC = predict(fitC$fit, newdata = list(x=0))+dfC$mean[1]
  abline(lxC, coef(fitC$fit), col="#081d58")

  par(pty=opar$pty, xpd=opar$xpd)

}



####-----------------------------------------------------------
#' expDegI = list(batch.name = "myBatch", treatment = "X.015.BY19", control = "X.015.uned")
computAngelFor1ExpDesign <- function(object, expDegI, var="volume", treatment.only=TRUE, plot=FALSE)
{
  ##------ for each model get fit ---------
  dt.fit= list()
  if(length(expDegI$treatment)>0)
  {
    for(ti in expDegI$treatment)
    {
      dt = getExperiment(object, ti, treatment.only = TRUE)
      dt.fit[[ti]] = .computSlopFun(dt$time, dt[,var])
    }
  }
  ###-------------
  dc.fit= list()
  if(length(expDegI$control)>0)
  {
    for(ci in expDegI$control)
    {
      dc = getExperiment(object, ci, treatment.only = TRUE)
      dc.fit[[ci]] = .computSlopFun(dc$time, dc[,var])
    }
  }

  DFx = getTimeVarData(object, expDegI, var=var, treatment.only=treatment.only)


  dfC = DFx[DFx$exp.type=="control",]
  dfT = DFx[DFx$exp.type=="treatment",]

  if(nrow(dfC)==0)
  {
    warning("Control have no data!")
    return(NA)
  }

  if(nrow(dfT)==0)
  {
    warning("treatment have no data!")
    return(NA)
  }
  fitC = .computSlopFun(dfC$time, dfC$mean)
  fitT = .computSlopFun(dfT$time, dfT$mean)
  angDiff = fitC$angel - fitT$angel

  if(plot==TRUE)
  { plot_Batch_angel(dt.fit, dc.fit, fitC, fitT,dfC, dfT)}

  return(angDiff)
}





computAngelFor1ExpDesign_old <- function(object, expDegI, var="volume", treatment.only=TRUE, plot=FALSE)
{
  DFx = getTimeVarData(object, expDegI, var=var, treatment.only=treatment.only)


  dfC = DFx[DFx$exp.type=="control",]
  dfT = DFx[DFx$exp.type=="treatment",]

  if(nrow(dfC)==0)
  {
    warning("Control have no data!")
    return(NA)
  }

  if(nrow(dfT)==0)
  {
    warning("treatment have no data!")
    return(NA)
  }
  fitC = .computSlopFun(dfC$time, dfC$mean)
  fitT = .computSlopFun(dfT$time, dfT$mean)

  angDiff = fitC$angel - fitT$angel

  if(plot==TRUE)
  { .plotAngelAndFit(dfT, dfC, fitT, fitC) }

  return(angDiff)
}


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
#' df = calculateAngle(object=pdxe, ExpDesign, var = "volume", treatment.only=TRUE, plot=TRUE)
#'
#' @param object The \code{Xeva} dataset
#' @param ExpDesign A list with batch.name, treatment and control
#' @param var Name of the variable, default \code{volume}
#' @return a \code{data.fram} with treatment, control and batch.name
setGeneric(name = "calculateAngle", def = function(object, ExpDesign=NULL, var="volume", treatment.only=TRUE, plot=TRUE)
{standardGeneric("calculateAngle")} )

#' @export
setMethod( f=calculateAngle,
           signature=c(object="XevaSet"),
           definition= function(object, ExpDesign=NULL, var="volume", treatment.only=TRUE, plot=TRUE)
           {
             if(!is.null(ExpDesign))
             { ExpDesign = list(ExpDesign) }

             if(is.null(ExpDesign))
             { ExpDesign = expDesignInfo(object) }

             angDiffX =list()
             for(I in 1:length(ExpDesign))
             {
               expDegI = ExpDesign[[I]]
               angDiffX[[expDegI$batch.name]] = computAngelFor1ExpDesign(object, expDegI, var=var,
                                                                     treatment.only=treatment.only, plot=plot)
             }
             return(angDiffX)
           })




#####================= setAngle ==================
#' setAngle generic
#'
#' Generic for setAngle method
#'
#' @examples
#' data(pdxe)
#' angles <- setAngle(object = pdxe)
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{list} with angle between treatment and control for each batch
setGeneric(name = "setAngle", def = function(object) {standardGeneric("setAngle")} )

#' @export
setMethod( f=setAngle, signature="XevaSet",
           definition=function(object)
           { calculateAngle(object, var = "volume", treatment.only=TRUE, plot=FALSE) } )


#' setAngle<- Generic
#'
#' Generic for setAngle replace method
#'
#' @examples
#' data(pdxe)
#' setAngle(pdxe) <- setAngle(object = pdxe)
#' @param object The \code{XevaSet} to replace drug info in
#' @return Updated \code{XevaSet}
setGeneric(name= "setAngle<-", def = function(object, value) {standardGeneric("setAngle<-")} )

#' @export
setMethod( f="setAngle<-",
           signature=c(object = "XevaSet"),
           definition=function(object, value)
           {
             #object@model = value
             for(I in 1:length(expDesignInfo(object)))
             {
               object@expDesign[[I]]$angle = value[[object@expDesign[[I]]$batch.name]]
             }
             return(object)
           } )





#####================= setSlop ==================
extractBetweenTags <- function(inVec, start.tag=0, end.tag=0)
{
  inVIndx= 1:length(inVec)
  stIndx = min(inVIndx[inVec!=start.tag])

  V2 = inVec[stIndx:length(inVec)]
  v2end = which(V2==end.tag)
  if(length(v2end)>0)
  {
    enIndx = min(v2end) -1
    enIndxR= enIndx + stIndx -1
  } else
    {enIndxR = length(inVec)}

  Vi = stIndx:enIndxR
  return(Vi)
}



##------- Compute slop for one model ------------
computeSlope <- function(object, model.id, treatment.only=TRUE)
{
  mod = getExperiment(object, model.id)
  time = mod$time; volume = mod$volume
  if(treatment.only==TRUE)
  {
    if(!is.null(mod$dose))
    {
      tretIndx = extractBetweenTags(mod$dose, start.tag=0, end.tag=0)
      time = time[tretIndx] ; volume=volume[tretIndx]
    }
  }
  rtx = .computSlopFun(time, volume)
  return(rtx)
}


#' setSlop generic
#'
#' Generic for setSlop method
#'
#' @examples
#' data(pdxe)
#' slops <- setSlop(object = pdxe, treatment.only=FALSE)
#' @param object The \code{XevaSet} to retrieve drug info from
#' @param treatment.only Default \code{TRUE}, take only first periode of treatment where dose is not zero
#' @details If dose column is not present treatment.only will be considered FALSE
#' @return a \code{list} with slop of each model.id
setGeneric(name = "setSlop", def = function(object, treatment.only=TRUE) {standardGeneric("setSlop")} )

#' @export
setMethod( f="setSlop", signature="XevaSet",
           definition=function(object, treatment.only=TRUE)
           {
             rtz = list()
             for(model.id in names(object@experiment))
             {
               fa = computeSlope(object, model.id, treatment.only=treatment.only)
               rtz[[model.id]] = fa$angel
             }
             return(rtz)
           })

#' setSlop<- Generic
#'
#' Generic for setSlop replace method
#'
#' @examples
#' data(pdxe)
#' setSlop(pdxe) <- setSlop(object = pdxe)
#' @param object The \code{XevaSet} to replace drug info in
#' @return Updated \code{XevaSet}
setGeneric(name= "setSlop<-", def = function(object, value) {standardGeneric("setSlop<-")} )

#' @export
setMethod( f="setSlop<-",
           signature=c(object = "XevaSet"),
           definition=function(object, value)
           {
             for(model.id in names(object@experiment))
             {
               object@experiment[[model.id]]$slop = value[[model.id]]
             }
             return(object)
           })

