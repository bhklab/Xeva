

.computSlopFun <- function(x,y, log.y=TRUE)
{
  #fit = lm(y~x)

  if(log.y==TRUE)
  {y = log(y)}

  #print(log(y))
  ##-------------------------------------------------------------
  data = data.frame(x=x, y=y)
  ##---- remove all non finite (Inf, NA, NaN) data --------------
  data = data[is.finite(data$x), ]
  data = data[is.finite(data$y), ]

  p= data$x[1]; q = data$y[1]
  fit <- (lm(I(y-q)~I(x-p) +0, data))
  ##------------------------------------------------
  #ang = atan(coef(fit)[["x"]]) *180 / pi
  ang = atan(coef(fit)[["I(x - p)"]]) *180/3.141593
  return(list(fit=fit, angel=ang, data=data))
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
plot_Batch_angel <- function(dt.fit, dc.fit, fitC, fitT, dfC, dfT, title="plot")
{
  xrng = range(sapply(dt.fit, function(i) range(i$data$x)),
               sapply(dc.fit, function(i) range(i$data$x)))

  yrng = range(sapply(dt.fit, function(i) range(i$data$y)),
               sapply(dc.fit, function(i) range(i$data$y)))

  opar <- par()      # make a copy of current settings
  par(pty="s",  xpd=TRUE)
  plot(NA, col="red",  pch=19, main=title,
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

  tPoint0 = mean(sapply(dt.fit, function(i) i$data$y[1]))

  lxT = predict(fitT$fit, newdata = list(x=0))+ tPoint0 # dfT$y[1]
  abline(lxT, coef(fitT$fit), col="#a50f15")

  cPoint0 = mean(sapply(dc.fit, function(i) i$data$y[1]))
  lxC = predict(fitC$fit, newdata = list(x=0))+ cPoint0 # dfC$y[1]
  abline(lxC, coef(fitC$fit), col="#081d58")

  par(pty=opar$pty, xpd=opar$xpd)

}

##--------------------------------------------------------------------

#' @import ggplot2
plot_Batch_angel_ggplot <- function(dt.fit, dc.fit, fitC, fitT, title="plot")
{

  ##-----make one DF ----------------------------------------
  dft = do.call(rbind, lapply(names(dt.fit), function(n)
                       {d=dt.fit[[n]]$data; d$model.id=n; d } ))
  dft$type = "treatment"

  dfc = do.call(rbind, lapply(names(dc.fit), function(n)
                       {d=dc.fit[[n]]$data; d$model.id=n; d } ))
  dfc$type = "control"

  DF = rbind(dft, dfc)
  tcCol <- c("control" = "#6baed6", "treatment" = "#fc8d59")

  plt <- ggplot(DF, aes_string(x="x", y="y", color= "type", group="model.id"))
  plt <- plt + geom_line(linetype = 2)+ geom_point()
  plt <- plt + scale_color_manual(values=tcCol)

  ##-------add lm line -----------------------------------------------------------
  addLMfitLines <- function(plt, fit, p0, color="black")
  {
    lx = predict(fit, newdata = list(x=0))+ p0
    plt + geom_abline(intercept = lx, slope = coef(fit)[1], color=color)
  }

  plt <- addLMfitLines(plt, fit = fitC$fit, p0= fitC$data$y[1], color="blue")
  plt <- addLMfitLines(plt, fit = fitT$fit, p0= fitT$data$y[1], color="red")
  ##-------------------------------------------------------------------------------
  plt <- plt + theme(aspect.ratio=1)
  plt <- .ggplotEmptyTheme(plt)
  plt <- plt + labs(title = title, x = "time", y = "log(volume)", colour = "")
  plt <- plt + theme(plot.title = element_text(hjust = 0.5))
  plt <- plt + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  #p + theme(legend.position='none')
  return(plt)
}






####-----------------------------------------------------------
#' data(lpdx); object=lpdx
#' expDegI  <- expDesign(lpdx, "PHLC111_P7")
computAngelFor1ExpDesign <- function(object, expDegI, var="volume", treatment.only=TRUE,
                                     plot=FALSE, log.y=TRUE)
{
  ##------ for each model get fit ---------
  dt.fit= list()
  if(length(expDegI$treatment)>0)
  {
    for(ti in expDegI$treatment)
    {
      dt = getExperiment(object, ti, treatment.only = TRUE)

      dt.fit[[ti]] = .computSlopFun(dt$time, dt[,var], log.y=log.y)
    }
  }
  ###-------------
  dc.fit= list()
  if(length(expDegI$control)>0)
  {
    for(ci in expDegI$control)
    {
      dc = getExperiment(object, ci, treatment.only = TRUE)
      dc.fit[[ci]] = .computSlopFun(dc$time, dc[,var], log.y=log.y)
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
  fitC = .computSlopFun(dfC$time, dfC$mean, log.y=log.y)
  fitT = .computSlopFun(dfT$time, dfT$mean, log.y=log.y)
  angDiff = fitC$angel - fitT$angel

  if(plot==TRUE)
  {
    #plot_Batch_angel(dt.fit, dc.fit, fitC, fitT,dfC, dfT, title=expDegI$batch.name)
    plt <- plot_Batch_angel_ggplot(dt.fit, dc.fit, fitC, fitT, title=expDegI$batch.name)
    #print(plt)
    return(list(angle= angDiff, plot=plt))
  }

  return(angDiff)
}






## calculate angle between control and treatment groups
#'
#' Calculate angle between control and treatment groups
#' @description
#' Given a batch (control and treatment model ids)
#' it will return angle (in Degree) between
#' the linear fit of treatment and control group
#'
#' @examples
#' data(pdxe)
#' # creat a experiment desing
#' ExpDesign = list(batch.name="myBatch", treatment=c("X.015.BY19"), control=c("X.015.uned"))
#' angl = calculateAngle(object=pdxe, ExpDesign, var = "volume", treatment.only=TRUE, plot=TRUE)
#' #print angle
#' print(angl$angle)
#' #print plot
#' print(angl$myBatch$plot)
#' #print without legend
#' print(angl$myBatch$plot+ theme(legend.position='none'))
#' @param object The \code{Xeva} dataset
#' @param ExpDesign A list with batch.name, treatment and control
#' @param var Name of the variable, default \code{volume}
#' @param treatment.only default {TRUE}. If TRUE only treatment periode will be considered
#' @param plot default \code{TRUE} will return a list with element \code{angle} value and plot
#' @param log.y default \code{TRUE} will take log of the values before linear regression. After taking log all non-numeric (-Inf, Inf, NaN) values will be ignored
#'
#' @details If plot==TRUE, this will return a list with angle and plot. The plot is made using ggplot and can be ggplot customized for aesthetics. See example for details
#'
#' @return a \code{data.fram} with treatment, control and batch.name
setGeneric(name = "calculateAngle",
           def = function(object, ExpDesign=NULL, var="volume", treatment.only=TRUE,
                          plot=TRUE, log.y=TRUE)
                          {standardGeneric("calculateAngle")} )

#' @export
setMethod( f=calculateAngle,
           signature=c(object="XevaSet"),
           definition= function(object, ExpDesign=NULL, var="volume", treatment.only=TRUE,
                                plot=TRUE, log.y=TRUE)
           {
             if(!is.null(ExpDesign))
             { ExpDesign = list(ExpDesign) }

             if(is.null(ExpDesign))
             { ExpDesign = expDesignInfo(object) }

             angDiffX =list()
             for(I in 1:length(ExpDesign))
             {
               expDegI = ExpDesign[[I]]
               angDiffX[[expDegI$batch.name]] = computAngelFor1ExpDesign(object, expDegI, var=var, treatment.only=treatment.only,
                                                                         plot=plot, log.y=log.y)
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
setGeneric(name = "setSlop", def = function(object, treatment.only=TRUE, verbose=TRUE) {standardGeneric("setSlop")} )

#' @export
setMethod( f="setSlop", signature="XevaSet",
           definition=function(object, treatment.only=TRUE, verbose=TRUE)
           {
             rtz = list()
             for(model.id in names(object@experiment))
             {
               if(verbose==TRUE){cat(sprintf("Calculating slope for %s\n", model.id))}
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

