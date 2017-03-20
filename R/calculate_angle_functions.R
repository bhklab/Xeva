

.computAngle <- function(a,b){
  if(length(a)>4 | length(b)>4)
  {
    msg <- sprintf("Only give start and end points of the lines. For example
                   a = c(start.point.x, start.point.y, end.point.x, end.point.y)")
    stop(msg)
  }
  theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
  as.numeric(theta)*180/base::pi
}


.computSlopFun <- function(x,y)
{
  data = data.frame(x=x, y=y)
  ##---- remove all non finite (Inf, NA, NaN) data --------------
  data = data[is.finite(data$x), ]
  data = data[is.finite(data$y), ]

  p= data$x[1]; q = data$y[1]
  fit <- (lm(I(y-q)~I(x-p) +0, data))
  ang <- atan(coef(fit)[["I(x - p)"]]) *180/base::pi

  return(list(fit=fit, angel=ang, data=data))
}

#' @import ggplot2
plot_Batch_angel_ggplot <- function(dt.fit, dc.fit, liC, liT,
                                    title="plot", log.y=FALSE)
{
  ##-----make one DF ----------------------------------------
  dft = do.call(rbind, lapply(names(dt.fit), function(n)
                       {d=dt.fit[[n]]$data; d$model.id=n; d } ))
  dft$type = "treatment"

  dfc = do.call(rbind, lapply(names(dc.fit), function(n)
                       {d=dc.fit[[n]]$data; d$model.id=n; d } ))
  dfc$type = "control"

  DF = rbind(dft, dfc)

  ##-----normalize data to 0 and 1 --------------------------------
  if(1==2){
  dfXmin <- min(DF$x); dfXmax <- max(DF$x)
  dfYmin <- min(DF$y); dfYmax <- max(DF$y)
  normalize01 <- function(v, vmin, vmax) { (v-vmin)/(vmax-vmin) }

  DF$x <- normalize01(DF$x, dfXmin, dfXmax)
  DF$y <- normalize01(DF$y, dfYmin, dfYmax)

  liC$x1 <- normalize01(liC$x1, dfXmin, dfXmax)
  liC$x2 <- normalize01(liC$x2, dfXmin, dfXmax)
  liC$y1 <- normalize01(liC$y1, dfYmin, dfYmax)
  liC$y2 <- normalize01(liC$y2, dfYmin, dfYmax)

  liT$x1 <- normalize01(liT$x1, dfXmin, dfXmax)
  liT$x2 <- normalize01(liT$x2, dfXmin, dfXmax)
  liT$y1 <- normalize01(liT$y1, dfYmin, dfYmax)
  liT$y2 <- normalize01(liT$y2, dfYmin, dfYmax)
  }
  ##----------------------------------------------------------------
  tcCol <- c("control" = "#6baed6", "treatment" = "#fc8d59")
  plt <- ggplot(DF, aes_string(x="x", y="y", color= "type", group="model.id"))
  plt <- plt + geom_line(linetype = 2)+ geom_point()
  plt <- plt + scale_color_manual(values=tcCol)

  ##-------add lm line -----------------------------------------------------------
  plt <- plt+ geom_segment(aes(x = liC$x1, y = liC$y1, xend = liC$x2, yend = liC$y2), color="blue")
  plt <- plt+ geom_segment(aes(x = liT$x1, y = liT$y1, xend = liT$x2, yend = liT$y2), color="red")

  ##-------------------------------------------------------------------------------
  plt <- plt + labs(title = title, x = "time", y = "volume", colour = "")

  if(log.y==TRUE)
  {
    plt <- plt+scale_y_continuous(trans='log')
    plt <- plt + labs(y = "log(volume)", colour = "")
    #plt <- plt+ coord_fixed(ratio = (max(DF$x)-min(DF$x)) / ( log(max(DF$y))- log(min(DF$y))))

  } else{
    plt <- plt+ coord_fixed(ratio = (max(DF$x)-min(DF$x)) / (max(DF$y)-min(DF$y)))
  }

  #plt <- plt + theme(aspect.ratio=1)
  plt <- .ggplotEmptyTheme(plt)
  plt <- plt + theme(plot.title = element_text(hjust = 0.5))
  plt <- plt + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  return(plt)
}

####-----------------------------------------------------------
#' data(lpdx); object=lpdx
#' expDegI  <- expDesign(lpdx, "PHLC111_P7")
computAngelFor1ExpDesign <- function(object, expDegI, var="volume", treatment.only=TRUE,
                                     plot=FALSE, log.y=FALSE)
{

  ##------ for each model get fit ---------
  dt.fit= list()
  if(length(expDegI$treatment)>0)
  {
    for(ti in expDegI$treatment)
    {
      dt = getExperiment(object, model.id=ti, treatment.only = TRUE)

      dt.fit[[ti]] = .computSlopFun(dt$time, dt[,var])
    }
  }
  ###-------------
  dc.fit= list()
  if(length(expDegI$control)>0)
  {
    for(ci in expDegI$control)
    {
      dc = getExperiment(object, model.id=ci, treatment.only = TRUE)
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

  ##--------------------------------------------------------------------
  .getLMfitLine <- function(fit, data)
  {
    x1 <- data$x[1]
    x2 <- data$x[length(data$x)]
    if(1==1){
    intercept <- predict(fit, newdata = list(x=0))+ data$y[1]
    y1 <- (as.numeric(coef(fit)) * x1) + intercept
    y2 <- (as.numeric(coef(fit)) * x2) + intercept}
    return(list(x1=x1, y1=y1, x2=x2,y2=y2))
  }
  liC <- .getLMfitLine(fit = fitC$fit, data=fitC$data)
  liT <- .getLMfitLine(fit = fitT$fit, data=fitT$data)
  ##--------------------------------------------------------------------

  angDiff <- .computAngle(c(liC$x1, liC$y1, liC$x2, liC$y2), c(liT$x1, liT$y1, liT$x2, liT$y2))

  plt <- NULL
  if(plot==TRUE)
  {
    plt <- plot_Batch_angel_ggplot(dt.fit, dc.fit, liC, liT, title=expDegI$batch.name,
                                   log.y=log.y)
  }
  return(list(angle= angDiff, plot=plt))
  #return(angDiff)
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
#' myDesign = list(batch.name="myBatch", treatment=c("X.015.BY19"), control=c("X.015.uned"))
#' angl = calculateAngle(object=pdxe, myDesign, var = "volume", treatment.only=TRUE, plot=TRUE)
#' print(angl$myBatch$angle)
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
                          plot=TRUE, log.y=FALSE)
                          {standardGeneric("calculateAngle")} )

#' @export
setMethod( f=calculateAngle,
           signature=c(object="XevaSet"),
           definition= function(object, ExpDesign=NULL, var="volume", treatment.only=TRUE,
                                plot=TRUE, log.y=FALSE)
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
  mod = getExperiment(object, model.id=model.id)
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

