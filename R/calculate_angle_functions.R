# This will comput angle from base line aka (0,0) (1,0)
#.computAngle_New(c(0,0), c(1,1))
#.computAngle_New(c(0,0), c(1,-1))
.computAngleOfLine <- function(startCod, endCod)
{
  theta <- atan2(endCod[2], endCod[1]) - atan2(startCod[2], startCod[1])
  thetaD <- as.numeric(theta)*180/base::pi
  return(thetaD)
}

## .computAngle(c(0,0), c(1,0), c(0,0), c(1,1))
## .computAngle(c(0,0), c(1,-1),c(0,0), c(1,1))
## anticlock wise angles are positive and clockwise are -ve
.computAngle <- function(l1StCo, l1EndCo, l2StCo, l2EndCo)
{
  # if(length(a)>4 | length(b)>4)
  # {
  #   msg <- sprintf("Only give start and end points of the lines. For example
  #                  a = c(start.point.x, start.point.y, end.point.x, end.point.y)")
  #   stop(msg)
  # }
  #theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
  #thetaD <- as.numeric(theta)*180/base::pi
  l1Ang <- .computAngleOfLine(l1StCo, l1EndCo)
  l2Ang <- .computAngleOfLine(l2StCo, l2EndCo)
  angDiff <- l2Ang - l1Ang
  return(angDiff)
}

.computSlopeFun <- function(x,y)
{
  df <- data.frame(time=x, volume=y)
  ##---- remove all non finite (Inf, NA, NaN) data --------------
  df <- df[is.finite(df$time), ]
  df <- df[is.finite(df$volume),]

  df$time  <- df$time- df$time[1]
  df$volume<- df$volume- df$volume[1]

  fit <- lm(volume~time +0, df)
  ang <- atan(coef(fit)[["time"]]) *180/base::pi

  #df <- df[is.finite(df$x), ]
  #df <- df[is.finite(df$y), ]
  #p <- df$x[1]; q <- df$y[1]
  #fit <- (lm(I(y-q)~I(x-p) +0, df))
  #ang <- atan(coef(fit)[["I(x - p)"]]) *180/base::pi

  #a=df$time; b=df$volume
  #theta <- acos(sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
  #theta*180/base::pi

  return(list(fit=fit, angel=ang)) # , df=df))
}

####----------------------------------------------------------------------------
#' Computes angle
#' Compute angle for given models or batch
#'
#' @examples
#' data(pdxe)
#' angle(object=pdxe, model.id=c("X.007.BG98", "X.6047.uned"))
#' angle(object=pdxe, batchName="X-6047.paclitaxel")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id} for which slope is required
#' @param batchName batch name
#' @param expDig Experiment design list
#' @param treatment.only default \code{FALSE}. Given full data treatment.only=\code{TRUE} will plot data only during treatment
#' @param max.time maximum time point of the plot, default \code{NULL} will plot complete data
#' @param impute.value default \code{TRUE} will impute missing values
#' @param return.fit default \code{FALSE}. If \code{TRUE} will return \code{lm} fit also
#'
#' @return a \code{data.fram} with treatment, control and batch.name
#' @export
angle <- function(object, model.id=NULL, batchName=NULL, expDig=NULL,
                  treatment.only=TRUE, max.time=NULL, impute.value=TRUE,
                  return.fit=FALSE)
{
  if(!is.null(model.id))
  {
    rtx <- list()
    for(mid in c(model.id))
    {
      d <- getExperiment(object, model.id= mid,treatment.only=treatment.only,
                         max.time = max.time)

      if(nrow(d)<2)
      {
        warning("too few data points to comput slope")
        return(NA)
      }

      ang <- .computSlopeFun(d$time, d$volume.normal)
      rtx[[mid]] <- ang
    }

    if(return.fit==FALSE)
    {
      rty <- sapply(rtx, "[[", "angel")
      return(rty)
    }
    return(rtx)
  }

  if(!is.null(batchName) | !is.null(expDig))
  {
    bt <- getTimeVarData(object, batchName=batchName, expDig=expDig,
                         treatment.only=treatment.only, impute.value=impute.value,
                         max.time=max.time)

    ###----compute angle -----
    rtx <- list()
    rty <- c(NA,NA,NA); names(rty) <- c("control.slope", "treatment.slope", "angle")

    if(sum(bt$exp.type=="control")>1)
    {
      crAng <- .computSlopeFun(bt$time[bt$exp.type=="control"],
                               bt$mean[bt$exp.type=="control"])
      rtx[["control.slope"]] <- crAng
      rty["control.slope"] <- crAng$angel
    }

    if(sum(bt$exp.type=="treatment")>1)
    {
      trAng <- .computSlopeFun(bt$time[bt$exp.type=="treatment"],
                               bt$mean[bt$exp.type=="treatment"])
      rtx[["treatment.slope"]] <- trAng
      rty["treatment.slope"] <- trAng$angel
    }

    if(!is.null(rtx$control.slope) & !is.null(rtx$treatment.slope))
    {
      rtx[["angle"]] <- rtx$control.slope$angel - rtx$treatment.slope$angel
      rty["angle"] <- rtx[["angle"]]
    }

    if(return.fit==FALSE)
    { return(rty) }

    return(rtx)
  }
}

####----------------------------------------------------------------------------
####----------------------------------------------------------------------------
####----------------------------------------------------------------------------
### data(lpdx); object=lpdx
### expDegI  <- expDesign(lpdx, "PHLC111_P7")
.computAngelFor1ExpDesign <- function(object, expDegI, var="volume", treatment.only=TRUE,
                                     plot=FALSE, log.y=FALSE)
{

  ##------ for each model get fit ---------
  dt.fit= list()
  if(length(expDegI$treatment)>0)
  {
    for(ti in expDegI$treatment)
    {
      dt = getExperiment(object, model.id=ti, treatment.only = treatment.only)

      dt.fit[[ti]] = .computSlopeFun(dt$time, dt[,var])
    }
  }
  ###-------------
  dc.fit= list()
  if(length(expDegI$control)>0)
  {
    for(ci in expDegI$control)
    {
      dc = getExperiment(object, model.id=ci, treatment.only = treatment.only)
      dc.fit[[ci]] = .computSlopeFun(dc$time, dc[,var])
    }
  }

  DFx = getTimeVarData(object, expDegI, var=var, treatment.only=treatment.only)
  dfC = DFx[DFx$exp.type=="control",]
  dfT = DFx[DFx$exp.type=="treatment",]

  if(nrow(dfC)==0)
  {
    warning("Control have no data!")
    #return(NA)
    return(list(angle= NA, plot=NULL))
  }

  if(nrow(dfT)==0)
  {
    warning("treatment have no data!")
    #return(NA)
    return(list(angle= NA, plot=NULL))
  }
  fitC = .computSlopeFun(dfC$time, dfC$mean)
  fitT = .computSlopeFun(dfT$time, dfT$mean)

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

  #angDiff <- .computAngle(c(liC$x1, liC$y1, liC$x2, liC$y2), c(liT$x1, liT$y1, liT$x2, liT$y2))
  angDiff <- .computAngle(c(liC$x1, liC$y1), c(liC$x2, liC$y2), c(liT$x1, liT$y1), c(liT$x2, liT$y2))

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
               angDiffX[[expDegI$batch.name]] = .computAngelFor1ExpDesign(object, expDegI, var=var, treatment.only=treatment.only,
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
setGeneric(name = "setAngle", def = function(object, treatment.only=FALSE) {standardGeneric("setAngle")} )

#' @export
setMethod( f=setAngle, signature="XevaSet",
           definition=function(object, treatment.only=FALSE)
           { calculateAngle(object, var = "volume", treatment.only=treatment.only, plot=FALSE) } )


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
             #for(I in 1:length(expDesignInfo(object)))
             #{
              # object@expDesign[[I]]$angle = value[[object@expDesign[[I]]$batch.name]]
             #}

             object@sensitivity$batch[ , "angle"] <- NA
             for(bn in names(value))
             {
               object@sensitivity$batch[bn, "angle"] <- value[[bn]][["angle"]]
             }
             return(object)
           } )


#####===========================================================================
#####=============================== setSlope ==================================

##------- Compute slope for one model ----------------------------
computeSlope <- function(object, model.id, treatment.only=TRUE)
{
  mod <- getExperiment(object, model.id=model.id, treatment.only=treatment.only)
  rtx <- .computSlopeFun(mod$time, mod$volume.normal)
  return(rtx)
}


#' setSlope generic
#'
#' Generic for setSlope method
#'
#' @examples
#' data(pdxe)
#' slope <- setSlope(object = pdxe, treatment.only=FALSE)
#' @param object The \code{XevaSet} to retrieve drug info from
#' @param treatment.only Default \code{TRUE}, take only first periode of treatment where dose is not zero
#' @details If dose column is not present treatment.only will be considered FALSE
#' @return a \code{list} with slope of each model.id
setGeneric(name = "setSlope", def = function(object, treatment.only=TRUE, verbose=TRUE) {standardGeneric("setSlope")} )

#' @export
setMethod( f="setSlope", signature="XevaSet",
           definition=function(object, treatment.only=TRUE, verbose=TRUE)
           {
             rtz = list()
             for(model.id in names(object@experiment))
             {
               if(verbose==TRUE){cat(sprintf("Calculating slope for %s\n", model.id))}
               fa = computeSlope(object, model.id, treatment.only=treatment.only)
               rtz[[model.id]] = fa$angel
             }
             rtz <- unlist(rtz) ##create flate list
             return(rtz)
           })

#' setSlope<- Generic
#'
#' Generic for setSlope replace method
#'
#' @examples
#' data(pdxe)
#' setSlope(pdxe) <- setSlope(object = pdxe)
#' @param object The \code{XevaSet} to replace drug info in
#' @return Updated \code{XevaSet}
setGeneric(name= "setSlope<-", def = function(object, value) {standardGeneric("setSlope<-")} )

#' @export
setMethod( f="setSlope<-",
           signature=c(object = "XevaSet"),
           definition=function(object, value)
           {
             object@sensitivity$model$slope <- NA
             ck1 <- symmetricSetDiff(rownames(object@sensitivity$model), names(value))
             if(length(ck1)!=0)
             {
               stop("Error model.id are not same in sensitivity and experiment slot\n%s", paste(ck1, collapse = "\n"))
             }
             object@sensitivity$model$slope <- unlist(value)[rownames(object@sensitivity$model)]
             return(object)
           })

