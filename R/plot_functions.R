.addlinetoplot <- function(dt, x, y, col='red', lty="dotted", alpha=1, size=0.5)
{
  list(
    geom_line( data=dt, aes_string(x=x, y=y), color=col, linetype=lty, alpha=alpha, size=size),
    #geom_point(data=dt, aes_string(x=x, y=y), color=col,size=0.5, shape=21, fill="white")
    geom_point(data=dt, aes_string(x=x, y=y), color=col, shape=20, alpha=alpha)
  )
}

#' @import ggplot2
.plotModelErrorBar <- function(dfp, control.col = "#6baed6", treatment.col="#fc8d59",
                              title="", xlab = "Time", ylab = "Volume",
                              log.y=FALSE, drgName="",
                              SE.plot = c("all","none","errorbar", "ribbon"),
                              modelLyt= "dotted",
                              aspect.ratio=c(1, NULL), minor.line.size=0.5,
                              major.line.size=0.7)
{
  SE.plot <- match.arg(SE.plot)
  aspect.ratio <- aspect.ratio[1]

  df <- dfp$mean
  df <- df[!is.na(df$mean), ]

  if(!is.null(df$upper) & !is.null(df$lower))
  {
    if(all(is.na(df$upper))==TRUE){ df$upper=NULL}
    if(all(is.na(df$lower))==TRUE){ df$lower=NULL}
  }

  if(nrow(df)==0)
  {
    stop("No data left after removing NA")
  }

  if(log.y==TRUE)
  {
    df$mean = log(df$mean)
    if(!is.null(df$upper)) {df$upper <- log(df$upper)}
    if(!is.null(df$upper)) {df$lower <- log(df$lower)}
  }

  plt <- ggplot(df, aes_string(x="time", y="mean", color= "exp.type"))
  plt <- plt + geom_line(linetype = 1, size=major.line.size)+ geom_point()

  if(SE.plot %in% c("errorbar", "ribbon"))
  {
    if(!is.null(df$upper) & !is.null(df$lower))
    {
      if(all(is.na(df$upper))==FALSE & all(is.na(df$lower))==FALSE)
      {
        if(SE.plot == "errorbar")
        {
          plt <- plt + geom_errorbar(aes_string(ymin = "lower", ymax = "upper"), width=0.25)
        }
        if(SE.plot == "ribbon")
        {
          plt <- plt + geom_ribbon(aes_string(ymin ="lower", ymax ="upper", fill ="exp.type"),
                                   linetype=0,  alpha = 0.25)
        }
      }
    }
  }

  if(SE.plot =="all")
  {
    if(!is.null(dfp$control))
    {
      for(ct in dfp$control)
      { plt <- plt + .addlinetoplot(ct, x="time", y="volume", col=control.col,
                                   lty=modelLyt, alpha=0.5, size = minor.line.size) }
    }
    if(!is.null(dfp$treatment))
    {
      for(tr in dfp$treatment)
      { plt <- plt + .addlinetoplot(tr, "time", "volume", col=treatment.col,
                                   lty=modelLyt, alpha=0.5,size = minor.line.size) }
    }
  }

  plt <- plt + geom_point(data=df, aes_string(x="time", y="mean", color= "exp.type"),
                          #size=1,
                          shape=21, fill="white")

  tcCol <- c("control" = control.col, "treatment" = treatment.col)
  plt <- plt + scale_color_manual(values=tcCol)
  if(SE.plot == "ribbon")
  { plt <- plt + scale_fill_manual(values=tcCol) }
  plt <- .ggplotEmptyTheme(plt)
  #drgName <- df[df$exp.type=="treatment", "drug.name"][1]

  plt <- plt + labs(title = title, x = xlab, y = ylab, colour = drgName, fill=drgName)
  plt <- plt + theme(plot.title = element_text(hjust = 0.5))
  plt <- plt + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  if(!is.null(aspect.ratio))
  {
    plt <- plt + theme(aspect.ratio=aspect.ratio)
  }

  ###----- to do ---------
  ###----- add line showing Dose at the bottam ---------
  ###----- look at the packege library(ggExtra) --------
  return(plt)
}

######--------------------------------------------------------------------------
#' @import ggplot2
.plotMultipalModels <- function(dfx, color=NULL, major.line.size=1, aspect.ratio=NULL)
{
  plt <- ggplot(dfx, aes_string(x="time", y="volume", color= "model.id"))
  plt <- plt + geom_line(linetype = 1, size=major.line.size)+ geom_point()

  if(!is.null(color))
  { plt <- plt+scale_color_manual(values = color) }

  plt <- .ggplotEmptyTheme(plt)

  if(!is.null(aspect.ratio))
  {
    plt <- plt + theme(aspect.ratio=aspect.ratio)
  }
  return(plt)

}

######--------------------------------------------------------------------------





#' Plot batch data
#'
#' Plot data for a batch id or experiment design
#'
#' @param object Xeva object
#' @param batch batch name or experiment design list
#' @param max.time maximum time point of the plot, default \code{NULL} will plot complete data
#' @param treatment.only default \code{FALSE}. Given full data treatment.only=\code{TRUE} will plot data only during treatment
#' @param vol.normal default \code{FALSE} . If TRUE volume will ne normalised
#' @param impute.value default \code{TRUE}, will impute values where missing
#' @param concurrent.time default \code{FALSE}. If TRUE will cut the batch data such that control and treatment will end at same time point
#' @param control.col color for control plots
#' @param treatment.col color for treatment plots
#' @param title title of the plot
#' @param xlab title of x axis
#' @param ylab title of y axis
#' @param log.y default \code{FALSE}, if \code{TRUE} y axis will be in log
#' @param drug default \code{NULL} will extract drug name from data
#' @param SE.plot plot type. Default \code{"all"} will plot all plots and average curves. Possible values are \code{"all"}, \code{"none"}, \code{"errorbar"}, \code{"ribbon"}
#' @param aspect.ratio default \code{1} will create equeal width and height plot
#' @param minor.line.size line size for minor lines default \code{0.5}
#' @param major.line.size line size for major lines default \code{0.7}
#'
#' @return A ggplot2 plot with control and treatment
#'
#' @examples
#' data(brca)
#' plotPDX(brca, model.id=c("X.6047.LJ16","X.6047.LJ16.trab"))
#'
#' plotPDX(brca, batch="X-1004.BGJ398", vol.normal=TRUE)
#' expDesign <- list(batch.name="myBatch", treatment=c("X.6047.LJ16","X.6047.LJ16.trab"),
#'              control=c("X.6047.uned"))
#' plotBatch(brca, batch=expDesign, vol.normal=T)
#' plotBatch(brca, batch=expDesign, vol.normal=F, SE.plot = "errorbar")
#' @export
plotPDX <- function(object, batch=NULL,
                    patient.id=NULL, drug=NULL, model.id=NULL, model.color=NULL,
                    control.name=NULL,
                    max.time=NULL, treatment.only=FALSE, vol.normal=FALSE,
                    impute.value=TRUE, concurrent.time=FALSE,
                    control.col = "#e41a1c", treatment.col = "#377eb8",
                    title="", xlab = "Time", ylab = "Volume",
                    log.y=FALSE, #drug.name=NULL,
                    SE.plot = c("all", "none", "errorbar", "ribbon"),
                    aspect.ratio=c(1, NULL),
                    minor.line.size=0.5, major.line.size=0.7)
{
  if(!is.null(model.id))
  {
    dfx <- getExperiment(object, model.id=model.id,
                         treatment.only=treatment.only, max.time=max.time,
                         vol.normal=vol.normal, return.list = FALSE,
                         concurrent.time = concurrent.time)

    .plotMultipalModels(dfx, color=model.color, major.line.size=major.line.size,
                        aspect.ratio=aspect.ratio)

  } else
  {
    plotBatch(object, batch=batch, patient.id=patient.id, drug=drug,
              control.name=control.name, max.time=max.time,
              treatment.only=treatment.only, vol.normal=vol.normal,
              impute.value=impute.value,
              concurrent.time=concurrent.time,
              control.col = control.col, treatment.col=treatment.col,
              title=title, xlab = xlab, ylab = ylab,
              log.y=log.y, #drug.name=drug.name,
              SE.plot =SE.plot,
              aspect.ratio=aspect.ratio,
              minor.line.size=minor.line.size, major.line.size=major.line.size)
  }

}

#' @rdname plotPDX
#### @export
plotBatch <- function(object, batch=NULL, patient.id=NULL, drug=NULL, control.name=NULL,
                      max.time=NULL, treatment.only=FALSE, vol.normal=FALSE,
                      impute.value=TRUE,
                      concurrent.time=FALSE,
                      control.col = "#6baed6", treatment.col="#fc8d59",
                      title="", xlab = "Time", ylab = "Volume",
                      log.y=FALSE,
                      SE.plot = c("all", "none", "errorbar", "ribbon"),
                      aspect.ratio=c(1, NULL),
                      minor.line.size=0.5, major.line.size=0.7)
{
  SE.plot <- match.arg(SE.plot)
  aspect.ratio <- aspect.ratio[1]

  dfp <- getExperiment(object, batch=batch,
                       patient.id=patient.id, drug=drug, control.name=control.name,
                       treatment.only=treatment.only, max.time=max.time,
                       vol.normal=vol.normal, return.list = TRUE,
                       impute.value=impute.value,
                       concurrent.time = concurrent.time)

  dfp$mean <- dfp$batch ##plot function uses mean as variable

  if(is.null(drug))
  {drug <- dfp$mean[dfp$mean$exp.type=="treatment", "drug.name"][1] }

  .plotModelErrorBar(dfp, control.col=control.col, treatment.col=treatment.col,
                     title=title, xlab = xlab, ylab = ylab,
                     log.y=log.y, drgName=drug, #.name,
                     SE.plot = SE.plot, aspect.ratio=aspect.ratio,
                     minor.line.size=minor.line.size,
                     major.line.size=major.line.size)
}
