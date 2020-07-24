#' @import ggplot2
.addlinetoplot <- function(dt, x, y, col='red', lty="dotted", alpha=1, size=0.5)
{
  list(
    geom_line( data=dt, aes_string(x=x, y=y), color=col, linetype=lty,
               alpha=alpha, size=size),
    geom_point(data=dt, aes_string(x=x, y=y), color=col, shape=20, alpha=alpha)
  )
}

#' @import ggplot2
.plotModelErrorBar <- function(dfp,control.col="#6baed6",treatment.col="#fc8d59",
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
          plt <- plt + geom_errorbar(aes_string(ymin = "lower", ymax = "upper"),
                                     width=0.25)
        }
        if(SE.plot == "ribbon")
        {
          plt <- plt + geom_ribbon(aes_string(ymin ="lower", ymax ="upper",
                                              fill ="exp.type"),
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

  plt <- plt + geom_point(data=df, aes_string(x="time", y="mean",
                                              color= "exp.type"),
                          shape=21, fill="white")

  tcCol <- c("control" = control.col, "treatment" = treatment.col)
  plt <- plt + scale_color_manual(values=tcCol)
  if(SE.plot == "ribbon")
  { plt <- plt + scale_fill_manual(values=tcCol) }
  plt <- .ggplotEmptyTheme(plt)
  plt <- plt + labs(title = title, x = xlab, y = ylab, colour = drgName,
                    fill=drgName)
  plt <- plt + theme(plot.title = element_text(hjust = 0.5))
  plt <- plt + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
  if(!is.null(aspect.ratio))
  {
    plt <- plt + theme(aspect.ratio=aspect.ratio)
  }
  return(plt)
}

#' @import ggplot2
.plotMultipalModels <- function(dfx, color=NULL, major.line.size=1,
                                aspect.ratio=NULL, model.lyt=1)
{
  plt <- ggplot(dfx, aes_string(x="time", y="volume", color= "model.id"))
  plt <- plt + geom_line(linetype = model.lyt, size=major.line.size)+ geom_point()

  if(!is.null(color))
  { plt <- plt+scale_color_manual(values = color) }

  plt <- .ggplotEmptyTheme(plt)

  if(!is.null(aspect.ratio))
  {
    plt <- plt + theme(aspect.ratio=aspect.ratio)
  }
  return(plt)

}


#' @import ggplot2
.plotDose <- function(do, point.shape=21, point.size=5, point.color="black",
                     line.size=4, line.color="black", modify.x.axis=TRUE)
{
  plt <- ggplot(do, aes_string(x="time", y="model.id"))
  plt <- plt + geom_point(size=0)
  plt <- plt + geom_hline(aes_string(yintercept = "model.n"), do,
                          size=line.size, color=line.color)
  plt <- plt + geom_point(fill=do$color,shape=point.shape, size=point.size,
                          color=point.color)

  if(modify.x.axis==TRUE)
  {
    unqTime <- unique(do$time)
    plt <- plt + scale_x_continuous(breaks=unqTime, labels=unqTime)
  }

  plt <- plt + theme_bw()
  return(plt)
}

#' plot dose data
#'
#' plot data for dose in model.id
#'
#' @param object Xeva object.
#' @param model.id one or multiple model.id
#' @param max.time Maximum time point of the plot. Default \code{NULL} will plot complete data
#' @param treatment.only Default \code{FALSE}. Given full data \code{treatment.only=TRUE} will plot data only during treatment
#' @param vol.normal Default \code{FALSE}. If \code{TRUE}, volume will be normalized
#' @param concurrent.time Default \code{FALSE}. If \code{TRUE}, cut the batch data such that control and treatment will end at the same time point
#' @param point.shape shape of the point
#' @param point.size size of the point
#' @param line.size size of the line
#' @param point.color color for point
#' @param line.color color for line
#' @param fill.col a vector with color to fill
#' @param modify.x.axis Default \code{FALSE}
#'
#' @return A ggplot2 plot
#'
#' @examples
#' data(brca)
#' dosePlot(brca, model.id=c("X.6047.LJ16","X.6047.LJ16.trab"),
#'          fill.col=c("#f5f5f5", "#993404"))
#' @export
dosePlot <- function(object, model.id, max.time=NULL, treatment.only=FALSE,
                     vol.normal=FALSE, concurrent.time=FALSE,
                     point.shape=21, point.size=3, line.size=4,
                     point.color="#878787", line.color="#bababa",
                     fill.col=c("#f5f5f5", "#E55100"),
                     modify.x.axis=FALSE)
{
  dfx <- getExperiment(object, model.id=model.id,
                       treatment.only=treatment.only, max.time=max.time,
                       vol.normal=vol.normal, return.list = FALSE,
                       concurrent.time = concurrent.time)
  if(is.null(dfx$dose))
  {
    warning("no dose information present! assuming dose = 1")
    dfx$dose <- 1
  }

  do <- dfx[, c("model.id", "time", "dose")]
  model.order <- unique(do$model.id)

  do$model.n <- as.numeric(factor(as.character(do$model.id), levels = model.order))
  do$color <- ifelse(do$dose==0, fill.col[1], fill.col[2])

  doplt <- .plotDose(do, point.shape, point.size, point.color, line.size,
                     line.color, modify.x.axis)
  return(doplt)
}




#' Plot batch data
#'
#' Plot data for a batch.id, experiment design or model.id
#'
#' @param object Xeva object.
#' @param batch Batch name or experiment design list.
#' @param patient.id Patient id from the \code{XevaSet}. Default \code{NULL}.
#' @param drug Name of the drug. Default \code{NULL}.
#' @param model.id One or multiple model.id. Default \code{NULL}.
#' @param model.color Color for \code{model.id}. Default \code{NULL}.
#' @param control.name Name of the control sample.
#' @param max.time Maximum time point of the plot. Default \code{NULL} will plot complete data.
#' @param treatment.only Default \code{FALSE}. Given full data \code{treatment.only=TRUE} will plot data only during treatment.
#' @param vol.normal Default \code{FALSE}. If \code{TRUE}, volume will be normalized.
#' @param impute.value Default \code{TRUE} will impute values if missing.
#' @param concurrent.time Default \code{FALSE}. If \code{TRUE}, cut the batch data such that control and treatment will end at the same time point.
#' @param control.col Color for control plots.
#' @param treatment.col Color for treatment plots.
#' @param title Title of the plot.
#' @param xlab Title of the x-axis.
#' @param ylab Title of the y-axis.
#' @param log.y Default \code{FALSE}. If \code{TRUE}, y-axis will be log-transformed.
#' @param SE.plot Plot type. Default \code{"all"} will plot all plots and average curves. Possible values are \code{"all"}, \code{"none"}, \code{"errorbar"}, and \code{"ribbon"}.
#' @param aspect.ratio Default \code{1} will create a plot of equal width and height.
#' @param minor.line.size Line size for minor lines. Default \code{0.5}.
#' @param major.line.size Line size for major lines. Default \code{0.7}.
#' @param model.lyt Line type for models. Default \code{"dotted"}.
#'
#' @return A ggplot2 plot with control and treatment batch data.
#'
#' @examples
#' data(brca)
#' plotPDX(brca, model.id=c("X.6047.LJ16","X.6047.LJ16.trab"))
#'
#' plotPDX(brca, batch="X-1004.BGJ398", vol.normal=TRUE)
#' expDesign <- list(batch.name="myBatch", treatment=c("X.6047.LJ16","X.6047.LJ16.trab"),
#'              control=c("X.6047.uned"))
#' plotBatch(brca, batch=expDesign, vol.normal=TRUE)
#' plotBatch(brca, batch=expDesign, vol.normal=FALSE, SE.plot = "errorbar")
#' @export
plotPDX <- function(object, batch=NULL,
                    patient.id=NULL, drug=NULL, model.id=NULL, model.color=NULL,
                    control.name=NULL,
                    max.time=NULL, treatment.only=FALSE, vol.normal=FALSE,
                    impute.value=TRUE, concurrent.time=FALSE,
                    control.col = "#e41a1c", treatment.col = "#377eb8",
                    title="", xlab = "Time", ylab = "Volume",
                    log.y=FALSE, SE.plot = c("all", "none", "errorbar", "ribbon"),
                    aspect.ratio=c(1, NULL),
                    minor.line.size=0.5, major.line.size=0.7,
                    model.lyt="dotted")
{
  if(!is.null(model.id))
  {
    dfx <- getExperiment(object, model.id=model.id,
                         treatment.only=treatment.only, max.time=max.time,
                         vol.normal=vol.normal, return.list = FALSE,
                         concurrent.time = concurrent.time)

    plt <- .plotMultipalModels(dfx, color=model.color,
                               major.line.size=major.line.size,
                               aspect.ratio=aspect.ratio,
                               model.lyt=model.lyt)
    return(plt)
  } else
  {
    plotBatch(object, batch=batch, patient.id=patient.id, drug=drug,
              control.name=control.name, max.time=max.time,
              treatment.only=treatment.only, vol.normal=vol.normal,
              impute.value=impute.value,
              concurrent.time=concurrent.time,
              control.col = control.col, treatment.col=treatment.col,
              title=title, xlab = xlab, ylab = ylab,
              log.y=log.y,
              SE.plot =SE.plot,
              aspect.ratio=aspect.ratio,
              minor.line.size=minor.line.size, major.line.size=major.line.size,
              model.lyt=model.lyt)
  }

}

#' @rdname plotPDX
#' @export
plotBatch <- function(object, batch=NULL, patient.id=NULL, drug=NULL, control.name=NULL,
                      max.time=NULL, treatment.only=FALSE, vol.normal=FALSE,
                      impute.value=TRUE,
                      concurrent.time=FALSE,
                      control.col = "#6baed6", treatment.col="#fc8d59",
                      title="", xlab = "Time", ylab = "Volume",
                      log.y=FALSE,
                      SE.plot = c("all", "none", "errorbar", "ribbon"),
                      aspect.ratio=c(1, NULL),
                      minor.line.size=0.5, major.line.size=0.7,
                      model.lyt="dotted")
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
                     modelLyt= model.lyt,
                     minor.line.size=minor.line.size,
                     major.line.size=major.line.size)
}
