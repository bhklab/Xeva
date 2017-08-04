.addlinetoplot <- function(dt, x, y, col='red', lty="dotted", alpha=1, size=0.5)
{
  list(
    geom_line( data=dt, aes_string(x=x, y=y), color=col, linetype=lty, alpha=alpha, size=size),
    #geom_point(data=dt, aes_string(x=x, y=y), color=col,size=0.5, shape=21, fill="white")
    geom_point(data=dt, aes_string(x=x, y=y), color=col, shape=20, alpha=alpha)
  )
}

#' ##df = readRDS("DATA-raw/toPlot_DF.Rda")
#' ##plotModelErrorBar(df)
#' @export
#' @import ggplot2
plotModelErrorBar <- function(dfp, control.col = "#6baed6", treatment.col="#fc8d59",
                              title="", xlab = "Time", ylab = "Volume",
                              log.y=FALSE, drgName="",
                              SE.plot = c("all","errorbar", "ribbon"),
                              modelLyt= "dotted",
                              aspect.ratio=c(1, NULL), minor.line.size=0.5,
                              major.line.size=0.7)
{
  SE.plot <- SE.plot[1]
  aspect.ratio <- aspect.ratio[1]

  df <- dfp$mean
  df <- df[!is.na(df$mean), ]

  if(!is.null(df$upper) & !is.null(df$lower))
  {
    if(all(is.na(df$upper))==TRUE){ df$upper=NULL}
    if(all(is.na(df$lower))==TRUE){ df$lower=NULL}

    #df <- df[!is.na(df$upper), ]
    #df <- df[!is.na(df$lower), ]
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


#' #' @export
#' #' @import ggplot2
#' plotDottedPDXCurves <- function(dfp, control.col = "#6baed6", treatment.col="#fc8d59",
#'                                 title="", xlab = "Time", ylab = "Volume",
#'                                 modelLyt= "dotted", meanLty="solid",
#'                                 log.y=FALSE,
#'                                 SE.plot = c("errorbar", "ribbon"), aspect.ratio=c(1, NULL))
#' {
#'
#'   plt <- ggplot()
#'   plt <- .ggplotEmptyTheme(plt)
#'   ##---add mean lines ----------------------------------------------------------
#'   plt <- plt + geom_line( data=dfp$mean, aes_string(x="time", y="mean", color= "type"))
#'   plt <- plt + geom_point(data=dfp$mean, aes_string(x="time", y="mean", color= "type"),
#'                           size=4, shape=21, fill="white")
#'
#'   tcCol <- c("control" = control.col, "treatment" = treatment.col)
#'   plt <- plt + scale_color_manual(values=tcCol)
#'
#'   ##-------------------------------------------------------------
#'   if(SE.plot == "errorbar")
#'   {
#'     #if(!is.null(dfp$controlMean))
#'     #{
#'       plt <- plt + geom_errorbar(data=dfp$mean,
#'                                  aes_string(x="time", ymin = "lower", ymax = "upper"),
#'                                  color=control.col, width=0.25)
#'
#'       plt <- plt + addlinetoplot(dfp$controlMean, x="time", y="mean",
#'                                  col=control.col, lty=meanLty)
#'     #}
#'     if(!is.null(dfp$treatmentMean))
#'     {
#'       plt <- plt + geom_errorbar(data=dfp$treatmentMean,
#'                                  aes_string(x="time", ymin = "lower", ymax = "upper"),
#'                                  color=treatment.col, width=0.25)
#'       plt <- plt + addlinetoplot(dfp$treatmentMean, x="time", y="mean",
#'                                  col=treatment.col, lty=meanLty)
#'     }
#'   }
#'
#'   ###---------------------------------------------------------------------------
#'   if(SE.plot =="ribbon")
#'   {
#'     if(!is.null(dfp$controlMean))
#'     {
#'       plt <- plt + geom_ribbon(data=dfp$controlMean,
#'                                aes_string(x="time", ymin ="lower", ymax ="upper"),
#'                                linetype=0, fill = "grey80", alpha = 0.4)
#'       plt <- plt + addlinetoplot(dfp$controlMean, x="time", y="mean",
#'                                  col=control.col, lty=meanLty)
#'     }
#'
#'     if(!is.null(dfp$treatmentMean))
#'     {
#'       plt <- plt + geom_ribbon(data=dfp$treatmentMean,
#'                                aes_string(x="time", ymin ="lower", ymax ="upper"),
#'                                linetype=0, fill = "grey80", alpha = 0.4)
#'       plt <- plt + addlinetoplot(dfp$treatmentMean, x="time", y="mean",
#'                                  col=treatment.col, lty=meanLty)
#'     }
#'   }
#'
#'   if(SE.plot =="all")
#'   {
#'     if(!is.null(dfp$control))
#'     {
#'       for(ct in dfp$control)
#'       { plt <- plt + addlinetoplot(ct, x="time", y="volume", col=control.col,
#'                                    lty=modelLyt, alpha=0.75) }
#'
#'       plt <- plt + addlinetoplot(dfp$controlMean, x="time", y="mean",
#'                                  col=control.col, lty=meanLty)
#'
#'     }
#'
#'     if(!is.null(dfp$treatment))
#'     {
#'       for(tr in dfp$treatment)
#'       { plt <- plt + addlinetoplot(tr, "time", "volume", col=treatment.col,
#'                                    lty=modelLyt, alpha=0.75) }
#'
#'       plt <- plt + addlinetoplot(dfp$treatmentMean, x="time", y="mean",
#'                                  col=treatment.col, lty=meanLty)
#'     }
#'   }
#'
#'
#'
#'
#'
#'
#'
#'   plt <- plt + labs(title = title, x = xlab, y = ylab)
#'   plt <- plt + theme(plot.title = element_text(hjust = 0.5))
#'   plt <- plt + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
#'   if(!is.null(aspect.ratio))
#'   {
#'     plt <- plt + theme(aspect.ratio=aspect.ratio)
#'   }
#'
#'
#' }


#batchName = "PHLC153_P6"
##
#' data(lpdx)
#' plotBatch(lpdx, "PHLC153_P6", treatment.only=FALSE, log.y=TRUE)
#' @export
plotBatch <- function(object, batchName=NULL, expDig =NULL, treatment.only=FALSE,
                      control.col = "#6baed6", treatment.col="#fc8d59",
                      title="", xlab = "Time", ylab = "Volume",
                      log.y=FALSE, drgName=NULL,
                      SE.plot = c("all","errorbar", "ribbon"),
                      aspect.ratio=c(1, NULL),
                      minor.line.size=0.5, major.line.size=0.7,
                      max.time=NULL)
{
  if(is.null(batchName) & is.null(expDig))
  {
    stop("please provide 'batchName' or 'expDig'")
  }

  if(!is.null(batchName))
  {
    expDig <- expDesign(object, batchName)
  }

  dfp <- list()
  if(!is.null(expDig$control) & length(expDig$control)>0)
  {
    dfp$control <- .getExperimentMultipalIDs(object, mids=expDig$control,
                                             treatment.only=treatment.only)

    if(!is.null(max.time))
    {
      dfp$control <- lapply(dfp$control, function(mi)
                            { mi[mi$time <= max.time , ] })
    }

    #controlMean <- .collapseRplicate(dfp$control, var = "volume")
    #controlMean$type <- "control"
  }

  if(!is.null(expDig$treatment) & length(expDig$treatment)>0)
  {
    dfp$treatment <- .getExperimentMultipalIDs(object, mids=expDig$treatment,
                                             treatment.only=treatment.only)

    if(!is.null(max.time))
    {
        dfp$treatment <- lapply(dfp$treatment, function(mi)
                                { mi[mi$time <= max.time , ] })
    }

    #treatmentMean <- .collapseRplicate(dfp$treatment, var = "volume")
    #treatmentMean$type <-  "treatment"
  }

  dfp$mean <- getTimeVarData(object, ExpDesign = expDig, treatment.only = treatment.only,
                             drug.name = TRUE)

  if(!is.null(max.time))
  {
    dfp$mean <- dfp$mean[dfp$mean$time<= max.time ,]
  }

  if(is.null(drgName))
  {drgName <- dfp$mean[dfp$mean$exp.type=="treatment", "drug.name"][1] }

  plotModelErrorBar(dfp, control.col=control.col, treatment.col=treatment.col,
                    title=title, xlab = xlab, ylab = ylab,
                    log.y=log.y, drgName=drgName,
                    SE.plot = SE.plot, aspect.ratio=aspect.ratio,
                    minor.line.size=minor.line.size,
                    major.line.size=major.line.size)
}








