#' ##df = readRDS("DATA-raw/toPlot_DF.Rda")
#' ##plotModelErrorBar(df)
#' @export
#' @import ggplot2
plotModelErrorBar <- function(df, control.col = "#6baed6", treatment.col="#fc8d59",
                              title="", xlab = "Time", ylab = "Volume", log.y=TRUE,
                              SE.plot = c("errorbar", "ribbon"), aspect.ratio=c(1, NULL))
{
  SE.plot <- SE.plot[1]
  aspect.ratio <- aspect.ratio[1]

  df <- df[!is.na(df$mean), ]

  if(!is.null(df$upper) & !is.null(df$lower))
  {
    if(all(is.na(df$upper))==TRUE){ df$upper=NULL}
    if(all(is.na(df$lower))==TRUE){ df$lower=NULL}
  }

  if(!is.null(df$upper) & !is.null(df$lower))
  #if( all(is.na(df$upper))==FALSE & all(is.na(df$lower))==FALSE)
  {
    df <- df[!is.na(df$upper), ]
    df <- df[!is.na(df$lower), ]
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
  plt <- plt + geom_line(linetype = 1)+ geom_point()

  if(!is.null(df$upper) & !is.null(df$lower))
  {
    if(all(is.na(df$upper))==FALSE & all(is.na(df$lower))==FALSE)
    {
      if(SE.plot == "errorbar")
      {
        plt <- plt + geom_errorbar(aes_string(ymin = "lower", ymax = "upper"))
      }
      if(SE.plot == "ribbon")
      {
        plt <- plt + geom_ribbon(aes_string(ymin ="lower", ymax ="upper"), linetype=0, fill = "grey80")
        plt <- plt + geom_line(aes_string(y = "mean"))+ geom_point()
      }
    }
  }
  tcCol <- c("control" = control.col, "treatment" = treatment.col)
  plt <- plt + scale_color_manual(values=tcCol)
  plt <- .ggplotEmptyTheme(plt)

  drgName <- df[df$exp.type=="treatment", "drug.name"][1]

  plt <- plt + labs(title = title, x = xlab, y = ylab, colour = drgName)
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



#batchName = "PHLC153_P6"
##
#' data(lpdx)
#' plotBatch(lpdx, "PHLC153_P6", treatment.only=FALSE, log.y=TRUE)
#' @export
plotBatch <- function(object, batchName, treatment.only=FALSE, ...)
{
  expDig <- expDesign(object, batchName)
  df <- getTimeVarData(object, ExpDesign = expDig, treatment.only = treatment.only,
                       drug.name = TRUE)
  plotModelErrorBar(df=df, ...)
}







