##-----
library(BBmisc)
waterfallPlot <- function(vx, col="#cc4c02", title="", yname = "volume")
{
  vx <- data.frame(t(as.matrix(vx)))
  vx$id <- rownames(vx)
  vx$col<- col
  colnames(vx) <- c("value", "id", "col")
  vx <- vx[!is.na(vx$value), ]

  vx <- BBmisc::sortByCol(vx, c("value", "id"), asc = FALSE)
  vx$id <- factor(vx$id, levels = vx$id)

  plt <- ggplot(vx, aes_string("id", "value") ) +
         geom_bar(stat = "identity", aes_string(fill = "col"))

  plt <- plt+scale_fill_manual(values=c(vx$col))

  plt <- plt +theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank())

  ## abb x axis line
  plt <- plt + geom_hline(yintercept=0)

  plt <- .ggplotEmptyTheme(plt)
  plt <- plt + labs(title = title, y = yname, colour = "")
  plt <- plt + theme(plot.title = element_text(hjust = 0.5))

  ry <- c(min(vx$value), max(vx$value)+1)

  #bry<- floor( seq(0, ry[2], length.out = 5) )
  bry<- floor( seq(0, ry[2], 25) )
  bry<- c(-bry, bry)
  bry<- bry[ bry> (ry[1]-25) ]

  #bry<- floor( seq(ry[1], ry[2], length.out = 10) )
  #bry<- sort(unique(c(bry, 0)))
  plt <- plt + scale_y_continuous(breaks=bry, limits = ry)
  plt + theme(legend.position="none")
}

#' Plot drug response waterfall
#'
#' @examples
#' data(pdxe)
#' drugWaterfall(pdxe, drug="binimetinib", value="best.avg.response_published",
#'               col="#E69F00", tumor.type = "CRC", title="Binimetinib",
#'               yname = "Change in tumor volume (%)")
#' @param object The \code{XevaSet}
#' @param drug Name of the drug
#' @param value Which value should be ploted
#'
#' @export
drugWaterfall <- function(object, drug, value, col="#cc4c02", group.by="patient.id",
                          tumor.type=NULL, title="", yname = "Response")
{

  df = summarizeResponse(pdxe, response.measure=value, group.by=group.by,
                         tumor.type=tumor.type)

  if((drug %in% rownames(df))==FALSE)
  {
    msg <- sprintf("Drug %s not present\n", drug)
    stop(msg)
  }
  vx = df[drug,]
  waterfallPlot(vx, col=col, title=title, yname = yname)
}
