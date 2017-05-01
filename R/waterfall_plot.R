#' Plot drug response waterfall
#'
#' @examples
#' data(pdxe)
#' drugWaterfall(pdxe, drug="binimetinib", value="best.avg.response_published",
#'               title="Binimetinib", yname = "Change in tumor volume (%)",
#'               tumor.type = "CRC")
#' @param object the \code{XevaSet}
#' @param drug name of the drug
#' @param value which value should be ploted
#' @param model.id which model.ids to plot. Default is \code{NULL} will plot all models
#' @param model.col color for model.ids. A data.frame with two colomns "color" and "type" can be provided
#' @export
#' @import ggplot2
drugWaterfall <- function(object, drug, value, model.ids=NULL, model.col="#cc4c02",
                          title="Waterfall plot", yname = "Response",
                          tumor.type=NULL) #, group.by="patient.id")
{
  if(is.null(tumor.type)){ warning("might take long time as tumor.type is NULL")}
  df = summarizeResponse(pdxe, response.measure=value, group.by="model.id",
                         tumor.type=tumor.type)

  if((drug %in% rownames(df))==FALSE)
  {
    msg <- sprintf("Drug %s not present\n", drug)
    stop(msg)
  }

  vl <- df[drug,]
  if(!is.null(model.ids)){ vl <- vl[, model.ids]}
  #plt <- .waterfallPlot(vx, model.col=model.col, title=title, yname = yname)

  ##----------------------------------------------------------------------------
  vx <- data.frame(t(as.matrix(vl)))
  vx$id <- rownames(vx)

  if( class(model.col) == "data.frame")
  {
    if(is.element("color", colnames(model.col))==FALSE)
    {
      msg <- sprintf("Column 'color' is missing in 'model.col'")
      stop(msg)
    }
    if(is.element("type" , colnames(model.col))==FALSE )
    {
      msg <- sprintf("Column 'type' is missing in 'model.col'")
      stop(msg)
    }
    vx$col <- model.col[rownames(vx), "color"]
    vx$type <- model.col[rownames(vx), "type"]
  }else
  {
    vx$col<- model.col
    vx$legends <- "model"
    #if(length(model.col)>1)
    #{ vx$col<- model.col[vx$id] } else
    #{ vx$col<- model.col }
  }

  colnames(vx) <- c("value", "id", "col", "type")
  vx <- vx[!is.na(vx$value), ]
  if(nrow(vx)==0)
  {stop("No data present for the drug")}

  vx <- BBmisc::sortByCol(vx, c("value", "id"), asc = FALSE)
  vx$id <- factor(vx$id, levels = vx$id)

  #plt <- ggplot(vx, aes_string(x="id", y= "value", fill="id") ) +
  #  geom_bar(stat = "identity")
  #plt <- plt + scale_fill_manual(values=c(vx$col))

  plt <- ggplot(vx, aes_string(x="id", y= "value", fill="type") ) +
    geom_bar(stat = "identity")

  colValX <- unique(vx[, c("col", "type")])
  colVal <- as.character(colValX[, "col"])
  names(colVal)  <- colValX$type
  plt <- plt + scale_fill_manual(values=colVal)

  plt <- plt +theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank())

  ##--------- add x axis line ----------------------------
  plt <- plt + geom_hline(yintercept=0, size =0.25)
  plt <- .ggplotEmptyTheme(plt)
  plt <- plt + labs(title = title, y = yname, colour = "")
  plt <- plt + theme(plot.title = element_text(hjust = 0.5))
  #plt <- plt + theme(legend.position="none")
  return(plt)
}
