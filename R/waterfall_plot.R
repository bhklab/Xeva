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
                          sortByType=FALSE,
                          title="Waterfall plot", yname = "Response",
                          tumor.type=NULL) #, group.by="patient.id")
{

  possibleValues <- colnames(slot(object, "sensitivity")[["model"]])
  if(is.element(value, possibleValues)==FALSE)
  {
    msg <- sprintf("Sensitivity value=%s not present in dataset. Possible sensitivity value are:\n%s",
                   value, paste(possibleValues, collapse=", "))
    stop(msg)
  }

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
    vx$type <- "model"
  }

  colnames(vx) <- c("value", "id", "col", "type")
  vx <- vx[!is.na(vx$value), ]
  if(nrow(vx)==0)
  {stop("No data present for the drug")}

  if(sortByType==TRUE)
  {
    vx <- BBmisc::sortByCol(vx, c("type", "value", "id"), asc = FALSE)
  } else
  {
    vx <- BBmisc::sortByCol(vx, c("value", "id"), asc = FALSE)
  }

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






###--------------------
##--- add oncoplot at the bottam of waterfall plot ---------
.add_OncoplotAt_bottam <- function()
{
  library(ggplot2)
  library(gridExtra)
  pMain <- ggplot(mtcars, aes(x = wt, y = mpg)) +
    geom_point()
  pTop <- ggplot(mtcars, aes(x = wt)) +
    geom_histogram()
  pRight <- ggplot(mtcars, aes(x = mpg)) +
    geom_histogram() + coord_flip()
  pEmpty <- ggplot(mtcars, aes(x = wt, y = mpg)) +
    geom_blank() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          panel.background = element_blank())

  grid.arrange(pMain, pRight,pTop, pEmpty,
               ncol = 2, nrow = 2, widths = c(3, 1), heights = c(3, 1))

}


