plotWaterFall <- function(x, y, type, color, title, yname, legend.name,
                          show.legend, sort=sort)
{

  dt <- data.frame(x=x, y=y, color=color, stringsAsFactors = FALSE)
  dt <- dt[!is.na(dt$y), ]

  dt[, legend.name] <- type

  if(sort==TRUE)
  { dt <- BBmisc::sortByCol(dt, c("y", legend.name, "x"), asc = FALSE) }
  dt$x <- factor(dt$x, levels = as.character(dt$x))

  plt <- ggplot(dt, aes_string(x="x", y="y", fill=legend.name))
  plt <- plt + geom_bar(stat = "identity")

  colValX <- unique(dt[, c("color", legend.name)])
  colVal <- as.character(colValX[, "color"])
  names(colVal)  <- colValX[, legend.name]

  plt <- plt + scale_fill_manual(values=colVal)
  plt <- plt +theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank())
  ##--------- add x axis line ----------------------------
  plt <- plt + geom_hline(yintercept=0, size =0.25)
  plt <- .ggplotEmptyTheme(plt)

  ##----remove x axis ------------------
  plt <- plt +theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(), axis.line.x = element_blank())

  plt <- plt + theme(plot.title = element_text(hjust = 0.5))
  plt <- plt + labs(title = title, y = yname)
  if(show.legend==FALSE)
  {
    plt <- plt + theme(legend.position="none")
  }
  return(plt)
}

#' waterfall plot
#' Creates waterfall plot for a given drug.
#'
#' @examples
#' data(brca)
#' waterfall(brca, drug="binimetinib", res.measure="best.avg.response_published")
#' ## example with model.type where we color the models by TP53 mutation type
#' mut <- summarizeMolecularProfiles(brca,drug = "binimetinib", mDataType="mutation")
#' model.type <- Biobase::exprs(mut)["TP53", ]
#' waterfall(brca, drug="binimetinib", res.measure="best.avg.response_published",
#'           tissue="BRCA", model.id=names(model.type), model.type= model.type)
#'
#' @param object The \code{XevaSet} object
#' @param res.measure PDX model drug response measure
#' @param drug Name of the drug
#' @param group.by Group drug response data
#' @param summary.stat How to summarize multiple values
#' @param tissue Tissue type
#' @param model.id Indicates which \code{model.id} to plot. Default \code{NULL} will plot all models
#' @param model.type Type of model, such as mutated or wild type
#' @param type.color A list with colors used for each type in the legend
#' @param legend.name Name of the legend
#' @param yname Name for the y-axis
#' @param title Title of the plot
#' @param sort Default \code{TRUE} will sort the data
#'
#' @export
#' @import ggplot2
waterfall <- function(object, res.measure, drug=NULL, group.by=NULL,
                      summary.stat = c(";", "mean", "median"),
                      tissue=NULL,
                      model.id=NULL, model.type= NULL, type.color="#cc4c02",
                      legend.name=NULL, yname = NULL, title=NULL, sort=TRUE)
{
  if(is.null(yname)){ yname <- res.measure}

  res <- summarizeResponse(object, response.measure = res.measure,
                           model.id=model.id,
                           group.by=group.by,
                           summary.stat=summary.stat,
                           tissue=tissue)

  if(is.null(drug))
  { drug <- unique(rownames(res)) }

  if(!(drug %in% rownames(res)))
  { stop(sprintf("drug %s not present in dataset (or tissue subset)", drug))}

  vl <- unlist(res[drug, ])
  validVl <- names(vl)[!is.na(vl)]
  vl <- vl[validVl]
  if(length(vl)==0)
  { stop(sprintf("No valid value of %s present in dataset (or tissue subset)",
                 res.measure))}

  if(class(vl)!="numeric")
  {stop(sprintf("%s is not a numeric response\n", res.measure))}

  vx <- data.frame(x=names(vl), y=vl, type=drug, col="#cc4c02",
                   stringsAsFactors = F)
  vx <- vx[!is.na(vx$y), ]; vx <- vx[!is.na(vx$x), ]
  rownames(vx) <- as.character(vx$x)

  if(!is.null(model.id))
  {
    vx <- vx[vx$x%in%model.id, ]
    if(nrow(vx)==0)
    { msg <- sprintf("given model.id are not present in the objcect\n") }
  }

  if(!is.null(model.type))
  {
    if(is.null(model.id))
    { stop("specifying 'model.id' is nesseary for 'model.type'") }
    vx[model.id, "type"] <- model.type
  }

  if(length(unique(vx$type))==1)
  { vx$col <- rep(type.color, nrow(vx))[1:nrow(vx)] }

  if(length(unique(vx$type))>1)
  {
    if(class(type.color)!="list")
    {
      type.color <- as.list(rainbow(length(unique(vx$type))))
      names(type.color) <- unique(vx$type)
    }
    vx$col <- unlist(type.color[vx$type])
  }

  show.legend <- TRUE
  if(is.null(legend.name))
  {
    legend.name <- "type"
    if(length(unique(vx$type))==1)
    {show.legend <- FALSE }
  }

  vx <- vx[!is.na(vx$x), ]; vx <- vx[!is.na(vx$y), ]
  plt <- plotWaterFall(x=vx$x, y=vx$y, type=vx$type, color=vx$col, title, yname,
                       legend.name, show.legend, sort=sort)
  return(plt)
}



##------------------------------------------------------------------------------
###--------------------
##--- add oncoplot at the bottam of waterfall plot ---------
.add_OncoplotAt_bottam <- function()
{

}


