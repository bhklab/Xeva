plotWaterFall <- function(x, y, type, color, title, yname, legend.name,
                          show.legend, sort=sort)
{

  dt <- data.frame(x=x, y=y, color=color, stringsAsFactors = FALSE)
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
#' creates waterfall plot for a given drug
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
#' @param object the \code{XevaSet}
#' @param drug name of the drug
#' @param res.measure PDX model response measure
#' @param group.by group response data
#' @param tissue tissue
#' @param model.id which model.id to plot. Default is \code{NULL} will plot all models
#' @param model.type type of model such as mutated or wild type
#' @param type.color a list with colors used for each type
#' @param legend.name name of the legend
#' @param yname name for y axis
#' @param title title of the plot
#' @param sort default TRUE will sort the data
#' @export
#' @import ggplot2
waterfall <- function(object, drug, res.measure, group.by=NULL, tissue=NULL,
                      model.id=NULL, model.type= NULL, type.color="#cc4c02",
                      legend.name=NULL, yname = NULL, title=NULL, sort=TRUE)
{
  if(is.null(yname)){ yname <- res.measure}
  if(is.null(title)){ title <- sprintf("waterfall plot for %s", drug)}

  res <- summarizeResponse(object, response.measure = res.measure,
                           model.id=model.id, batch.id=batch.id,
                           group.by=group.by,
                           tissue=tissue)
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
  rownames(vx) <- as.character(vx$x)

  if(!is.null(model.id))
  {
    vx <- vx[model.id,]
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
  {
    vx$col <- rep(type.color, nrow(vx))[1:nrow(vx)]
  }

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

  plt <- plotWaterFall(x=vx$x, y=vx$y, type=vx$type, color=vx$col, title, yname,
                       legend.name, show.legend, sort=sort)
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


