#' @import ggplot2
waterfall.plot <- function(dt, x, y, type, type.color, sort=TRUE, na.value="#878787")
{
  if(sort==TRUE)
  { dt <- BBmisc::sortByCol(dt, c(y, x), asc = c(FALSE, FALSE)) }
  dt[, x] <- factor(dt[, x], levels = as.character(dt[, x]))
  
  plt <- ggplot(dt, aes_string(x=x, y=y, fill=type))
  plt <- plt + geom_bar(stat = "identity")
  plt <- plt + scale_fill_manual(values=type.color, na.value=na.value)
  
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
  return(plt)
}

getColPal <- function(N)
{
  pl <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f')
  if(N > length(pl)) 
  { pl <- rainbow(N, alpha=0.9) }
  return(pl[1:N])
}


getDataForWaterFall <- function(object, res.measure=NULL, drug=NULL, group.by=NULL,
                                summary.stat = c(";", "mean", "median"),
                                type=NULL, type.color=NULL)
{
  res <- summarizeResponse(object, response.measure = res.measure,
                           group.by=group.by,
                           summary.stat=summary.stat,
                           other.col = type,
                           return.type = "data.frame")
  if(!(drug %in% res$drug))
  { stop(sprintf("drug %s not present in dataset (or tissue subset)", drug))}
  
  dt <- res[!is.na(res$drug) & res$drug == drug, ]
  
  x <- ifelse("model.id" %in% colnames(dt), "model.id", "batch.name")
  y <- res.measure
  
  if(is.null(type))
  {
    type <- "type"
    dt[, "type"] <- res.measure
  }
  
  if(is.null(type.color)) 
  { 
    type.color <- getColPal(length(unique(dt[, type]))) 
  }
  
  if(is.null(names(type.color)))
  { names(type.color) <-  unique(dt[, type]) }
  dl <- list(dt=dt, x=x, y=y, type=type, type.color=type.color)
  return(dl)
}
##---------------------------------------

#' waterfall plot
#' Creates waterfall plot for a given drug.
#'
#' @examples
#' data(brca)
#' ##plot best.average.response for all models tested with binimetinib
#' waterfall(brca, drug="binimetinib", res.measure="best.average.response")
#' 
#' ##plot same by taking mean of multiple models of each patient
#' waterfall(brca, drug="binimetinib", res.measure="best.average.response", 
#'           group.by = "patient.id", summary.stat = "mean")
#'
#' ## plot by specifing color by mutation type
#' ## extract mutation information
#' mut <- summarizeData(brca,drug = "binimetinib", mDataType="mutation")
#' model.type <- Biobase::exprs(mut)["TP53", ]
#' ## extract data.frame of response
#' df <- summarizeResponse(brca, response.measure = "best.average.response", return.type = "data.frame")
#' df <- df[df$drug == "binimetinib", ]
#' ## add values to data.frame 
#' df$TP53 <- model.type[df$model.id]
#' ## now plot the data 
#' waterfall(df, x="model.id", y="best.average.response", type="TP53")
#'
#' @param object The \code{XevaSet} object or a data.frame. See \code{Details}.
#' @param res.measure PDX model drug response measure
#' @param drug Name of the drug
#' @param group.by Group drug response data
#' @param summary.stat How to summarize multiple values. Options are ";", "mean" or "median". See \code{Details}.
#' @param x,y If object is data.frame, x and y indicates column names of x and y axis
#' @param type Type for each bar in waterfall (such as mutated or wild type). See Details.
#' @param type.color A color vector for type
#' @param yname Name for the y-axis
#' @param title Title of the plot
#' @param sort Default \code{TRUE} will sort the data
#' @param na.value Color for NA values. Default "#878787"
#'
#' @return waterfall plot in ggplot2
#' 
#' @details The function waterfall can plot from a XevaSet or from a data.frame . 
#' If a data.frame is specified, x and y parameters must also be specified. 
#'
#' @export
#' @import ggplot2
waterfall <- function(object, res.measure=NULL, drug=NULL, group.by=NULL,
                      summary.stat = c(";", "mean", "median"),
                      x = NULL, y = NULL, 
                      type=NULL, type.color=NULL, na.value="#878787",
                      yname = NULL, title=NULL, sort=TRUE)
{
  summary.stat <- match.arg(summary.stat)
 
  if(is(object, "XevaSet"))
  {
    dl <- getDataForWaterFall(object, res.measure, drug, group.by,
                                    summary.stat, type, type.color)
    dt <- dl$dt 
    x <- dl$x
    y <- dl$y
    type <- dl$type
    type.color <- dl$type.color
  }
  
  if(is(object, "data.frame"))
  {
    dt <- object
    if(is.null(x)) { stop("for data.frame x value must be defined")}
    if(is.null(y)) { stop("for data.frame y value must be defined")}
    
    if(is.null(type))
    {
      type <- "type"
      dt[, "type"] <- y
    }
    
    if(is.null(type.color)) 
    { type.color <- getColPal(length(unique(dt[, type]))) }
    if(is.null(names(type.color)))
    { names(type.color) <-  unique(dt[, type]) }
  }
  
  checkCol <- setdiff(c(x, y, type), colnames(dt))
  if(length(checkCol)>0)
  { stop(sprintf("%s is not present in data\n", checkCol)) }
  
  if(length(type.color) < length(unique(dt[, type])))
  { stop(sprintf("type.color should be length %d", length(unique(dt[, type])))) }
  
  plt <- waterfall.plot(dt, x, y, type, type.color, sort, na.value)
  if(length(type.color)==1)
  { plt <- plt + theme(legend.position="none") }
  if(!is.null(title)) { plt <- plt + labs(title = title)}
  if(!is.null(title)) { plt <- plt + labs(y = yname) }
  plt
}

