.summarizePerModelResponse <- function(object, response.measure, model.id,
                                       group.by, summary.stat, tissue)
{
  if(is.element(response.measure, colnames(slot(object, "sensitivity")[["model"]]) )==FALSE)
  { stop(sprintf("'%s' is not present in sensitivity slot\n", response.measure)) }

  dfVal <- slot(object, "sensitivity")[["model"]] [,c("model.id", response.measure)]
  df <- modelInfo(object)
  df[, response.measure] <- dfVal[df$model.id, response.measure]

  if(!is.null(model.id))
  {
    df <- df[model.id, ]
    if(nrow(df)==0)
    { stop(sprintf("given model.id not present in the Xeva object\n")) }
  }

  if(!is.null(tissue))
  {
    df <- df[df$tissue==tissue,]
    if(nrow(df)==0)
    { stop(sprintf("given tissue not present in the Xeva object\n")) }
  }

  if(is.null(group.by)){group.by <- "model.id"}

  if(is.element(group.by, colnames(df))==FALSE)
  { stop(sprintf("'group.by' %s not present in model\n", group.by)) }

  mat <- .castDataFram(df, row.var="drug", col.var = group.by,
                       value=response.measure, collapse = summary.stat)
  return(mat)
}


.summarizePerBatchResponse <- function(object, response.measure = NULL, batch.name=NULL)
{
  rtx <- slot(object, "sensitivity")[["batch"]]
  if(!is.null(response.measure))
  { rtx <- rtx[, c("batch.name", response.measure)] }

  if(!is.null(batch.name))
  {
    bn2take <- batch.name[batch.name %in% rtx$batch.name]
    if(length(bn2take)==0)
    {
      msg <- sprintf("No batch.name present in dataset. Please check the batch.name")
      stop(msg)
    }
    rtx <- rtx[rtx$batch.name %in% bn2take, ]
  }
  return(rtx)
}


.checkResMes <- function(object, response.measure)
{
  rm.type <- NULL
  if(response.measure %in% colnames(slot(object, "sensitivity")[["model"]]))
  { rm.type <- "model" }

  if(response.measure %in% colnames(slot(object, "sensitivity")[["batch"]]))
  { rm.type <- "batch" }

  if(is.null(rm.type))
  {
    msg <- sprintf("valid response.measure values are\nFor model: %s\n\nFor batch: %s\n",
                   paste0(colnames(slot(object, "sensitivity")[["model"]]), collapse = ", "),
                   paste0(colnames(slot(object, "sensitivity")[["batch"]]), collapse = ", ")
                   )
    stop(msg)
  }

  return(rm.type)
}
#####================= summarizeResponse ==================
#' Summarize Response of PDXs
#'
#' Summarize Response of PDXs.
#'
#' @param object The \code{XevaSet}
#' @param response.measure \code{character} . Which response measure to use? Use the responseMeasures function to find out what measures are available for each Xeva set.
#' @param group.by default \code{patient.id}. How the models should be grouped togather. See details
#' @param summary.stat which summary method to use if multipal ids were found
#' @param batch.name a vector of batch names. Default NULL will return all batchs
#' @return a \code{matrix} with rows as drug names, coulmn as \code{group.by} and each cell contains \code{response.measure} for the pair.
#'
#' @details
#' There can be two types of response measure
#' \itemize{
#' \item{per model response : One response value for each Model, e.g. mRECIST_recomputed for each model}
#' \item{per batch response : One response value for each Batch, e.g. angle between treatment and control groups}
#' }
#' In case of \code{per model response} output columns will be \code{model.id} (or group.by).
#' For \code{per batch response} \code{group.by} value can be \code{"batch.name"} .
#'
#' @examples
#' data(brca)
#' brca.mR <- summarizeResponse(brca, response.measure = "mRECIST", group.by="patient.id")
#' @export
summarizeResponse <- function(object, response.measure = "mRECIST",
                              model.id=NULL, batch.id=NULL,
                              group.by="patient.id",
                              summary.stat=c(";", "mean", "median"),
                              tissue=NULL)
{
  summary.stat <- c(summary.stat)[1]

  rm.type <- .checkResMes(object, response.measure)

  if(rm.type=="model")
  {
    mat <- .summarizePerModelResponse(object, response.measure=response.measure,
                                      model.id=model.id, group.by=group.by,
                                      summary.stat=summary.stat,
                                      tissue=tissue)
    return(mat)
  }

  if(rm.type=="batch")
  {
    mat <- .summarizePerBatchResponse(object, response.measure = response.measure,
                                      batch.name=batch.id)
    return(mat)
  }
}






