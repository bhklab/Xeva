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


.summarizePerBatchResponse <- function(object, response.measure = NULL,
                                       batch.name=NULL, summarize=FALSE)
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

  if(summarize==TRUE)
  {
    if(is.null(response.measure)){stop("specify batch response measure")}

    binf <- batchInfo(object, batch = rtx[,"batch.name"],patient.id = TRUE,
                      drug = TRUE)
    rtx$patient.id <- rtx$drug <- NA
    for(i in 1:nrow(rtx))
    {
      b <- binf[[rtx[i,"batch.name"]]]
      if(!is.null(b[["patient.id"]]))
      {
        rtx$patient.id[i] <- paste0(b[["patient.id"]], collapse = ";")
      } else {rtx$patient.id[i] <- NA }

      if(!is.null(b[["drug"]]))
      {
        rtx$drug[i] <- paste0(b[["drug"]], collapse = ";")
      } else { rtx$drug[i] <- NA }
    }

    rtx <- rtx[!is.na(rtx$drug), ]
    rtx <- rtx[!is.na(rtx$patient.id), ]

    meltDF <- matrix(data = NA, nrow = length(unique(rtx$drug)),
                     ncol = length(unique(rtx$patient.id)),
                     dimnames = list(unique(rtx$drug), unique(rtx$patient.id)))

    for(i in 1:nrow(rtx))
    { meltDF[rtx$drug[i], rtx$patient.id[i]] <- rtx[i, response.measure] }

    rtx <- meltDF[, !apply(meltDF, 2, function(i)all(is.na(i)))]
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
#' This function summarizes the drug response information of PDXs.
#'
#' @param object The \code{XevaSet} object.
#' @param response.measure \code{character} indicating which response measure to use. Use the \code{responseMeasures} function to find out what measures are available for each \code{XevaSet}.
#' @param model.id The \code{model.id} for which data is required.
#' @param batch.id A \code{vector} of batch names. Default \code{NULL} will return all batches.
#' @param group.by Default \code{patient.id}. Dictates how the models should be grouped together. See details below.
#' @param summary.stat Dictates which summary method to use if multiple IDs are found.
#' @param tissue Name of the tissue. Default \code{NULL}
#' @param summarize.batch If TRUE, batch response will be summarized as drug patient matrix. Default \code{FALSE}
#'
#' @return A \code{matrix} with rows as drug names, column as \code{group.by}. Each cell contains \code{response.measure} for the pair.
#'
#' @details
#' There can be two types of drug response measure.
#' \itemize{
#' \item{Per model response: One response value for each Model, eg. \code{mRECIST_recomputed} for each model.}
#' \item{Per batch response: One response value for each Batch, eg. \code{angle} between treatment and control groups.}
#' }
#' For the \code{per model response} output, columns will be \code{model.id} (or \code{group.by}).
#' For the \code{per batch response} output, the \code{group.by} value can be \code{"batch.name"}.
#' \code{summarize.batch = TRUE} will return response as drug patient matrix.
#'
#' @examples
#' data(brca)
#' brca.mR <- summarizeResponse(brca, response.measure = "mRECIST", group.by="patient.id")
#'
#' @export
summarizeResponse <- function(object, response.measure = "mRECIST",
                              model.id=NULL, batch.id=NULL,
                              group.by="patient.id",
                              summary.stat=c(";", "mean", "median"), tissue=NULL,
                              summarize.batch=FALSE)
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
                                      batch.name=batch.id,
                                      summarize=summarize.batch)
    return(mat)
  }
}
