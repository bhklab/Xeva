.summarizePerModelResponse <-
  function(object,
           response.measure,
           model.id,
           group.by,
           summary.stat,
           other.col = NULL,
           tissue = NULL)
  {
    if (is.element(response.measure, colnames(slot(object, "sensitivity")[["model"]])) ==
        FALSE)
    {
      stop(sprintf("'%s' is not present in sensitivity slot\n", response.measure))
    }
    
    dfVal <-
      slot(object, "sensitivity")[["model"]] [, c("model.id", response.measure)]
    df <- modelInfo(object)
    df[, response.measure] <- dfVal[df$model.id, response.measure]
    
    if (!is.null(model.id))
    {
      df <- df[model.id, ]
      if (nrow(df) == 0)
      {
        stop(sprintf("given model.id not present in the Xeva object\n"))
      }
    }
    
    if (!is.null(tissue))
    {
      df <- df[df$tissue == tissue,]
      if (nrow(df) == 0)
      {
        stop(sprintf("given tissue not present in the Xeva object\n"))
      }
    }
    
    if (is.null(group.by)) {
      group.by <- "model.id"
    }
    
    if (is.element(group.by, colnames(df)) == FALSE)
    {
      stop(sprintf("'group.by' %s not present in model\n", group.by))
    }
    
    ##------
    dfcols <-
      unique(c("model.id", response.measure, "drug", "patient.id", group.by))
    if (!is.null(other.col))
    {
      dfcols <- unique(c(dfcols, other.col))
    }
    
    df <- df[, dfcols]
    
    mat <- .castDataFram(
      df,
      row.var = "drug",
      col.var = group.by,
      value = response.measure,
      collapse = summary.stat
    )
    mat <- as.matrix(mat)
    
    ###---------------
    if (summary.stat %in% c("mean", "median"))
    {
      df2 <-
        reshape2::melt(mat,
                       varnames = c("drug", group.by),
                       value.name = response.measure)
      df2[] <-
        lapply(df2, function(x)
          if (is.factor(x))
            as.character(x)
          else
            x)
      df2 <- df2[!is.na(df2[, response.measure]),]
      
      otherCols <- setdiff(colnames(df), colnames(df2))
      if (length(otherCols) > 0)
      {
        df2[, otherCols] <- NA
        for (i in 1:nrow(df2))
        {
          dr <- df2[i, "drug"]
          pt <- df2[i, group.by]
          dv <-
            df[!is.na(df$drug) &
                 df$drug == dr & df[, group.by] == pt & !is.na(df[, group.by]),]
          if (nrow(dv) > 0)
          {
            for (ci in otherCols)
            {
              df2[i, ci] <- paste0(unique(dv[, ci]), collapse = ";")
            }
          }
        }
      }
      df2 <- df2[, colnames(df)]
      df <- df2
    }
    return(list(mat = as.matrix(mat), df = df))
  }

.ModelResponsePerBatch <-
  function(object, response.measure, batch.id)
  {
    batch.name <- batchInfo(object)
    if (!batch.id %in% batch.name) {
      stop("batch.id not present, See batchInfo(object)")
    }
    
    batch.design <- batchInfo(object, batch = batch.id)
    
    rtx <- data.frame()
    if (length(batch.design[[1]]$control) > 0)
    {
      rcn <-
        sensitivity(object, type = "model")[batch.design[[1]]$control, ]
      rcn$type <- "control"
      rtx <- rbind(rtx, rcn)
    }
    
    if (length(batch.design[[1]]$treatment) > 0)
    {
      rtn <-
        sensitivity(object, type = "model")[batch.design[[1]]$treatment, ]
      rtn$type <- "treatment"
      rtx <- rbind(rtx, rtn)
    }
    if (nrow(rtx) > 0)
    {
      rtx$batch.name <- batch.id
      rtx <-
        rtx[, c("model.id", "batch.name", "type", response.measure)]
    }
    
    return(rtx)
  }

.summarizePerBatchResponse <- function(object,
                                       response.measure,
                                       other.col = c("patient.id", "drug"),
                                       group.by = "patient.id",
                                       batch.name = NULL)
{
  rtx <- slot(object, "sensitivity")[["batch"]]
  rtx <- rtx[, c("batch.name", response.measure)]
  
  if (!is.null(batch.name))
  {
    bn2take <- batch.name[batch.name %in% rtx$batch.name]
    if (length(bn2take) == 0)
    {
      msg <-
        sprintf("No batch.name present in dataset. Please check the batch.name")
      stop(msg)
    }
    rtx <- rtx[rtx$batch.name %in% bn2take,]
  }
  
  other.col <- unique(c("patient.id", "drug", other.col, group.by))
  binf <- batchInfo(object,
                    batch = rtx[, "batch.name"],
                    other.col = other.col,
                    return.df = TRUE)
  cols2add <- setdiff(colnames(rtx), "batch.name")
  binf[, cols2add] <- NA
  for (i in 1:nrow(binf))
  {
    binf[i, cols2add] <-
      rtx[rtx$batch.name == binf$batch.name[i], cols2add]
  }
  
  rt <- list(df = binf)
  
  ##------convert NA values to "NA" for matrix ----
  binf$drug[is.na(binf$drug)] <- "NA"
  binf[, group.by][is.na(binf[, group.by])] <- "NA"
  colsName <- unique(binf[, group.by])
  rowsName <- unique(binf$drug)
  mat <-
    matrix(
      data = NA,
      nrow = length(rowsName),
      ncol = length(colsName),
      dimnames = list(rowsName, colsName)
    )
  for (i in 1:nrow(binf))
  {
    mat[binf[i, "drug"], binf[i, group.by]] <-
      binf[i, response.measure]
  }
  
  rt$mat <- mat
  
  return(rt)
}


.checkResMes <- function(object, response.measure)
{
  rm.type <- NULL
  if (response.measure %in% colnames(slot(object, "sensitivity")[["model"]]))
  {
    rm.type <- "model"
  }
  
  if (response.measure %in% colnames(slot(object, "sensitivity")[["batch"]]))
  {
    rm.type <- "batch"
  }
  
  if (is.null(rm.type))
  {
    msg <-
      sprintf(
        "valid response.measure values are\nFor model: %s\n\nFor batch: %s\n",
        paste0(colnames(slot(
          object, "sensitivity"
        )[["model"]]), collapse = ", "),
        paste0(colnames(slot(
          object, "sensitivity"
        )[["batch"]]), collapse = ", ")
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
#' @param other.col Names of other columns to add in the return, if return.type is data.frame
#' @param return.type Specify return type, allowed values are 'matrix' or 'data.frame'. Default \code{'matrix'}
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
#'
#'
#' @examples
#' data(brca)
#' brca.mR <- summarizeResponse(brca, response.measure = "mRECIST", group.by="patient.id")
#'
#' ## to get batch level response
#' br <- summarizeResponse(brca, response.measure = "abc", return.type="data.frame")
#'
#' ## to get model level response for a batch
#' r <- summarizeResponse(brca, response.measure = "mRECIST", batch.id = "X-1004.BGJ398")
#'
#' @export
summarizeResponse <- function(object,
                              response.measure = "mRECIST",
                              model.id = NULL,
                              batch.id = NULL,
                              group.by = "patient.id",
                              summary.stat = c(";", "mean", "median"),
                              tissue = NULL,
                              #summarize.batch=FALSE,
                              other.col = NULL,
                              return.type = c("matrix", "data.frame"))
{
  summary.stat <- match.arg(summary.stat)
  return.type <- match.arg(return.type)
  
  rm.type <- .checkResMes(object, response.measure)
  
  if (rm.type == "model")
  {
    if (!is.null(batch.id))
    {
      mat <- .ModelResponsePerBatch(object, response.measure, batch.id)
      return(mat)
    } else {
      rt <-
        .summarizePerModelResponse(
          object,
          response.measure = response.measure,
          model.id = model.id,
          group.by = group.by,
          summary.stat = summary.stat,
          other.col = other.col,
          tissue = tissue
        )
      
      if (return.type == "matrix")    {
        return(rt$mat)
      }
      if (return.type == "data.frame") {
        rownames(rt$df) <- NULL
        return(rt$df)
      }
    }
  }
  
  if (rm.type == "batch")
  {
    if (!is.null(other.col))
    {
      otcol <- unique(c("patient.id", "drug", other.col))
    } else {
      otcol <- c("patient.id", "drug")
    }
    
    rt <-
      .summarizePerBatchResponse(
        object,
        response.measure = response.measure,
        other.col = otcol,
        group.by = group.by,
        batch.name = batch.id
      )
    
    if (return.type == "matrix")    {
      return(rt$mat)
    }
    if (return.type == "data.frame") {
      rownames(rt$df) <- NULL
      return(rt$df)
    }
  }
}
