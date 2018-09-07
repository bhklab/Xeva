#' Get batch names/ids
#'
#' Get batch names/ids from a Xeva dataset. If \code{model.id} is specified, will return all batch names conting that model.id
#' @examples
#' data(brca)
#' batchNames(brca)
#' batchNames(brca, model.id="X.6047.uned")
#' @param object A \code{XevaSet}
#' @param model.id default \code{NULL}. if specified it will return batch names conting that model.id
#'
#' @return A \code{Vector} with batch names
setGeneric(name= "batchNames", def = function(object, model.id=NULL)
  {standardGeneric("batchNames")} )
#' @export
setMethod( f="batchNames",
           signature=c(object = "XevaSet"),
           definition=function(object, model.id=NULL)
           {
            if(is.null(model.id))
            {
             rtx <- names(expDesign(object))
             return(rtx)
            } else
            {
              rtx <- list()
              for(ed in slot(object, "expDesign"))
              {
                if(is.element(model.id, ed$treatment) | is.element(model.id, ed$control))
                { rtx <- .appendToList(rtx, ed$batch.name) }
              }
              return(unlist(rtx))
            }
            })



#' Given a batch.name get batch
#'
#' Given a batch.name get batch from a Xeva dataset
#' @examples
#' data(brca)
#' expDesign(brca, batch.name = "X-6047.paclitaxel")
#' @param object \code{XevaSet}
#' @param object \code{batch.name}. If NULL will return all batch in the dataset
#' @return A \code{Vector} with all batch.name
setGeneric(name= "expDesign", def = function(object, batch.name=NULL)
  {standardGeneric("expDesign")} )

#' @export
setMethod(f="expDesign", signature=c(object = "XevaSet"),
          definition=function(object, batch.name=NULL)
          {
            if(is.null(batch.name))
            {
              return(slot(object, "expDesign"))
            } else
            {
             btRTX <- list()
             for(bn in c(batch.name))
             {
               bt <- slot(object, "expDesign")[[bn]]
               if(is.null(bt))
               {
                 msg <- sprintf("batch name %s not present\nuse batchNames(object) to see all batch names", batch.name)
                 stop(msg)
               }
               btRTX[[bn]] <- bt
             }
             return(btRTX)
            }
          })

##------------------------------------------------------------------------------


getBatchFormatted <- function(object, batch=NULL, patient.id=NULL, drug=NULL, control.name=NULL)
{
  if(!is.null(batch))
  {
    if(is.character(batch))
    {
      bt <- slot(object, "expDesign")[[batch]]
      if(is.null(bt))
      {
        msg <- sprintf("batch name %s not present in object", batch)
        stop(msg)
      }
      return(bt)
    }

    if(is.list(batch))
    {
      if(!"batch.name" %in% names(batch))
      { stop(sprintf("'batch.name' is required")) }

      if(is.null(batch$treatment) & is.null(batch$control))
      { stop(sprintf("'treatment' and 'control' both can't be NULL")) }

      allMod <- c(batch$treatment, batch$control)
      modNotPr <- setdiff(allMod, rownames(modelInfo(object)))
      if(length(modNotPr)>0)
      {
        stop(sprintf("model.id not present in dataset: %s", paste(modNotPr, collapse = ", ")))
      }

      return(batch)
    }
  } else {
    mid <- modelInfo(object)
    mid <- mid[mid$patient.id==patient.id, ]
    rtx <- list(name=patient.id)
    rtx$treatment <- mid[mid$drug==drug, "model.id"]
    if(!is.null(control.name))
    { rtx$control <- mid[mid$drug==control.name, "model.id"]}
    return(rtx)
  }
}
