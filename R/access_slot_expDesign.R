# batchNames<-function(object, model.id=NULL)
# {
# if(is.null(model.id))
# {
#  #rtx <- names(expDesign(object))
#   rtx <- names(slot(object, "expDesign"))
#  return(rtx)
# } else
# {
#   rtx <- list()
#   for(ed in slot(object, "expDesign"))
#   {
#     if(is.element(model.id, ed$treatment) | is.element(model.id, ed$control))
#     { rtx <- .appendToList(rtx, ed$batch.name) }
#   }
#   return(unlist(rtx))
# }
# }
# expDesign<- function(object, batch.name=NULL)
# {
#   if(is.null(batch.name))
#   {
#     return(slot(object, "expDesign"))
#   } else
#   {
#    btRTX <- list()
#    for(bn in c(batch.name))
#    {
#      bt <- slot(object, "expDesign")[[bn]]
#      if(is.null(bt))
#      {
#        msg <- sprintf("batch name %s not present\nuse batchNames(object) to see all batch names", batch.name)
#        stop(msg)
#      }
#      btRTX[[bn]] <- bt
#    }
#    return(btRTX)
#   }
# }



##------------------------------------------------------------------------------

#' Get batch information
#'
#' Get batch information from a Xeva dataset. By default it will return the names of all the batch present in the data-set.
#' \code{model.id} is specified, will return all batch names conting that model.id
#' @param object xeva object
#' @param batch name of the batch. Default \code{NULL}
#' @param model.id model id for which need to be searched in the batches. Default \code{NULL}
#' @param model.id.type type of the model id in a batch. See details
#'
#' @details By default it will return the names of all the batch present in the
#' data-set. If batch specified it will retun the experiment design (control
#' and treatment model ids) of that perticular batch. If model.id is specified
#' it will return names of all the batches where this perticuler model.id is present.
#'
#' For model.id.type the default value \code{'any'} will return all batch ids
#' where given model id is present in any arm (control or treatment) of the
#' batch. It can be restriceted to look only for treatment (or control) arm by
#' specifying the type.
#'
#' @examples
#' data(brca)
#' ##to get all the batch names
#' batch.name <- batchInfo(brca)
#'
#' ##to get an specific batch
#' batch.design <- batchInfo(brca, batch=batch.name[1])
#'
#' ##to get all the batchs where a model.id is present
#' batchInfo(brca, model.id="X.6047.uned")
#'
#' @return A \code{Vector} with batch names
setGeneric(name= "batchInfo", def = function(object, batch=NULL, model.id=NULL,
                                             model.id.type=c("any", "control", "treatment"))
{standardGeneric("batchInfo")} )
#' @export
setMethod( f="batchInfo",
           signature=c(object = "XevaSet"),
           definition=function(object, batch=NULL, model.id=NULL,
                               model.id.type=c("any", "control", "treatment"))
           {
             if(is.null(batch) & is.null(model.id))
             {
               rtx <- names(slot(object, "expDesign"))
               return(rtx)
             }

             if(is.null(batch) & !is.null(model.id))
             {
               model.id.type <- match.arg(model.id.type)
               rtx <- list()
               for(ed in slot(object, "expDesign"))
               {
                 if(model.id.type=="any")
                 {
                   if(is.element(model.id, ed$treatment) |
                      is.element(model.id, ed$control))
                   { rtx <- .appendToList(rtx, ed$batch.name) }
                 }

                 if(model.id.type=="control" & is.element(model.id, ed$control))
                 { rtx <- .appendToList(rtx, ed$batch.name) }

                 if(model.id.type=="treatment" & is.element(model.id, ed$treatment))
                 { rtx <- .appendToList(rtx, ed$batch.name) }
               }

               return(unique(unlist(rtx)))
             }

             if(!is.null(batch) & is.null(model.id))
             {
               btRTX <- list()
               for(bn in c(batch))
               {
                 bt <- slot(object, "expDesign")[[bn]]
                 if(is.null(bt))
                 {
                   msg <- sprintf("batch name %s not present\nuse batchInfo(object) to see all batch names", batch.name)
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
