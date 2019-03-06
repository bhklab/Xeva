#' Get batch information
#'
#' Get batch information from a Xeva dataset.
#'
#' @param object The Xeva object from which batch information is obtained.
#' @param batch Name of the batch. Default \code{NULL}.
#' @param model.id Model ID for which need to be searched in the batches. Default \code{NULL}.
#' @param model.id.type Type of the model ID in a batch. See the Details section below.
#'
#' @details By default this function will return the names of all the batches present in the
#' dataset. If a batch specified, it will return the experiment design (control
#' and treatment model IDs) of that particular batch. If \code{model.id} is specified,
#' it will return the names of all the batches where this particuler \code{model.id} is present.
#' If both \code{batch} and \code{model.id} are specified, \code{batch} will take precedent.
#'
#' For \code{model.id.type}, the default value \code{'any'} will return all the batch IDs
#' where the given model ID is present in any arm (ie. control or treatment) of the
#' batch. It can also be restricted to look only for treatment (or control) arm by
#' specifying the type.
#'
#'
#' @examples
#' data(brca)
#' ##to get all the batch names
#' batch.name <- batchInfo(brca)
#'
#' ##to get a specific batch
#' batch.design <- batchInfo(brca, batch=batch.name[1])
#'
#' ##to get all the batches where a model.id is present
#' batchInfo(brca, model.id="X.6047.uned")
#'
#' @return A \code{Vector} with batch names.
#'
#' @name batchInfo
setGeneric(name= "batchInfo",
           def = function(object, batch=NULL, model.id=NULL,
                          model.id.type=c("any", "control", "treatment"))
{standardGeneric("batchInfo")} )

#' @rdname batchInfo
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
                   msg <- sprintf("batch name %s not present\nuse batchInfo(object) to see all batch names", bn)
                   stop(msg)
                 }
                 btRTX[[bn]] <- bt
               }
               return(btRTX)
             }
           })
