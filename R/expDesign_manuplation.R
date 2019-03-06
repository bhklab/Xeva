.generateNewNameForBatch <- function(inVec, txt)
{
  for(I in (length(inVec)+1):(3*length(inVec)))
  {
    tn = sprintf("%s.%d", txt,I)
    if(is.element(tn, inVec)==FALSE){return(tn)}
  }
}


#' Add a new experimental design
#'
#'
#' Add a new experimental design in the \code{expDesign} slot.
#' @examples
#' data(brca)
#' brca <- addExperimentalDesign(object=brca, treatment=c("X.6047.LL71"),
#'         control=c("X.6047.uned"), batch.id="new.batch", replace=FALSE)
#'
#' @param object The \code{Xeva} dataset.
#' @param treatment The \code{model.id} of treatment.
#' @param control The \code{model.id} of control.
#' @param batch.id The \code{batch.id} for a new batch.
#' @param replace If \code{TRUE}, replace an old batch with new values.
#' @return Returns \code{Xeva} dataset with new experimental design added.
setGeneric(name = "addExperimentalDesign",
           def = function(object, treatment=NULL, control=NULL, batch.id=NULL,replace=FALSE)
                          {standardGeneric("addExperimentalDesign")} )


#' @rdname addExperimentalDesign
#' @export
setMethod( f=addExperimentalDesign,
           signature=c(object="XevaSet"),
           definition= function(object, treatment=NULL, control=NULL,
                                batch.id=NULL, replace=FALSE)
             {
              if(is.null(treatment) & is.null(control) )
              { stop("treatment and control both can't be NULL") }

              for(mi in unlist(c(treatment,control)))
              {
                if(!(mi %in% names(slot(object, "experiment")) ))
                {
                  stop(sprintf("model.id %s not present in dataset", mi))
                }
              }

             allBatchName <- vapply(slot(object, "expDesign"), function(x)
                                   {x$batch.name}, FUN.VALUE = character(1))
             if(!is.null(batch.id))
             {
               if(is.element(batch.id, allBatchName)==TRUE)
               {
                 if(replace==FALSE){
                 msg = sprintf("\nbatch.id %s already exist\nPlease give a new name\n", batch.id)
                 stop(msg)}

                 if(replace==TRUE){
                   msg = sprintf("\nThis will replace the old batch.id %s\n", batch.id)
                   cat(msg)}
               }
             }
             slot(object, "expDesign")[[batch.id]] <- list(batch.name=batch.id,
                                                   treatment=c(treatment),
                                                   control=c(control))
           return(object)
           })

