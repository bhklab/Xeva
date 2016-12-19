.generateNewNameForBatch <- function(inVec, txt)
{
  for(I in (length(inVec)+1):(3*length(inVec)))
  {
    tn = sprintf("%s.%d", txt,I)
    if(is.element(tn, inVec)==FALSE){return(tn)}
  }
}


#' add a new experimental design
#'
#' @examples
#' data(pdxe)
#' @param object The \code{Xeva} dataset
#' @param treatment The \code{model.id} of treatment
#' @param control The \code{model.id} of control
#' @param batch.name The \code{batch.name} for new batch
#' @param replace If TRUE will replace the old batch with new values
#' @return returns \code{Xeva} dataset with new experimental design added
setGeneric(name = "addExperimentalDesign", def = function(object, treatment, control=NULL, batch.name=NULL,replace=FALSE) {standardGeneric("addExperimentalDesign")} )

#' @export
setMethod( f=addExperimentalDesign,
           signature=c(object="XevaSet"),
           definition= function(object, treatment, control=NULL, batch.name=NULL, replace=FALSE)
           {

             allBatchName = sapply(object@expDesign, '[[', "batch.name")
             if(!is.null(batch.name))
             {
               if(is.element(batch.name, allBatchName)==TRUE)
               {
                 if(replace==FALSE){
                 msg = sprintf("\nbatch.name %s already exist\nPlease give a new name\n", batch.name)
                 stop(msg)}

                 if(replace==TRUE){
                   msg = sprintf("\nThis will replace the old batch.name %s\n", batch.name)
                   cat(msg)}
               }
             } else
             { batch.name = .generateNewNameForBatch(allBatchName, "batch")}

             object@expDesign = .appendToList(object@expDesign,
                                              list(batch.name=batch.name,
                                                   treatment=c(treatment),
                                                   control=c(control) ))
             return(object)
           })






addExperiment <- function()
{
  ##this will add an Experiment in dataset
}

