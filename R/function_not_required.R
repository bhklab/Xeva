.foo<- function()
{
if(1==2)
{

.getModelIdIndexInExpDesign <- function(object, model.id)
{
  cntrL = c(); tretL = c()
  for(i in 1:length(object@expDesign))
  {
    ed = object@expDesign[[i]]
    if(is.element(model.id, ed$control))
    { cntrL = c(cntrL, i) }

    if(is.element(model.id, ed$treatment))
    { tretL = c(cntrL, i) }
  }
  return(list(control.indx=cntrL, treat.indx=tretL))
}

#' Get experiment type (treatment or control) for a given model.id
#'
#' @examples
#' data(pdxe)
#' # get experiment type for model.id
#' experimentType(object=pdxe, model.id="X.1655.LE11.biib")
#' experimentType(object=pdxe, model.id="X.1655.uned")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return returns \code{treatment} or \code{control}
setGeneric(name = "experimentType", def = function(object, model.id) {standardGeneric("experimentType")} )

#' @export
setMethod( f=experimentType,
           signature="XevaSet",
           definition= function(object,model.id)
           {
             ct.indx = .getModelIdIndexInExpDesign(object, model.id)
             if(length(ct.indx$control.indx)>0 & length(ct.indx$treat.indx)==0)
             {return("control")}

             if(length(ct.indx$control.indx)==0 & length(ct.indx$treat.indx)>0)
             {return("treatment")}

             if(length(ct.indx$control.indx)>0 & length(ct.indx$treat.indx)>0)
             {return("control and treatment")}
             return(NA)
           })






## checks if a variable exist in Experiment slot
checkExperimentSlotVariable <- function(object, value)
{
  vars = lapply(slot(object, "experiment"), names)
  obj.vars = unique(unlist(vars))
  if(is.element(value, obj.vars)==FALSE)
  {
    msg = sprintf("%s is not part of experiment slot.
                  Avalible variables are\n%s", value, paste(obj.vars, collapse="\n"))
    stop(msg)
  }
  return(TRUE)
}



}
}
