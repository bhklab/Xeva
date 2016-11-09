

#experimentInformation

#' @export
PharmacoPSet <- setClass( "PharmacoPSet",
                          slots = c( molecularProfiles = "list",
                                     model = "data.frame",
                                     drug = "data.frame",
                                     experiment = "list" )
                          )


#' @export
formatExperimentSlot <- function(experiment, model=NULL, drug=NULL)
{
  ##-- the table experiment must have these colomns
  ##-- model.id
  requredCols = c("model.id", "drug", "time", "volume","body.weight")
  #if(colnames(experiment))

  Model.IDS = unique(experiment$model.id)
  unqExperiments = unique(experiment[, c("model.id","drug")])


}

#' if model and drug slot is NULL it will try to infer it
#' otherwise will use the given model and drug slot
#' @export
creatPharmacoPXSet <- function(molecularProfiles,
                               experiment,
                               model=NULL,
                               drug=NULL)
{

  molecularProfiles = list(a=1:10)
  model = data.frame(1:5, 1:5)
  drug = data.frame(1:5, 1:5)
  experiment = list(a=1:5)
  pxset = PharmacoPSet(molecularProfiles = molecularProfiles,
                       model = model,
                       drug  = drug,
                       experiment = experiment)
  return(pxset)
}




#' Add together two numbers.
#'
#' @param x A number.
#' @param y A number.
#' @param Whatever u put here
#' @return The sum of \code{x} and \code{y} and \code{z} . z is nothing
#' @examples
#' add(1, 1)
#' This will be also included
#' @export
documentationExample <- function(x, y) {
  x + y
}

