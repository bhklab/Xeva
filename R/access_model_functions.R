accessModel <- function(object, model.id)
{
  model.id = c(model.id)
  return(object@model[object@model$model.id %in% model.id, ])
}
