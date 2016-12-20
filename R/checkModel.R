
.checkModel <- function(model, expSlot)
{
  reqColName = c("model.id", "biobase.id")
  if(all(reqColName %in% colnames(model))==FALSE)
  {
    msg = sprintf("The required colmns for model are\n%s", paste(reqColName, collapse = ', '))
    stop(msg)
  }

  ##----cheack if all expSlot model.id is present in model
  for(I in expSlot)
  {
    if(is.element(I$model.id, model$model.id)==FALSE)
    {
      msg = sprintf("No informaton present in Model datafram about model.id =%s", I$model.id)
      stop(msg)
    }

  }
  if(("patient.id" %in% colnames(model)) ==FALSE)
  {model$patient.id = model$biobase.id}

  rownames(model) = NULL

  return(model)
}
