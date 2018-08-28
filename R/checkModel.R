
.checkModel <- function(model, expSlot)
{
  reqColName <- c("model.id", "patient.id")
  if(all(reqColName %in% colnames(model))==FALSE)
  {
    msg <- sprintf("The required colmns for model are\n%s", paste(reqColName, collapse = ', '))
    stop(msg)
  }

  ##----cheack if all expSlot model.id is present in model slot
  for(I in expSlot)
  {

    #if(is.element(I$model.id, model$model.id)==FALSE)
    if(is.element(slot(I, "model.id"), model$model.id)==FALSE)
    {
      msg = sprintf("No informaton present in Model datafram about model.id =%s", I$model.id)
      stop(msg)
    }
  }

  ##----- set model slot rownames as model.id ----------------
  mdup <- model$model.id[duplicated(model$model.id)]
  if(length(mdup)>0)
  {
    msg <- sprintf("duplicated model.id in model slot:\n%s\n", paste(mdup, collapse = "\n"))
    stop(msg)
  }
  rownames(model) <- as.character(model$model.id)
  return(model)
}
