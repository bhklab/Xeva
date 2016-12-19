
.checkModel <- function(model, expSlot)
{
  reqColName = c("model.id", "biobase.id")
  if(all(reqColName %in% colnames(model))==FALSE)
  {
    msg = sprintf("The required colmns for model are\n%s", paste(reqColName, collapse = ', '))
    stop(msg)
  }

  if(("patient.id" %in% colnames(model)) ==FALSE)
  {model$patient.id = model$biobase.id}

  ##---- add drug column ---------------
  #model$drug.join.name = sapply(model$model.id, function(x){ expSlot[[x]]$drug$join.name})
  return(model)
}
