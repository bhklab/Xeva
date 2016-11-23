

checkModel <- function(model, expSlot)
{
  reqColName = c("model.id", "batch", "exp.type", "biobase.id")#, "drug.join.name")
  if(all(reqColName %in% colnames(model))==FALSE)
  {
    msg = sprintf("The required colmns for model are\n%s", paste(reqColName, collapse = ', '))
    stop(msg)
  }

  ##---- add drug column ---------------
  model$drug.join.name = sapply(model$model.id, function(x){ expSlot[[x]]$drug$join.name})
  return(model)
}

creatExperimentDesign <- function(model, expSlot)
{
  drgNames  = unique(sapply(expSlot, "[[", c("drug", "join.name")))
  tretBatch = unique(model[,  "batch"])

  rtx=list()
  for(drgI in drgNames)
  {
    for(batI in tretBatch)
    {
      Lx = list(drug.join.name = drgI, batch = batI)
      Lx$treatment = unique(model[model$drug.join.name == drgI &
                                    model$batch == batI &
                                    model$exp.type == "treatment", "model.id"] )

      Lx$control = unique(model[model$batch == batI & model$exp.type == "control", "model.id"] )

      if(length(Lx$treatment) > 0)
      {
        namx = sprintf("%s.%s", drgI, batI)
        rtx[[namx]]= Lx
      }
    }
  }

  trG = sapply(rtx, "[[", "treatment", simplify = TRUE)
  cnG = sapply(rtx, "[[", "control", simplify = TRUE)
  modInExpDes = unique(unlist(list(trG,cnG), recursive=TRUE))
  modInExpSlot= sapply(expSlot, "[[", "model.id")
  modInExpSlotNotPre = modInExpSlot[!modInExpSlot %in% modInExpDes]

  if(length(modInExpSlotNotPre)>0)
  {
    msg = sprintf("These model ids are not present in experiment design\n%s", paste(modInExpSlotNotPre, collapse = ', '))
    stop(msg)
  }

  return(rtx)
}

##================================================================================================


#' @export
PharmacoPSet <- setClass( "PharmacoPSet",
                          slots = list( molecularProfiles = "list",
                                     model = "data.frame",
                                     drug = "data.frame",
                                     experiment = "list",
                                     expDesign = "list")
                          )




#' if model and drug slot is NULL it will try to infer it
#' otherwise will use the given model and drug slot
#' @export
creatPharmacoPXSet <- function(name,
                               molecularProfiles = list(),
                               experiment = data.frame(),
                               model = data.frame(),
                               drug  = data.frame())
{

  expSlot = experimentSlotfromDf(experiment)

  model = checkModel(model, expSlot)

  expDesign = creatExperimentDesign(model, expSlot)

  pxset = PharmacoPSet(molecularProfiles = molecularProfiles,
                       model = model,
                       drug  = drug,
                       experiment= expSlot,
                       expDesign = expDesign )
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

