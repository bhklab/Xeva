

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


#' An S4 class for PharmacoPSet
#'
#' @slot molecularProfiles List of molecular profiles
#' @slot model A datafram containing model infromation
PharmacoPSet <- setClass( "PharmacoPSet",
                          slots = list(annotation = "list",
                                       molecularProfiles = "list",
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

  annotation <- list( name = as.character(name),
                      dateCreated = date(),
                      sessionInfo = sessionInfo())

  expSlot = experimentSlotfromDf(experiment)

  model = checkModel(model, expSlot)

  expDesign = creatExperimentDesign(model, expSlot)

  pxset = PharmacoPSet(annotation=annotation,
                       molecularProfiles = molecularProfiles,
                       model = model,
                       drug  = drug,
                       experiment= expSlot,
                       expDesign = expDesign )
  return(pxset)
}



#' A method to display object
#' for "show" setGeneric is already defined
#' @export
setMethod(f="show",
          signature="PharmacoPSet",
          definition= function(object)
          {
            slotsName = paste(slotNames(object), collapse = "\n")
            cat(sprintf("Slots are:\n%s\n", slotsName))
          }
          )



##--------------------------------------------------------------------------

##----- get drugInfo -------------
#' drugInfo Generic
#' Generic for drugInfo method
#'
#' @examples
#' data(pdxe)
#' drugInfo(pdxe)
#' @param object The \code{PharmacoPSet} to retrieve drug info from
#' @return a \code{data.frame} with the drug annotations
setGeneric(name = "drugInfo", def = function(object) {standardGeneric("drugInfo")} )

#### @describeIn PharmacoSet Returns the annotations for all the drugs tested in the PharmacoSet
#' @export
setMethod( f=drugInfo, signature="PharmacoPSet", definition=function(object){ dim(object@drug) } )



#' drugInfo<- Generic
#' Generic for drugInfo replace method
#' @examples
#' data(pdxe)
#' drugInfo(pdxe) <- drugInfo(pdxe)
#' @param object The \code{PharmacoPSet} to replace drug info in
#' @param value A \code{data.frame} with the new drug annotations
#' @return Updated \code{PharmacoPSet}
setGeneric(name= "drugInfo<-", def = function(object, value) {standardGeneric("drugInfo<-")} )

###### @describeIn PharmacoSet Update the drug annotations
#' @export
setMethod( f="drugInfo<-",
           signature="PharmacoPSet",
           definition=function(object, value)
           {
             object@annotation$drugInfo = value
             #object@drugInfo = value  ##This will not work as slot drugInfo already have to be present
             object
           } )













#' Add together two numbers.
#'
#' @param x A number.
#' @param y A number.
#' @param Whatever u put here
#' @return The sum of \code{x} and \code{y} and \code{z} . z is nothing
#' @examples
#' documentationExample(1, 1)
#' This will be also included
#' @export
documentationExample <- function(x, y) {
  x + y
}

