
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


##================================================================================================


#' An S4 class for XenoSet
#'
XenoSet <- setClass( "XenoSet",
                          slots = list(annotation = "list",
                                       molecularProfiles = "list",
                                       model = "data.frame",
                                       drug = "data.frame",
                                       experiment = "list",
                                       expDesign = "list")
                          )

######################################################################
#' Creat XenoR class object
#'
#' \code{creatPharmacoPxSet} returns XenoR class object
#'
#' @param name A \code{character} string detailing the name of the dataset
#' @param molecularProfiles A \code{list} of ExpressionSet objects containing molecular profiles
#' @param experiment A \code{data.frame} containg all experiment information
#' @param model A \code{data.frame} containg the annotations for all models used in the experiment
#' @param drug A \code{data.frame} containg the annotations for all the drugs profiled in the data set, across all data types
#'
#' @return  Returns XenoR object
#'
#' @examples
#' geoExp = readRDS("DATA-raw/Geo_Exp.Rda")
#' pdxe = creatPharmacoPxSet(name = "PDXE",
#'                           molecularProfiles = list(RNASeq = geoExp$RNASeq),
#'                           experiment = geoExp$experiment,
#'                           model = geoExp$model,
#'                           drug  = geoExp$drug )
#' save(pdxe, file = "data/pdxe.rda")
#' data("pdxe")
#' @export
creatPharmacoPxSet <- function(name,
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

  expDesign = creatExperimentDesign(expSlot)

  pxset = XenoSet(annotation=annotation,
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
          signature="XenoSet",
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
#' @param object The \code{XenoSet} to retrieve drug info from
#' @return a \code{data.frame} with the drug annotations
setGeneric(name = "drugInfo", def = function(object) {standardGeneric("drugInfo")} )

#### @describeIn PharmacoSet Returns the annotations for all the drugs tested in the PharmacoSet
#' @export
setMethod( f=drugInfo, signature="XenoSet", definition=function(object){ dim(object@drug) } )



#' drugInfo<- Generic
#' Generic for drugInfo replace method
#' @examples
#' data(pdxe)
#' drugInfo(pdxe) <- drugInfo(pdxe)
#' @param object The \code{XenoSet} to replace drug info in
#' @param value A \code{data.frame} with the new drug annotations
#' @return Updated \code{XenoSet}
setGeneric(name= "drugInfo<-", def = function(object, value) {standardGeneric("drugInfo<-")} )

###### @describeIn PharmacoSet Update the drug annotations
#' @export
setMethod( f="drugInfo<-",
           signature="XenoSet",
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

