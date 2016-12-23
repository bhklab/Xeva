

##=====================================================================
#' An S4 class for XevaSet
#'
XevaSet <- setClass( "XevaSet",
                          slots = list(annotation = "list",
                                       molecularProfiles = "list",
                                       model = "data.frame",
                                       drug = "data.frame",
                                       experiment = "list",
                                       expDesign = "list")
                          )

######################################################################
#' Creat Xeva class object
#'
#' \code{creatPharmacoPxSet} returns Xeva class object
#'
#' @param name A \code{character} string detailing the name of the dataset
#' @param molecularProfiles A \code{list} of ExpressionSet objects containing molecular profiles
#' @param experiment A \code{data.frame} containg all experiment information
#' @param model A \code{data.frame} containg the annotations for all models used in the experiment
#' @param drug A \code{data.frame} containg the annotations for all the drugs profiled in the data set, across all data types
#'
#' @return  Returns Xeva object
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
creatXevaSet <- function(name,
                         molecularProfiles = list(),
                         experiment = data.frame(),
                         expDesign  = list(),
                         model = data.frame(),
                         drug  = data.frame())
{

  annotation <- list( name = as.character(name),
                      dateCreated = date(),
                      sessionInfo = sessionInfo())

  expSlot = experimentSlotfromDf(experiment)

  model = .checkModel(model, expSlot)

  expDesign = .checkExperimentDesign(expDesign)

  ##----check if drug present in both drug slot and expSlot

  pxset = XevaSet(annotation=annotation,
                       molecularProfiles = molecularProfiles,
                       model = model,
                       drug  = drug,
                       experiment= expSlot,
                       expDesign = expDesign )
  return(pxset)
}



##===============================================================================
#' A method to display object
#' for "show" setGeneric is already defined
#' @export
setMethod(f="show",
          signature="XevaSet",
          definition= function(object)
          {
            slotsName = paste(slotNames(object), collapse = "\n")
            cat(sprintf("Slots are:\n%s\n", slotsName))
          }
          )

