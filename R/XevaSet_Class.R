<<<<<<< HEAD
##An S4 class for XevaSet
##
=======


##=====================================================================
#' An S4 class for XevaSet
#'
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
XevaSet <- setClass( "XevaSet",
                     slots = list(annotation = "list",
                                  model = "data.frame",
                                  drug = "data.frame",
                                  sensitivity="list",
                                  expDesign = "list",
                                  experiment = "list",
                                  molecularProfiles = "list",
                                  modToBiobaseMap = "data.frame"
                                  )
                     )

######################################################################
#' Creat Xeva class object
#'
#' \code{creatXevaSet} returns Xeva class object
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
<<<<<<< HEAD
#' @import methods
=======
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
creatXevaSet <- function(name,
                         model = data.frame(),
                         drug  = data.frame(),
                         experiment = data.frame(),
                         expDesign  = list(),
                         modelSensitivity= data.frame(),
                         batchSensitivity= data.frame(),
                         molecularProfiles = list(),
                         modToBiobaseMap = data.frame())
{
  annotation <- list( name = as.character(name),
                      dateCreated = date(),
                      sessionInfo = sessionInfo())

  expSlot <- experimentSlotfromDf(experiment)
  model <- .checkModel(model, expSlot)
  expDesign <- .checkExperimentDesign(expDesign)
  sensitivity<- .creatSensitivitySlot(modelSensitivity, batchSensitivity, expSlot, expDesign)
  drug <- .checkDrugSlot(drug)
<<<<<<< HEAD
  modToBiobaseMap <- .checkmodToBiobaseMapSlot(modToBiobaseMap, molecularProfiles)
=======

  ##----check if drug present in both drug slot and expSlot
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4

  pxset = XevaSet(annotation = annotation,
                  model = model,
                  drug  = drug,
                  sensitivity = sensitivity,
                  expDesign = expDesign,
                  experiment= expSlot,
                  molecularProfiles = molecularProfiles,
                  modToBiobaseMap= modToBiobaseMap)
  return(pxset)
}



##===============================================================================
#' A method to display object
#' for "show" setGeneric is already defined
#' @export
<<<<<<< HEAD
#' @import methods
=======
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
setMethod(f="show",
          signature="XevaSet",
          definition= function(object)
          {
<<<<<<< HEAD
            msg <- sprintf(
"Xeva-set name: %s
Creation date: %s
Number of models: %d
Number of drugs: %d
Moleculer dataset: %s\n",
object@annotation$name, object@annotation$dateCreated,
length(object@experiment), dim(object@drug)[1],
paste(names(object@molecularProfiles), collapse = ", ")
)
            cat(msg)
          }
)



=======
            msg <- sprintf("Xeva-set name: %s\nCreation date: %s\nNumber of models: %d\nNumber of drugs: %d",
                           object@annotation$name, object@annotation$dateCreated,
                           length(object@experiment), dim(object@drug)[1])
            cat(msg)
            #slotsName <- paste(slotNames(object), collapse = "\n")
            #cat(sprintf("Slots are:\n%s\n", slotsName))
          }
)

>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
