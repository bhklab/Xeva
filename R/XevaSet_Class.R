##An S4 class for XevaSet
##
XevaSet <- setClass(
  "XevaSet",
  slots = list(
    annotation = "list",
    model = "data.frame",
    drug = "data.frame",
    sensitivity = "list",
    expDesign = "list",
    experiment = "list",
    molecularProfiles = "list",
    modToBiobaseMap = "data.frame"
  )
)


#' XevaSet constructor
#'
#' A constructor to create XevaSet. Only objects returned by this constructor
#' are expected to work with the XevaSet methods.
#'
#' @param name A \code{character} string detailing the name of the dataset.
#' @param model A \code{data.frame} containing the annotations for all the models used
#'   in the experiment.
#' @param drug A \code{data.frame} containing the annotations for all the drugs
#'   profiled in the dataset, across all data types.
#' @param experiment A \code{data.frame} containing all experiment information.
#' @param expDesign
#' @param modelSensitivity
#' @param batchSensitivity
#' @param molecularProfiles A \code{list} of \code{ExpressionSet} objects containing
#'   different molecular profiles.
#' @param modToBiobaseMap
#'
#' @return  Returns Xeva object.
#'
#' @examples
#' NULL
#'
#' @export
#' @import methods
createXevaSet <- function(name,
                         model = data.frame(),
                         drug  = data.frame(),
                         experiment = data.frame(),
                         expDesign  = list(),
                         modelSensitivity = data.frame(),
                         batchSensitivity = data.frame(),
                         molecularProfiles = list(),
                         modToBiobaseMap = data.frame())
{
  annotation <- list(
    name = as.character(name),
    dateCreated = date(),
    sessionInfo = sessionInfo()
  )

  expSlot <- experimentSlotfromDf(experiment)
  model <- .checkModel(model, expSlot)
  expDesign <- .checkExperimentDesign(expDesign)
  sensitivity <-
    .creatSensitivitySlot(modelSensitivity, batchSensitivity, expSlot, expDesign)
  drug <- .checkDrugSlot(drug)
  modToBiobaseMap <-
    .checkmodToBiobaseMapSlot(modToBiobaseMap, molecularProfiles)

  pxset <- XevaSet(
    annotation = annotation,
    model = model,
    drug  = drug,
    sensitivity = sensitivity,
    expDesign = expDesign,
    experiment = expSlot,
    molecularProfiles = molecularProfiles,
    modToBiobaseMap = modToBiobaseMap
  )
  return(pxset)
}



##==============================================================================
#setGeneric("show", function(object, model.id) { standardGeneric("show")})

#' A method to display object
#' for "show" setGeneric is already defined.
#' @import methods
setMethod(
  f = "show",
  signature = "XevaSet",
  definition = function(object)
  {
    msg <-
      sprintf(
        "Xeva-set name: %s\nCreation date: %s\nNumber of models: %d\nNumber of drugs: %d\nMoleculer dataset: %s\n",
        slot(object, "annotation")$name,
        slot(object, "annotation")$dateCreated,
        length(slot(object, "experiment")),
        nrow(slot(object, "drug")),
        paste(names(slot(
          object, "molecularProfiles"
        )), collapse = ", ")
      )
    cat(msg)
  }
)
