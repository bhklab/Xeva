##--------------------------------------------------------------------------
##------------ To create Sensitivity Slot ----------------------------------
.checkUnqLength <- function(inVec)
{ length(inVec)== length(unique(inVec)) }

.creatSensitivitySlot <- function(modelSensitivity, batchSensitivity, expSlot,
                                  expDesign)
{
  ##------------ for modelSensitivity ------------------------------------
  if(nrow(modelSensitivity)==0)
  {
    modelSensitivity <-data.frame(model.id= names(expSlot), stringsAsFactors = FALSE)
  }

  for(mid in names(expSlot))
  {
    if(is.element(mid, modelSensitivity$model.id) ==FALSE)
    {
      msg <- sprintf("provide modelSensitivity for all models\nmodelSensitivity missing for %s", mid)
      stop(msg)
    }
  }

  if( .checkUnqLength(modelSensitivity$model.id)==FALSE)
  {stop("model.ids are not unique")}
  rownames(modelSensitivity) <- as.character(modelSensitivity$model.id)
  modelSensitivity <- modelSensitivity[names(expSlot), ,drop=FALSE]
  ##--------------------------------------------------------------------------
  ##------------ for Batch Sensitivity ---------------------------------------

  if(nrow(batchSensitivity)>0)
  {
    if(is.element("batch.name", colnames(batchSensitivity)) ==FALSE)
    {
      stop("in 'batchSensitivity' datafram one column must be 'batch.name'")
    }

    if(nrow(batchSensitivity)!= length(names(expDesign)))
    {
      batchSensitivity <- .reorderCol(batchSensitivity, "batch.name", 1)
      missingId <- unique(setdiff(names(expDesign), batchSensitivity$batch.name))
      bsN <- data.frame(matrix(NA, nrow = length(missingId),
                               ncol = ncol(batchSensitivity)))
      colnames(bsN) <- colnames(batchSensitivity)
      bsN$batch.name <- missingId
      batchSensitivity <- rbind(batchSensitivity, bsN)
    }
  }

  if(nrow(batchSensitivity)==0)
  {
    batchSensitivity <- data.frame(batch.name= names(expDesign),
                                   stringsAsFactors = FALSE)
  }

  if( .checkUnqLength(batchSensitivity$batch.name)==FALSE)
  {stop("batch names are not unique")}
  rownames(batchSensitivity) <- as.character(batchSensitivity$batch.name)
  batchSensitivity <- batchSensitivity[names(expDesign), ,drop=FALSE]

  if(is(modelSensitivity, "data.frame")==FALSE | is(batchSensitivity, "data.frame")==FALSE)
  {stop("slot class error")}

  rtx <- list(model = modelSensitivity,
              batch = batchSensitivity)
  return(rtx)
}



##--------------------------------------------------------------------------
##------------ To check the input parameters--------------------------------
.checkModel <- function(model, expSlot)
{
  reqColName <- c("model.id", "patient.id")
  if(all(reqColName %in% colnames(model))==FALSE)
  {
    msg <- sprintf("The required colmns for model are\n%s", paste(reqColName, collapse = ', '))
    stop(msg)
  }

  for(I in expSlot)
  {
    if(is.element(slot(I, "model.id"), model$model.id)==FALSE)
    {
      msg = sprintf("No informaton present in Model datafram about model.id =%s", I$model.id)
      stop(msg)
    }
  }

  mdup <- model$model.id[duplicated(model$model.id)]
  if(length(mdup)>0)
  {
    msg <- sprintf("duplicated model.id in model slot:\n%s\n", paste(mdup, collapse = "\n"))
    stop(msg)
  }
  rownames(model) <- as.character(model$model.id)
  return(model)
}


.checkDrugSlot <- function(drf)
{
  if(is(drf, "data.frame"))
  {
    if(!"drug.id" %in% colnames(drf))
    {
      stop("drug data.frame must have column drug.id")
    }
    drf <- data.frame(apply(drf, 2, as.character), stringsAsFactors = FALSE)
  } else {stop("drug not in data.frame")}
  return(drf)
}

.checkmodToBiobaseMapSlot <- function(modToBiobaseMap, molecularProfiles)
{
  if(!is.null(modToBiobaseMap) & nrow(modToBiobaseMap) > 0 & length(molecularProfiles)>0)
  {
    rqdCol <- c("model.id", "biobase.id", "mDataType")
    for(cx in rqdCol)
    {
      if(is.element(cx, colnames(modToBiobaseMap))==FALSE)
      {
        msg <- sprintf("column %s not present is modToBiobaseMap\nmodToBiobaseMap must have the columns\n%s\n",
                       cx, paste(rqdCol, collapse = "\n"))
        stop(msg)
      }
    }

    mbDataTypes <- unique(as.character(modToBiobaseMap$mDataType))
    w <- names(molecularProfiles)[!(names(molecularProfiles) %in% mbDataTypes)]
    if(length(w)>0)
    {
      msg <- sprintf("Id mapping for molecular data type %s not present in modToBiobaseMap", paste(w, collapse = "\n"))
      warning(msg)
    }
  } else
  {
    if(nrow(modToBiobaseMap) > 0)
    { warning("modToBiobaseMap not present")}

    if(length(molecularProfiles)>0)
    { warning("molecularProfiles not present")}
  }

  return(modToBiobaseMap)
}


##-------------------------------------------------------------------------
##--------- An S4 class for XevaSet ---------------------------------------
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
#' @param expDesign A list containing name of the batch, control and treatment model.id
#' @param modelSensitivity A \code{data.frame} containing sensitivity for each model
#' @param batchSensitivity A \code{data.frame} containing sensitivity for each batch
#' @param molecularProfiles A \code{list} of \code{ExpressionSet} objects containing
#'   different molecular profiles.
#' @param modToBiobaseMap A \code{data.frame} containing model.id corresponding Biobase object id and name of the molecularProfiles
#'
#' @return  Returns Xeva object
#'
#' @details This function creates a XevaSet object. It takes different model
#' infromation and genomic data as input. For detailed discription of all
#' varaibles please see Xeva vignette section \strong{"Creating new Xeva object"}
#'
#' @examples
#' ## read raw data files containg PDX experiment information and genomic data
#' model = read.csv(system.file("extdata", "model.csv", package = "Xeva"))
#' drug = read.csv(system.file("extdata", "drug.csv", package = "Xeva"))
#' experiment= read.csv(system.file("extdata", "experiments.csv", package = "Xeva"))
#' expDesign=readRDS(system.file("extdata", "batch_list.rds", package = "Xeva"))
#' RNASeq=readRDS(system.file("extdata", "rnaseq.rds", package = "Xeva"))
#' modToBiobaseMap=read.csv(system.file("extdata", "modelToExpressionMap.csv", package = "Xeva"))
#'
#' ## create Xeva object
#' xeva.set = createXevaSet(name="example xevaSet", model=model, drug=drug,
#'                          experiment=experiment, expDesign=expDesign,
#'                          molecularProfiles=list(RNASeq = RNASeq),
#'                          modToBiobaseMap = modToBiobaseMap)
#' print(xeva.set)
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

  cat(sprintf("##--- making experiment slot -------\n"))
  expSlot <- experimentSlotfromDf(experiment)

  cat(sprintf("##--- making model information -------\n"))
  model <- .checkModel(model, expSlot)

  cat(sprintf("##--- making experiment design slot -------\n"))
  expDesign <- .checkExperimentDesign(expDesign)

  cat(sprintf("##--- making sensitivity slot -------\n"))
  sensitivity <-
    .creatSensitivitySlot(modelSensitivity, batchSensitivity, expSlot, expDesign)

  cat(sprintf("##--- making drug slot -------\n"))
  drug <- .checkDrugSlot(drug)

  cat(sprintf("##--- making id mapping slot -------\n"))
  modToBiobaseMap <-
    .checkmodToBiobaseMapSlot(modToBiobaseMap, molecularProfiles)

  cat(sprintf("##--- creating XevaSet -------\n"))
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



#' A method to display object
#' for "show" setGeneric is already defined.
#' @param object A Xeva object
#' @import methods
#' @noRd
setMethod(
  f = "show",
  signature = "XevaSet",
  definition = function(object)
  {
    modSen <- colnames(slot(object, "sensitivity")$model)
    modSentxt <- ""
    if(length(modSen)>1)
    { modSentxt <- paste0(modSen[modSen!="model.id"], collapse = ", ") }

    batSen <- colnames(slot(object, "sensitivity")$batch)
    batSentxt <- ""
    if(length(batSen)>1)
    { batSentxt <- paste0(batSen[batSen!="batch.name"], collapse = ", ") }

    msg <-
      sprintf(
        "XevaSet\nname: %s\nCreation date: %s\nNumber of models: %d\nNumber of drugs: %d\nMoleculer dataset: %s\nSensitivity\nmodel:%s\nbatch:%s\n",
        slot(object, "annotation")$name,
        slot(object, "annotation")$dateCreated,
        length(slot(object, "experiment")),
        nrow(slot(object, "drug")),
        paste(names(slot(
          object, "molecularProfiles"
        )), collapse = ", "),
        modSentxt,
        batSentxt
      )
    cat(msg)
  }
)





#' print Xeva object
#'
#' \code{print} displays Xeva object information or model or batch information
#'
#' @param object \code{Xeva} object
#' @param id default \code{NULL}, id can be \code{model.id} or \code{batch.name}
#'
#' @return  Prints object, model or batch information.
#'
#' @examples
#' \dontrun{
#' data(brca)
#' print(brca)
#'
#' # to print a model
#' print(brca, id = "X.1004.BG98")
#'
#' # to print a batch
#' print(brca, id = "X-1004.BGJ398")
#' }
#' @keywords internal
#' @noRd
print.XevaSet <- function(object, id=NULL)
{
  if(is.null(id))
  {
    show(object)
  } else
  {
    if(is.character(id)==FALSE)
    {
      msg <- sprintf("id should be character type")
      stop(msg)
    }
    mod <- slot(object, "experiment")[[id]]
    if(!is.null(mod))
    {show(mod)} else
    {
      mod <- slot(object, "expDesign")[[id]]
      print(mod)
    }
  }
}
