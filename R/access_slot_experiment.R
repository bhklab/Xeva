
.subsetExperimentSlotForDrug <- function(object, drugName, exact.match=TRUE)
{
  drgNames = stringr::str_trim(strsplit(drugName, "\\+")[[1]])
  expList = c()
  for(Ix in object@experiment)
  {
    if(exact.match==TRUE)
    {
      if(length(drgNames) == length(Ix$drug$names))
      {
        if(all(drgNames %in% Ix$drug$names)==TRUE)
        {expList = c(expList, Ix$model.id)}
      }
    }
    ###------------------------------------------------
    if(exact.match==FALSE)
    {
      if(any(drgNames %in% Ix$drug$names)==TRUE)
      {expList = c(expList, Ix$model.id)}
    }
    ###------------------------------------------------
  }
  rdx = data.frame(model.id=unique(expList), drug=drugName, stringsAsFactors = FALSE)
  return(rdx)
}

##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
##----- select model.id based on drug, tumor.type -----------------------------------------------------
#' To select model ids based on drug name and/or tumor type
#' @examples
#' data(pdxe)
#' selectModelIds(pdxe, drug="paclitaxel", drug.match.exact=TRUE, tumor.type="BRCA")
#' @param object The \code{XevaSet}
#' @param drug Name of the \code{drug}
#' @param drug.match.exact Default \code{TRUE}
#' @param tumor.type Tumor type. Default \code{NULL}
#'
#' @return a \code{vector} with the matched model.ids
#'
setGeneric(name = "selectModelIds",
           def = function(object,
                          drug=NULL, drug.match.exact=TRUE,
                          tumor.type=NULL)
           {standardGeneric("selectModelIds")} )

#' @export
setMethod( f=selectModelIds, signature="XevaSet",
           definition=function(object,
                               drug=NULL, drug.match.exact=TRUE,
                               tumor.type=NULL)
           {
             if(is.null(drug) & is.null(tumor.type))
             {stop("drug and tumor.type both NULL, Please provide atleast one")}

             ExpIdsDrug = NULL
             if(!is.null(drug))
             {
               drug = c(drug)
               ExpIdsDrug = .subsetExperimentSlotForDrug(object, drug, exact.match=drug.match.exact)
             }

             ExpIdsTumor = NULL
             if(!is.null(tumor.type))
             {
               tumor.type = c(tumor.type)
               ExpIdsTumor = mapModelSlotIds(object, id = tumor.type, id.name = "tumor.type", map.to="all")
             }


             if(!is.null(drug) & is.null(tumor.type))
             { return(ExpIdsDrug$model.id) }

             if(is.null(drug) & !is.null(tumor.type))
             { return(ExpIdsTumor$model.id) }

             if(!is.null(drug) & !is.null(tumor.type))
             {
               rtx = intersect(ExpIdsDrug$model.id, ExpIdsTumor$model.id)
               return(rtx)
             }
           })



###########################################################################
###-----------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------
.getExperimentDataFromAExpID <- function(object, model.id, treatment.only)
{
  #modX = .subsetExperimentSlot(object, id=model.id, id.type="model.id")
  #if(length(modX)==1)
  #{ mod = modX[[1]] }else
  #{return(NULL)}

  mod <- slot(object, "experiment")[[model.id]]
  if(is.null(mod))
  {
    msg <- sprintf("model.id '%s' not present in object", model.id)
    stop(msg)
  }

  mod.data <- mod$data
  mod.data$model.id <- mod$model.id
  mod.data$drug.join.name <- mod$drug$join.name

  mod.data = .removeNAcol(mod.data)
  mod.data = .reorderCol(mod.data, "model.id", 1)
  mod.data = .reorderCol(mod.data, "drug.join.name", 2)

  if(treatment.only==TRUE & !is.null(mod.data$dose))
  {
    tretIndx = extractBetweenTags(mod.data$dose, start.tag=0, end.tag=0)
    mod.data = mod.data[tretIndx, ]
  }

  mod.data$volume.normal <- NA
  if(mod.data$volume[1] > 0)
  { mod.data$volume.normal <- (mod.data$volume - mod.data$volume[1])/mod.data$volume[1]}
  return(mod.data)
}


.getCombinedDFforMultipalIDs <- function(object, model.ids, treatment.only)
{
  rv <- lapply(c(model.ids), function(modID)
  {.getExperimentDataFromAExpID(object=object, modID, treatment.only=treatment.only)})
  rv = .rbindListOfDataframs(rv)
}


.getExperimentForBatchName <- function(object, batch.name, treatment.only)
{
  rtxTr = NULL; rtxCn=NULL
  expDesign <- expDesign(object, batch.name)
  if(length(expDesign$treatment)>0)
  {
    rtxTr <- .getCombinedDFforMultipalIDs(object, expDesign$treatment, treatment.only)
    rtxTr$exp.type = "treatment"
  }

  if(length(expDesign$control)>0)
  {
    rtxCn <- .getCombinedDFforMultipalIDs(object, expDesign$control, treatment.only)
    rtxCn$exp.type = "control"
  }

  rv = rbind(rtxTr, rtxCn)
  return(rv)
}

##--------------------------------------------------------------------------------------------------
##----- get experiment data in flat data.fram ------------------------------------------------------
#' For a given  model.id, it will return a data.fram
#' containing all data stored in experiment slot
#'
#' @examples
#' data(pdxe)
#' getExperiment(pdxe, model.id="X.1004.pael", treatment.only=TRUE)
#' @param object The \code{XevaSet}
#' @param model.id The \code{model.id} for which data is required
#' @param batch.name The \code{batch.name} for which data is required
#' @param treatment.only Default \code{FALSE}. If TRUE give data only for non-zero dose periode (if dose data avalible)
#' @return a \code{data.fram} will all the the values stored in experiment slot
setGeneric(name = "getExperiment", def = function(object, model.id=NULL, batch.name=NULL,
                                                  treatment.only=FALSE)
  {standardGeneric("getExperiment")} )

#' @export
setMethod( f=getExperiment,
           signature="XevaSet",
           definition=function(object, model.id=NULL, batch.name=NULL, treatment.only=FALSE)
           {
             if(is.null(model.id) & is.null(batch.name))
             {
               msg = sprintf("model.id and batch.name both NULL")
               stop(msg)
             }
             if(!is.null(model.id))
             {
               model.ids <- unique(c(model.id))
               rtz <- .getCombinedDFforMultipalIDs(object, model.ids, treatment.only)
             }

             if(is.null(model.id) & !is.null(batch.name))
             {
               rtz <- .getExperimentForBatchName(object, batch.name, treatment.only)
             }
             return(rtz)
           })




##-----------------------------------------------------------------------

## checks if a variable exist in Experiment slot

checkExperimentSlotVariable <- function(object, value)
{
  vars = lapply(slot(object, "experiment"), names)
  obj.vars = unique(unlist(vars))
  if(is.element(value, obj.vars)==FALSE)
  {
    msg = sprintf("%s is not part of experiment slot.
                  Avalible variables are\n%s", value, paste(obj.vars, collapse="\n"))
    stop(msg)
  }
  return(TRUE)
}


##---------------------------------------------------------------------------------------
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
#' data(pdxe)
#' # to print object information
#' print(pdxe)
#'
#' # to print a model
#' model.id = modelInfo(pdxe)$model.id[1]
#' print(pdxe, id = model.id)
#'
#' # to print a batch
#' batch.id = batchNames(pdxe)[1]
#' print(pdxe, id = batch.id)
#' @export
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
    if(is.null(mod))
    {
      mod <- slot(object, "expDesign")[[id]]
    }
    print(mod)
  }
}

##--------------------------------------------------------------------------------

