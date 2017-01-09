.subsetExperimentSlot <- function(object, id, id.type)
{
  rtx = list()
  for(Ix in object@experiment)
  {
    if(is.element(id, Ix[[id.type]])==TRUE)
    {rtx = .appendToList(rtx, Ix) }
  }
  return(rtx)
}

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
.getExperimentDataFromAExpID <- function(object, model.id)
{
  #modX = .subsetExperimentSlot(object, id=experiment.id, id.type="experiment.id")
  modX = .subsetExperimentSlot(object, id=model.id, id.type="model.id")
  if(length(modX)==1)
  { mod = modX[[1]] }else
  {return(NULL)}

  mod.data = mod$data
  mod$data = NULL
  modDf = reshape2::melt(mod)
  collpsCol = sort(colnames(modDf)[colnames(modDf)!="value"])
  collpsVal = apply(modDf[,collpsCol], 1, pasteWithoutNA, collapse="." )
  modDf[,"value.name"] = collpsVal
  for(I in 1:nrow(modDf))
  { mod.data[, modDf[I, "value.name"]] = as.character(modDf[I, "value"]) }

  mod.data = .removeNAcol(mod.data)
  mod.data = .reorderCol(mod.data, "model.id", 1)
  mod.data = .reorderCol(mod.data, "drug.join.name", 2)
  return(mod.data)
}

##--------------------------------------------------------------------------------------------------
##----- get experiment data in flat data.fram ------------------------------------------------------
#' For a given  model.id, it will return a data.fram
#' containing all data stored in experiment slot
#'
#' @examples
#' data(pdxe)
#' getExperiment(pdxe, model.id="X.1004.pael")
#' @param object The \code{XevaSet}
#' @param model.id The \code{model.id} for which data is required
#' @return a \code{data.fram} will all the the values stored in experiment slot
setGeneric(name = "getExperiment", def = function(object, model.id){standardGeneric("getExperiment")} )

#' @export
setMethod( f=getExperiment,
           signature="XevaSet",
           definition=function(object, model.id)
           {
             model.ids = unique(c(model.id))
             rtx = lapply(model.ids, .getExperimentDataFromAExpID, object=object)
             rtz = .rbindListOfDataframs(rtx)
             return(rtz)
           })






