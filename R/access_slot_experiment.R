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
        {expList = c(expList, Ix$experiment.id)}
      }
    }
    ###------------------------------------------------
    if(exact.match==FALSE)
    {
      if(any(drgNames %in% Ix$drug$names)==TRUE)
      {expList = c(expList, Ix$experiment.id)}
    }
    ###------------------------------------------------
  }
  rdx = data.frame(experiment.id=unique(expList), drug=drugName, stringsAsFactors = FALSE)
  return(rdx)
}

##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
##----- get experiment Ids -----------------------------
#' getExperimentIds Generic
#' Generic for getExperimentIds method
#'
#' @examples
#' data(pdxe)
#' getExperimentIds(pdxe, drug="paclitaxel", drug.match.exact=TRUE, tumor.type="BRCA")
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{list} with the all experiment designs
setGeneric(name = "getExperimentIds",
           def = function(object,
                          drug=NULL, drug.match.exact=TRUE,
                          tumor.type=NULL)
                  {standardGeneric("getExperimentIds")} )

#' @export
setMethod( f=getExperimentIds, signature="XevaSet",
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
  { return(ExpIdsDrug$experiment.id) }

  if(is.null(drug) & !is.null(tumor.type))
  { return(ExpIdsTumor$experiment.id) }

  if(!is.null(drug) & !is.null(tumor.type))
  {
    rtx = intersect(ExpIdsDrug$experiment.id, ExpIdsTumor$experiment.id)
    return(rtx)
  }
})



###-----------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------
.getExperimentDataFromAExpID <- function(object, experiment.id)
{
  modX = .subsetExperimentSlot(object, id=experiment.id, id.type="experiment.id")
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
#' getExperiment Generic
#' Generic for getExperiment method
#'
#' @examples
#' data(pdxe)
#' getExperiment(pdxe, experiment.id="X.1004.pael.paclitaxel")
#' getExperiment(pdxe, experiment.id=c("X.1286.pael.paclitaxel", "X.1298.pael.paclitaxel") )
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{list} with the all experiment designs
setGeneric(name = "getExperiment", def = function(object, experiment.id){standardGeneric("getExperiment")} )

#' @export
setMethod( f=getExperiment,
           signature="XevaSet",
           definition=function(object, experiment.id)
           {
             experiment.ids = unique(c(experiment.id))
             rtx = lapply(experiment.ids, .getExperimentDataFromAExpID, object=object)
             rtz = .rbindListOfDataframs(rtx)
             return(rtz)
           })






