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
  modList = c()
  for(Ix in object@experiment)
  {
    if(exact.match==TRUE)
    {
      if(length(drgNames) == length(Ix$drug$names))
      {
        if(all(drgNames %in% Ix$drug$names)==TRUE)
        {modList = c(modList, Ix$model.id)}
      }
    }
    ###------------------------------------------------
    if(exact.match==FALSE)
    {
      if(any(drgNames %in% Ix$drug$names)==TRUE)
      {modList = c(modList, Ix$model.id)}
    }
    ###------------------------------------------------
  }
  rdx = data.frame(model.id=unique(modList), drug=drugName, stringsAsFactors = FALSE)
  return(rdx)
}

##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
##----- get model Ids -----------------------------
#' getModelIds Generic
#' Generic for getModelIds method
#'
#' @examples
#' data(pdxe)
#' getModelIds(pdxe, drug="paclitaxel", drug.match.exact=TRUE, tumor.type="BRCA")
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{list} with the all experiment designs
setGeneric(name = "getModelIds",
           def = function(object,
                          drug=NULL, drug.match.exact=TRUE,
                          tumor.type=NULL)
                  {standardGeneric("getModelIds")} )

#' @export
setMethod( f=getModelIds, signature="XevaSet",
           definition=function(object,
                               drug=NULL, drug.match.exact=TRUE,
                               tumor.type=NULL)
{
  if(is.null(drug) & is.null(tumor.type))
  {stop("drug and tumor.type both NULL, Please provide atleast one")}

  ModIdsDrug = NULL
  if(!is.null(drug))
  {
    drug = c(drug)
    ModIdsDrug = .subsetExperimentSlotForDrug(object, drug, exact.match=drug.match.exact)
  }

  ModIdsTumor = NULL
  if(!is.null(tumor.type))
  {
    tumor.type = c(tumor.type)
    ModIdsTumor = mapModelSlotIds(object, id = tumor.type, id.name = "tumor.type", map.to="all")
  }


  if(!is.null(drug) & is.null(tumor.type))
  { return(ModIdsDrug$model.id) }

  if(is.null(drug) & !is.null(tumor.type))
  { return(ModIdsTumor$model.id) }

  if(!is.null(drug) & !is.null(tumor.type))
  {
    rtx = intersect(ModIdsDrug$model.id, ModIdsTumor$model.id)
    return(rtx)
  }
})



###-----------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------
.getExperimentDataFromAModelID <- function(object, model.id)
{
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
#' getExperiment Generic
#' Generic for getExperiment method
#'
#' @examples
#' data(pdxe)
#' getExperiment(pdxe, model.id="X.1655.LE11.biib")
#' getExperiment(pdxe, model.id=c("X.1655.LE11.biib", "X.1298.pael") )
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{list} with the all experiment designs
setGeneric(name = "getExperiment", def = function(object, model.id){standardGeneric("getExperiment")} )

#' @export
setMethod( f=getExperiment, signature="XevaSet",
           definition=function(object, model.id)
           {
             model.ids = unique(c(model.id))
             rtx = lapply(model.ids, .getExperimentDataFromAModelID, object=object)
             rtz = .rbindListOfDataframs(rtx)
             return(rtz)
           })








