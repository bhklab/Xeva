
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



################################################################################
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

