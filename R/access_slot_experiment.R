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
  rdx = data.frame(model.id=unique(modList), drug=drugName)
  return(rdx)
}



#' #' get experiment time series data for a given biobase.id or model.id
#' #'
#' #' @examples
#' #' data(pdxe)
#' #' # extract controls for a given model.id
#' #' getExperimentData(object=pdxe, id ="X.007.trab", id.type="model.id")
#' #' getExperimentData(object=pdxe, id ="X-007", id.type="biobase.id")
#' #' @param object The \code{Xeva} dataset
#' #' @param model.id The \code{model.id}
#' #' @return a \code{vector} with control model.id
#' setGeneric(name = "getExperimentData", def = function(object, id, id.name, exact.match) {standardGeneric("getExperimentData")})
#'
#' #' @export
#' setMethod( f=getExperimentData,
#'            signature=c(object="XevaSet"),
#'            definition= function(object, id, id.name, exact.match=TRUE)
#'            {
#'
#'              if(id.name=="drug")
#'              { modId = .subsetExperimentSlotForDrug(object, drugName, exact.match=TRUE) }
#'
#'              if(id.name!="model.id")
#'              {
#'                if(is.element(id.name, colnames(object@model))==TRUE)
#'                {
#'                  #mdf = getModelIds(object, id=id, id.name=id.name)
#'                  #mapModelSlotIds(object, id, id.name, map.to)
#'                }
#'              }
#'
#'
#'              lRt = list()
#'              for(I in 1:dim(mdf)[1])
#'              {
#'                model.id= mdf[I, "model.id"]
#'                map.id  = mdf[I, id.type]
#'
#'                expSubSet = .subsetExperimentSlot(object, id=model.id, id.type="model.id")
#'                if(length(expSubSet)>0)
#'                {
#'                  dfy = do.call(rbind, lapply(expSubSet, "[[", "data"))
#'                  dfy[, "model.id"] = model.id
#'                  dfy[, id.name] = map.id
#'                  lRt = .appendToList(lRt, dfy)
#'                }
#'              }
#'              dfRt = .rbindListOfDataframs(lRt)
#'              dfRt = dfRt[, !apply(is.na(dfRt), 2, all)]
#'              return(dfRt)
#'            })

##----------------------------------------------------



#' get experiment time series data for any id (biobase.id, model.id or tumor.type)
#'
#' @examples
#' data(pdxe)
#' # extract controls for a given model.id
#' GCdf = getModelTimeSeries(object=pdxe, id ="GC", id.name ="tumor.type")
#' GCdf = getModelTimeSeries(object=pdxe, id ="GC", id.name ="tumor.type",
#'                           drug="binimetinib", drug.match.exact=TRUE, include=NULL)
#' GCdf = getModelTimeSeries(object=pdxe, id ="GC", id.name ="tumor.type",
#'                           drug="binimetinib", drug.match.exact=TRUE,
#'                           include=c("mRECIST"))
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return a \code{vector} with control model.id
getModelTimeSeries <- function(object,id=NULL, id.name=NULL,
                       drug=NULL, drug.match.exact=TRUE, include=NULL)
{

  include= c(include, id.name)

  IDfrmMod = NULL
  if(!is.null(id))
  {
    if(is.null(id.name)){stop("For id id.name must be specifyied")}
    IDfrmMod = mapModelSlotIds(object, id = id, id.name = id.name, map.to="all")
  }

  IDfrmDrug = NULL
  if(!is.null(drug))
  {
    IDfrmDrug = .subsetExperimentSlotForDrug(object, drug, exact.match=drug.match.exact)
  }

  ##---- get the intersect of model.id ---------------------------------
  if(!is.null(IDfrmMod) & !is.null(IDfrmDrug))
  {
    mid = intersect(IDfrmMod$model.id, IDfrmDrug$model.id)
  } else
  {
    if(!is.null(IDfrmMod) ){ mid = IDfrmMod$model.id }
    if(!is.null(IDfrmDrug)){ mid = IDfrmDrug$model.id}
  }

  if(!is.null(include))
  {
    if(is.vector(include))
    {
      includeLt = as.list(include)
      names(includeLt) = include
    }

    ##----for nested values creat list -----------------
    #if(grep("^best.response$", includeLt))
  } else{includeLt=NULL}

  rtx = list()
  for(mx in mid)
  {
    expSubSetX= .subsetExperimentSlot(object, id=mx, id.type="model.id")
    expSubSet = expSubSetX[[1]]
    dfx = expSubSet[["data"]]
    dfx[, "model.id"] = mx
    dfx[, "drug.name"] = expSubSet[["drug"]][["join.name"]]
    if(!is.null(includeLt))
    {
      ##--- change for nasted values ---------------
      for(incl in includeLt)
      {
        dfx[, incl] =  expSubSet[[incl]]
      }
    }
    rtx = .appendToList(rtx, dfx)
  }

  dfRt = .rbindListOfDataframs(rtx)

  if(!is.null(includeLt))
  {
    for(incl in includeLt)
    {
      if(is.element(incl, colnames(dfRt))==FALSE)
      {
        if(is.element(incl, colnames(object@model))==TRUE)
        {
          clx = mapModelSlotIds(object, id = dfRt$model.id,
                                id.name = "model.id", map.to=incl)
          dfRt = merge(dfRt, clx, by.x = "model.id", by.y = "model.id")
        }
      }
    }
  }

  ##---- remove NA cols ----------------------
  dfRt = dfRt[, !apply(is.na(dfRt), 2, all)]

  return(dfRt)
}

















