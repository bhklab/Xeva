
.subsetExperimentSlotForDrug <- function(object, drugName, exact.match=TRUE)
{
  dnSplit <- strsplit(drugName, "\\+")[[1]]
  rdx <- data.frame()
  for(Ix in slot(object, "experiment"))
  {
    if(exact.match==TRUE)
    {
      if(drugName==slot(Ix, "drug")[["join.name"]])
      {
        rdx <- rbind(rdx, data.frame(model.id=slot(Ix, "model.id"), drug=drugName,
                                     stringsAsFactors = FALSE))
      }
    }

    if(exact.match==FALSE)
    {
      if(any(dnSplit %in% slot(Ix, "drug")[["names"]])==TRUE)
      {
        rdx <- rbind(rdx, data.frame(model.id=slot(Ix, "model.id"),
                                     drug=slot(Ix, "drug")[["join.name"]],
                                     stringsAsFactors = FALSE))
      }
    }
  }
  return(rdx)
}

##----- select model.id based on drug, tissue -----------------------------------------------------
#' To select model IDs based on drug name and/or tissue type.
#' @examples
#' data(brca)
#' df = selectModelIds(brca, drug="trastuzumab", drug.match.exact=TRUE, tissue="BRCA")
#' head(df)
#' df2 = selectModelIds(brca, drug="trastuzumab", drug.match.exact=FALSE)
#' head(df2)
#'
#' @param object The \code{XevaSet}.
#' @param drug Name of the \code{drug}.
#' @param drug.match.exact Default \code{TRUE}.
#' @param tissue Tumor type. Default \code{NULL}.
#'
#' @return A \code{vector} with the matched \code{model.id}s.
#'
setGeneric(name = "selectModelIds",
           def = function(object,
                          drug=NULL, drug.match.exact=TRUE,
                          tissue=NULL)
           {standardGeneric("selectModelIds")} )

#' @rdname selectModelIds
#' @export
setMethod( f=selectModelIds, signature="XevaSet",
           definition=function(object,
                               drug=NULL, drug.match.exact=TRUE,
                               tissue=NULL)
           {
             if(is.null(drug) & is.null(tissue))
             {stop("drug and tissue both NULL, Please provide atleast one")}

             ExpIdsDrug <- NULL
             if(!is.null(drug))
             {
               ExpIdsDrug <- .subsetExperimentSlotForDrug(object, drug, exact.match=drug.match.exact)
             }

             ExpIdsTumor <- NULL
             if(!is.null(tissue))
             {
               ExpIdsTumor <- mapModelSlotIds(object, id = tissue, id.name = "tissue", map.to="all")
             }

             if(!is.null(drug) & is.null(tissue))
             { return(ExpIdsDrug) }

             if(is.null(drug) & !is.null(tissue))
             { return(ExpIdsTumor) }

             if(!is.null(drug) & !is.null(tissue))
             {
               cmid <- intersect(ExpIdsDrug$model.id, ExpIdsTumor$model.id)
               ExpIdsDrug <- ExpIdsDrug[ExpIdsDrug$model.id %in% cmid,]
               ExpIdsTumor<- ExpIdsTumor[ExpIdsTumor$model.id%in% cmid,]
               rtx <- merge(ExpIdsDrug, ExpIdsTumor, by="model.id")
               return(rtx)
             }
           })



