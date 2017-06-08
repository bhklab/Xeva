#####================= getMolecularProfiles ==================
#' Get Molecular Profiles
#'
#' Get Molecular Profiles
#'
#' @param object The \code{XevaSet}
#' @param data.type \code{character}, which one of the molecular data types is needed
#' @return a \code{ExpressionSet} where sample names are \code{biobase.id} of model
#' @examples
#' data(pdxe)
#' pdxe_RNA <- getMolecularProfiles(pdxe, data.type="RNASeq")
#' @export
getMolecularProfiles <- function(object, data.type)
{
  if(is.element(data.type, names(slot(object, "molecularProfiles")))==FALSE)
  {
    msg = sprintf("available molecular data are\n%s\n",
                  paste(names(object@molecularProfiles), collapse ="\n"))
    stop(msg)
  }
  expset <- slot(object, "molecularProfiles")[[data.type]]
  return(expset)
}





#####================= Summarize Molecular Profiles ==================
#' summarizeMolecularProfiles
#'
#' summarizeMolecularProfiles
#'
#' @param object The \code{XevaSet}
#' @param drug Name of the drug
#' @param mDataType \code{character}, which one of the molecular data types is needed
#' @param tumor.type default \code{NULL} will return all across all tumor.type
#' @param sensitivity.measure default \code{NULL} will return all sensitivity measure
#' @param unique.model default TRUE will return only one sequncing id, in case where one model id mapes to several sequencing ids
#' @return A \code{ExpressionSet} where sample names are model.id and sensitivity measure will be present in pData
#' @examples
#' data(pdxe)
#' pacRNA <- summarizeMolecularProfiles(pdxe, drug="paclitaxel", mDataType="RNASeq",
#'                                      tumor.type= "BRCA", sensitivity.measure="mRECIST")
#' print(pacRNA)
#' @details
#' \itemize{
#' \item {If a sequencing sample belong to multipal models, summarizeMolecularProfiles
#' will creat saperate column for each model. }
#' \item {All the models without the moleculer data will be removed from the output expression set.}
#' }
#' @export
summarizeMolecularProfiles <- function(object, drug, mDataType, tumor.type=NULL,
                                       sensitivity.measure=NULL, unique.model=TRUE)
{
  modIn <- modelInfo(object, mDataType = mDataType)
  if(unique.model==TRUE)
  {
    modIn <- modIn[unique(modIn$model.id), ]
  }
  modIn <- modIn[modIn$drug %in% c(drug), ]
  if(nrow(modIn)==0)
  {
    msg <- sprintf("No model present with drug %s ", drug)
    stop(msg)
  }

  if(!is.null(tumor.type))
  {
    modIn <- modIn[modIn$tumor.type %in% c(tumor.type), ]
    if(nrow(modIn)==0)
    {
      msg <- sprintf("No model present with tumor.type %s ", tumor.type)
      stop(msg)
    }
  }

  bioName <- sprintf("biobase.id.%s", mDataType)
  modIn <- modIn[ !is.na(modIn[, bioName]), ]
  if(nrow(modIn)==0)
  {
    msg <- sprintf("No model present for drug %s with molecular data type %s", drug, mDataType)
    if(!is.null(tumor.type))
    {
      msg <- sprintf("No model present for drug %s and tumor type %s with molecular data type %s",
                     drug, tumor.type, mDataType)
    }
    stop(msg)
  }

  sm <- sensitivity(object, type = "model", sensitivity.measure = sensitivity.measure)

  if(is.null(sensitivity.measure))
  {
    sensitivity.measure <- colnames(sm)[colnames(sm)!="model.id"]
  }

  modIn[,c(sensitivity.measure)] <- sm[modIn$model.id, c(sensitivity.measure)]
  # if(!is.null(sensitivity.measure))
  # {
  #   sm <- sm[modIn$model.id, ]
  #   if(nrow(sm)==0)
  #   {
  #     msg <- sprintf("No model persent in sensitivity")
  #     stop(msg)
  #   }
  #
  #   for(s in c(sensitivity.measure))
  #   {
  #     modIn[,s] <- sm[modIn$model.id, s]
  #   }
  # }

  molP <- getMolecularProfiles(object, mDataType)
  molP <- molP[, modIn[, bioName]]
  colnames(molP) <- rownames(modIn)
  rownames(Biobase::pData(molP)) <- rownames(modIn)

  pd <- Biobase::pData(molP)
  for(i in colnames(modIn))
  { pd[,i] <- modIn[,i] }
  Biobase::pData(molP) <- pd
  return(molP)
}







