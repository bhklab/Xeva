##---- functions to acess experiment slot --------

drugMatchFun <- function(drug, objPDX, exact.match)
{
  objDrugNames = toupper(objPDX$drug$names)
  if(exact.match == TRUE &
     length(intersect(drug, objDrugNames))== length(objDrugNames) &
     length(intersect(drug, objDrugNames))== length(drug) )
  { return(TRUE) }

  if(exact.match == FALSE & length(intersect(drug, objDrugNames))>0 )
  { return(TRUE) }

  return(FALSE)
}

tumorTypeMatchFun<- function(tumor.type, objPDX)
{

  if (toupper(tumor.type) == toupper(objPDX$tumor.type))
  {return(TRUE)}

  return(FALSE)
}


#' @export
getModels <- function(expSlot, drug=NULL, drug.exact.match=TRUE, tumor.type=NULL)
{
  objIndx = list()
  if(!(is.null(drug)))
  {
    drug = c(toupper(drug))
    objIndx[["Drug"]] = sapply(expSlot, drugMatchFun, drug=drug, exact.match=drug.exact.match)
  }

  if(!(is.null(tumor.type)))
  {
    objIndx[["tumor.type"]] = sapply(expSlot, tumorTypeMatchFun, tumor.type=tumor.type)
  }

  rtIndx = apply( do.call(cbind.data.frame, objIndx), 1, all )
  #rtName = sapply(expSlot, "[[", "model.id")[rtIndx]
  rtName = names(expSlot)[rtIndx]

  #return(expSlot[rtIndx])
  return(rtName)
}


getTreatmentControlForModel <- function(model.idx)
{
  batchID = model[model$model.id ==model.idx, "batch"]
  return(model[model$batch==batchID, c("model.id", "batch", "exp.type")])
}


getTreatmentControlX <- function(expSlot, objNames, model)
{
  drgNames  = unique(sapply(expSlot[objNames], "[[", c("drug", "join.name")))
  tretBatch = unique(model[model$drug.join.name%in% drgNames,  "batch"])

  rtx=list()
  for(drgI in drgNames)
  {
    for(batI in tretBatch)
    {
      Lx = list(drug.join.name = drgI,
                batch = batI)
      Lx$treatment = unique(model[model$drug.join.name == drgI &
                        model$batch == batI &
                        model$exp.type == "treatment", "model.id"] )

      Lx$control = unique(model[model$batch == batI & model$exp.type == "control", "model.id"] )

      namx = sprintf("%s.%s", drgI, batI)
      rtx[[namx]]= Lx
    }
  }


  rdf = data.frame()
  for(model.idx in objNames)
  {
    tc = getTreatmentControlForModel(model.idx)
    tc$drug.join.name = sapply(expSlot[tc$model.id], "[[", c("drug", "join.name"))
    rdf = rbind(rdf, tc)
  }

  rdf = unique(rdf)
  #tretID = rdf[rdf$exp.type=="treatment", "model.id"]
  drgeID = unique(rdf[rdf$exp.type=="treatment", "drug.join.name"])
  rtx = list()
  for(drI in drgeID)
  {
    rdf[rdf$drug.join.name==drI, ]
  }


  return(rtx)
}


#' @export
getExpDesign <- function(objNames, expDesign)
{
  expDesIndx = sapply(expDesign, function(x){
                      if( length(intersect(objNames, c(x$treatment,x$control) ))>0)
                      {return(TRUE)}
                      return(FALSE) })
  expDesName = names(expDesign)[expDesIndx]
  return(expDesign[expDesName])
}



getExperimentValues <- function(object, name = c("experiment.id"))
{
  #object = pdxe
  #name = c("drug", "join.name")
  rt  = sapply(object@experiment, "[[", c(name))
  return(rt)
}









