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
  #treatment= model[model$batch==batchID & model$exp.type=="treatment", "model.id"]
  #control  = model[model$batch==batchID & model$exp.type=="control"  , "model.id"]
  #return(list(treatment=model.idx, control=control))
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



getExpDesign <- function(objNames, expDesign)
{
  expDesIndx = sapply(expDesign, function(x){
                      if( length(intersect(objNames, c(x$treatment,x$control) ))>0)
                      {return(TRUE)}
                      return(FALSE) })
  expDesName = names(expDesign)[expDesIndx]
  return(expDesign[expDesName])
}




gPx = readRDS("data/Gao_PharPx_obj.Rda")
expSlot  = gPx$experiment
model = gPx$model
expDesign = gPx$expDesign






drug = c( "LCL161", "paclitaxel")
#drug = "paclitaxel"
objNames = getModels(expSlot, drug=drug, drug.exact.match = TRUE)
length(objNames)

unique(sapply(expSlot[objNames], "[[" , "tumor.type" ))

expDesX = getExpDesign(objNames, expDesign)

tosave = list(expSlot=expSlot, expDesX=expDesX)
saveRDS(tosave, file = "data/toPlot.Rda")



drug = c( "BKM120", "binimetinib")
tumor.type = "Colorectal Cancer"
objNames = getModels(expSlot, drug=drug, drug.exact.match = TRUE, tumor.type=tumor.type)
length(objNames)

unique(sapply(expSlot[objNames], "[[" , "tumor.type" ))
unique(sapply(expSlot[objNames], "[[" , c("drug", "join.name")))

objPDX = getControlExp(expSlot, objNames, model)










#tumor.type = "Breast Cancer"
tumor.type = "Non-small Cell Lung Carcinoma"
z = getModels(expSlot, drug=drug, drug.exact.match = FALSE, tumor.type=tumor.type)
length(z)
unique(sapply(z , "[[" , "tumor.type" ))
unique(lapply(z, `[[`, c("drug", "join.name")))

#sapply(z , "[[" , "exp.type" )
#k=lapply(expSlot, `[[`, c("drug", "join.name"))


drug = "paclitaxel"
tumor.type = "Breast Cancer"
#tumor.type = "Non-small Cell Lung Carcinoma"
z = getModels(expSlot, drug=drug, drug.exact.match = TRUE, tumor.type=tumor.type)

length(z)
unique(sapply(z , "[[" , "tumor.type" ))
unique(lapply(z, `[[`, c("drug", "join.name")))


drug = "encorafenib"
tumor.type = "Cutaneous Melanoma"
z = getModels(expSlot, drug=drug, drug.exact.match = TRUE, tumor.type=tumor.type)
length(z)


objPDX = z[[1]]
expDesign = bIDdf
getControlExp(objPDX, expDesign)






colors <- rainbow(length(z));
plot(c(1,100), c(1,1000),type="n", xlab="Time", ylab="Tumor volume" )
colIndx = 1
for(i in z)
{
  lines(i$data$time, i$data$volume, type = "o", pch=19, col=colors[colIndx])
  colIndx = colIndx+1
  #print(i$drug)
}

drgName = sapply(z, function(x)x$drug$join.name)

legend("topleft", legend=drgName, cex=0.8, col=colors, pch=19, lty=1, title="Drug")


















