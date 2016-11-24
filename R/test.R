if(1==2){ print("this is running....")

library(ggplot2)


geoExp = readRDS("DATA/Geo_Exp.Rda")
experiment = geoExp$experiment
model = geoExp$model

#expSlot = experimentSlotfromDf(experiment)

#model = checkModel(model, expSlot)

#expDesign = creatExperimentDesign(model, expSlot)



pdxe = creatPharmacoPXSet("PDXE",
                          molecularProfiles = list(RNASeq = geoExp$RNASeq),
                          experiment = geoExp$experiment,
                          model = geoExp$model,
                          drug  = geoExp$drug
                          )









if(1==2){



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





}















}

