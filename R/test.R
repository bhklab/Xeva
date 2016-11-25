if(1==2){ print("this is running....")

library(ggplot2)


geoExp = readRDS("DATA-raw/Geo_Exp.Rda")
#experiment = geoExp$experiment
#model = geoExp$model

pdxe = creatPharmacoPXSet(name = "PDXE",
                          molecularProfiles = list(RNASeq = geoExp$RNASeq),
                          experiment = geoExp$experiment,
                          model = geoExp$model,
                          drug  = geoExp$drug
                          )
save(pdxe, file = "data/pdxe.rda")

data("pdxe")



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



##================================================================

  setClass(
    Class = "Circle",
    representation = representation(
      radius = "numeric",
      diameter = "numeric"
    ) )

  a <- new("Circle")
  a@radius = 1
  radius(a) ##gives error could not find function "radius"
  radius(a) <- 2 ##gives error could not find function "radius<-"

  # define radius
  setGeneric("radius", function(self) standardGeneric("radius"))
  setMethod("radius",
            signature(self = "Circle"),
            function(self) {
              self@radius
            } )

  radius(a) ##this works
  radius(a) <- 2  ##gives error could not find function "radius<-"

  setGeneric("radius<-", function(self, value) standardGeneric("radius<-"))
  setReplaceMethod("radius",
                   "Circle",
                   function(self, value) {
                     self@radius <- value
                     self
                   } )


  radius(a) <- 2  ##Now this works
  radius(a) ##this will give 2



  ## to get the diameter
  setGeneric("diameter", function(self) standardGeneric("diameter"))
  setMethod("diameter",
            signature(self = "Circle"),
            function(self) {
              self@diameter
            }
  )

  ## to set the diameter
  setGeneric("diameter<-", function(self, value) standardGeneric("diameter<-"))
  setReplaceMethod("diameter", "Circle",
                   function(self, value) {
                     self@diameter <- value
                     self
                   } )


  # Method that calculates one value from another
  setGeneric("calc_diameter", function(self) { standardGeneric("calc_diameter")})
  setMethod("calc_diameter",
            signature(self = "Circle"),
            function(self) {
              self@diameter <- self@radius * 2
              self
            }
  )






}


