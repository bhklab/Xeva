processRawData <- function()
{
source("../PT/src/cxpMeta_mysql_functions.R")
DBcon = getCXPMetaDB_Conn()
qtxt = "SELECT * FROM PDX_Gao2015_NatureMed_RAW;"
rs <- dbSendQuery(DBcon, qtxt)

experimentX <- fetch(rs, n=-1)
closeAllDBconn()



##------------------------------------------------
modTr = unique(experimentX[, c("Model", "Treatment.1", "Treatment.2", "Treatment.3", "Treatment")])

library(stringi)
creatNewModId <- function(modId, Tx)
{
  Tx = c(Tx)
  Tx = Tx[!is.na(Tx)]
  drgName = paste(sapply(Tx, function(x)
                            { S1 = stri_sub(x,1,2); S2 = stri_sub(x,-2,-1)
                              sprintf("%s%s", S1, S2)}
                         ), collapse = ".")

  modIdx = sprintf("%s.%s", modId, drgName)
  return(modIdx)
}

model.id = apply(modTr, 1, function(x){
                           Tx = c(x["Treatment.1"], x["Treatment.2"], x["Treatment.3"])
                           creatNewModId(x["Model"], Tx) })

#model.id[duplicated(model.id)]
model.id = make.names(model.id, unique=TRUE)
modTr$model.id = model.id



experimentX$model.id = apply(experimentX, 1, function(x){ modTr[modTr$Model== x["Model"] &
                                                                  modTr$Treatment== x["Treatment"], "model.id"] })

##-------------------------------------------------------------------------
#modTr = unique(experimentX[, c("Model", "Treatment")])
#mid = lapply( unique(modTr$Model), function(x){
#                                   mx = modTr[modTr$Model==x, "Model"]
#                                   paste(mx, 1:length(mx), sep=".")})
#
#modTr$model.id = unlist(mid)
#experimentX$model.id = apply(experimentX, 1, function(x){ modTr[modTr$Model== x["Model"] &
#                                                                modTr$Treatment== x["Treatment"], "model.id"] })


experiment = data.frame(model.id = experimentX$model.id,
                        drug.1 = experimentX$Treatment.1,
                        drug.2 = experimentX$Treatment.2,
                        drug.3 = experimentX$Treatment.3,
                        time  = experimentX$Days.Post.T0,
                        volume= experimentX$`Volume.(mm3)`,
                        body.weight= experimentX$`body.weight.(g)`,
                        stringsAsFactors = FALSE)

tumor.type = experimentX$Tumor.Type
tumor.type[tumor.type=="BRCA"]= "Breast Cancer" #BRCA, breast carcinoma
tumor.type[tumor.type=="CM"]  = "Cutaneous Melanoma" #CM, cutaneous melanoma
tumor.type[tumor.type=="CRC"] = "Colorectal Cancer"  #CRC, colorectal cancer
tumor.type[tumor.type=="GC"]  = "Gastric Cancer"     #GC, Gastric Cancer
tumor.type[tumor.type=="NSCLC"]="Non-small Cell Lung Carcinoma" #NSCLC, non-small cell lung carcinoma;
tumor.type[tumor.type=="PDAC" ]="Pancreatic Ductal Carcinoma"   #PDAC, pancreatic ductal carcinoma;

experiment$tumor.type = as.character(tumor.type)

#experiment$batch  = sapply(strsplit(experiment$model.id, "[.]"), `[[`, 1)
experiment$batch  = experimentX$Model

exp.type = experiment$drug.1
exp.type[exp.type=="untreated"] = "control"
exp.type[exp.type!="control"  ] = "treatment"
experiment$exp.type = exp.type

##====================================
####----- create model matrix ---------------

model = unique(experiment[, c("model.id", "batch")])
colnames(model) = c("model.id", "biobase.id")
model$patient.id = model$biobase.id
#exp.type = experimentX$Treatment.1
#exp.type[exp.type=="untreated"] = "control"
#exp.type[exp.type!="control"  ] = "treatment"
#edf = experimentX[, c("model.id", "Model")]
#colnames(edf) = c("model.id", "batch")
#edf$exp.type = exp.type
#model = unique(edf)

#seqObjId  = sapply(strsplit(model$model.id, "[.]"), `[[`, 1)
#model$biobase.id = seqObjId
#model$biobase.id = experiment$batch

geoExp = list(experiment=experiment, model = model)

##====== get the drug datafram =============================================

drgNames = apply(geoExp$experiment[, c("drug.1","drug.2","drug.3")], 1, Xeva::pasteWithoutNA)
drgNames = unique(drgNames)
drug = data.frame(drug.id = drgNames,
                 standard.name = drgNames)

geoExp$drug = drug


##====================================
## Create design matrix ---------------------------------

#head(experiment)
experiment$drug = Xeva::pasteColTogather(experiment[, c("drug.1","drug.2","drug.3")], collapse = " + ")

dsg = unique(experiment[, c("model.id","drug","batch","exp.type")])

expDesign = list()
for(b in unique(dsg$batch))
{
  dx = dsg[dsg$batch==b, ]
  ##--- for a batch all ids with "control" will be considered as control
  controlIDs = c(unique( dx[dx$exp.type=="control", "model.id"]))
  if(length(controlIDs)==0){controlIDs=c(NULL)}
  tretDrugs = unique( dx[dx$exp.type=="treatment", "drug"] )
  for(drugI in tretDrugs)
  {
    #v = list(batch=b, drug=drugI)
    v = list(batch.name = sprintf("%s.%s", b, drugI))
    v$treatment = c(dx[dx$drug==drugI & dx$exp.type=="treatment", "model.id"])
    v$control = controlIDs
    expDesign[[length(expDesign)+1]] = v
  }
}

geoExp$expDesign = expDesign

##------------------------------------------------------
#library(Biobase)
rseq = readRDS("DATA-raw/Geo_RNAseq_fpkm.rdata")
assaydata = as.matrix(rseq[1:50,]) ##only 50 genes

sampleID = colnames(assaydata)
rnaseqMeta = data.frame(sampleID=sampleID)
rownames(rnaseqMeta) = rnaseqMeta$sampleID
phenodata   <- Biobase::AnnotatedDataFrame(data = rnaseqMeta)
#phenodata   <- new("AnnotatedDataFrame", data = rnaseqMeta)

featureDF = data.frame(geneName = rownames(assaydata), ensembl.id = NA)
rownames(featureDF) = featureDF$geneName
#featuredata <- new("AnnotatedDataFrame", data = featureDF)
featuredata <- Biobase::AnnotatedDataFrame(data = featureDF)


rnaseq <- Biobase::ExpressionSet(assayData=assaydata,
                                 phenoData=phenodata,
                                 featureData=featuredata)



geoExp$RNASeq = rnaseq

saveRDS(geoExp, file = "DATA-raw/Geo_Exp.Rda")


}


creatXevaObject <- function()
{
  processRawData()
  geoExp = readRDS("DATA-raw/Geo_Exp.Rda")
  pdxe = creatXevaSet(name = "PDXE",
               molecularProfiles = list(RNASeq = geoExp$RNASeq),
               experiment = geoExp$experiment,
               expDesign  = geoExp$expDesign,
               model = geoExp$model,
               drug  = geoExp$drug)

  setmRECIST(pdxe)<- setmRECIST(pdxe)
  save(pdxe, file = "data/pdxe.rda")

  data(pdxe)
}



#creatXevaObject()
