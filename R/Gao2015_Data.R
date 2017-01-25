read_curveMetrics <- function()
{

  rd = readRDS("~/CXP/XG/Data/Gao_2015_NatureMed/nm.3954-S2_PCT_raw_data.Rda")
  rdv = unique(rd[, c("Model","Tumor.Type")])
  #naModID = dzv[is.na(dzv$Tumor.Type), "Model"]
  #dzv[dzv$Model%in%naModID,]
  ##--- remove NA
  rdv = rdv[!is.na(rdv$Tumor.Type),]
  ##--- read and merge the curve matrix ---------------------------------------
  fl = "~/CXP/XG/Data/Gao_2015_NatureMed/nm.3954-S2_PCT_curve_metrics.Rda"
  cvm = readRDS(fl)

  cvm = merge(cvm, rdv, by.x = "Model", by.y = "Model")
  pubLung = unique(cvm[cvm$Tumor.Type=="NSCLC", c("Model","Treatment")])
  table(pubLung$Treatment)

  data(pdxe)
  dfx = getmRECIST(pdxe)
  dfMap = mapModelSlotIds(object=pdxe, id=dfx$model.id, id.name="model.id",
                          map.to="tumor.type", unique=TRUE)

  dfx = merge(dfx, dfMap, by.x = "model.id", by.y = "model.id")
  lungDf = dfx[dfx$tumor.type=="NSCLC",]

  pubLung[!(pubLung$Model %in% lungDf$biobase.id),]

  lungDf[!(lungDf$biobase.id %in% pubLung$Model), "biobase.id"]






}



processRawData <- function()
{
source("~/CXP/PT/src/cxpMeta_mysql_functions.R")
DBcon = getCXPMetaDB_Conn()
qtxt = "SELECT * FROM PDX_Gao2015_NatureMed_RAW;"
rs <- dbSendQuery(DBcon, qtxt)

experimentX <- fetch(rs, n=-1)
closeAllDBconn()

##-----------------------------------------------------------------------

curMFile <- "~/CXP/XG/Data/Gao_2015_NatureMed/nm.3954-S2_PCT_curve_metrics.Rda"
curM <- readRDS(curMFile)
curM$Treatment = gsub('figitumumab\"', 'figitumumab', curM$Treatment)

curM$mod.treat = paste(curM$Model, curM$Treatment, sep = "_")
experimentX$mod.treat = paste(experimentX$Model, experimentX$Treatment, sep = "_")

expMerge = merge(experimentX, curM, by.x = "mod.treat", by.y = "mod.treat")

colNm <- c("mod.treat", "Model.x", "Tumor.Type", "Treatment.x", "Volume.(mm3)",
           "body.weight.(g)", "Days.Post.T0", "percent.TVol.Difference",
           "percent.BW.Difference", "N.Treatment", "Treatment.1", "Treatment.2",
           "Treatment.3", "Treatment.target", "Treatment.type", "BestResponse",
           "Day_BestResponse", "BestAvgResponse", "Day_BestAvgResponse",
           "TimeToDouble", "Day_Last", "ResponseCategory")

expMerge = expMerge[, colNm]
colnames(expMerge) = gsub("Model.x", "Model", colnames(expMerge))
colnames(expMerge) = gsub("Treatment.x", "Treatment", colnames(expMerge))

experimentX = expMerge
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

experiment = data.frame(model.id = experimentX$model.id,
                        drug.1 = experimentX$Treatment.1,
                        drug.2 = experimentX$Treatment.2,
                        drug.3 = experimentX$Treatment.3,
                        time  = experimentX$Days.Post.T0,
                        volume= experimentX$`Volume.(mm3)`,
                        body.weight= experimentX$`body.weight.(g)`,
                        treatment.target = experimentX$Treatment.target,
                        treatment.type = experimentX$Treatment.type,
                        best.response_published = experimentX$BestResponse,
                        time.best.response_published = experimentX$Day_BestResponse,
                        best.avg.response_published = experimentX$BestAvgResponse,
                        time.best.avg.response_published = experimentX$Day_BestAvgResponse,
                        timeToDouble_published = experimentX$TimeToDouble,
                        time.last_published = experimentX$Day_Last,
                        mRECIST_published = experimentX$ResponseCategory,
                        stringsAsFactors = FALSE)

experimentMin = experiment
geoExp = list(experiment=experimentMin)

##---------------------------------------------------------------------------------------------------------
##=========================================================================================================
####----- create model matrix ------------------------
tumor.type = experimentX$Tumor.Type
tumor.type.name = tumor.type
tumor.type.name[tumor.type=="BRCA"]= "Breast Cancer" #BRCA, breast carcinoma
tumor.type.name[tumor.type=="CM"]  = "Cutaneous Melanoma" #CM, cutaneous melanoma
tumor.type.name[tumor.type=="CRC"] = "Colorectal Cancer"  #CRC, colorectal cancer
tumor.type.name[tumor.type=="GC"]  = "Gastric Cancer"     #GC, Gastric Cancer
tumor.type.name[tumor.type=="NSCLC"]="Non-small Cell Lung Carcinoma" #NSCLC, non-small cell lung carcinoma;
tumor.type.name[tumor.type=="PDAC" ]="Pancreatic Ductal Carcinoma"   #PDAC, pancreatic ductal carcinoma;
experiment$tumor.type = as.character(tumor.type)
experiment$tumor.type.name = as.character(tumor.type.name)
experiment$batch  = experimentX$Model
exp.type = experiment$drug.1
exp.type[exp.type=="untreated"] = "control"
exp.type[exp.type!="control"  ] = "treatment"
experiment$exp.type = exp.type

model = unique(experiment[, c("model.id", "batch","tumor.type", "tumor.type.name")])
colnames(model) = c("model.id", "biobase.id", "tumor.type", "tumor.type.name")
model$patient.id = model$biobase.id

geoExp$model = model

##====== get the drug datafram ===================================================================

drgNames = apply(geoExp$experiment[, c("drug.1","drug.2","drug.3")], 1, Xeva::pasteWithoutNA)
drgNames = unique(drgNames)
drug = data.frame(drug.id = drgNames,
                 standard.name = drgNames)

geoExp$drug = drug


##=================================================================================================
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
assaydata = as.matrix(rseq)#[1:50,]) ##only 50 genes
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

##---------------------------------------------------------

##--- do once ---------------------------------------------
if(1==2){
mutX = readRDS("DATA-raw/Geo_Mut.rdata")
rNam = unique(mutX$Gene); cNam = unique(mutX$Sample)
mut = data.frame(matrix(0, nrow = length(rNam), ncol = length(cNam)))
rownames(mut) = rNam ; colnames(mut)= cNam

for(I in 1:nrow(mutX))
{
  vx = mutX[I, ]
  if(mut[vx$Gene, vx$Sample]==0){ px = vx$Category}
  if(mut[vx$Gene, vx$Sample]!=0)
  {
    pOld = mut[vx$Gene, vx$Sample]
    pNew = vx$Category
    p1 = c(pNew, unlist(strsplit(pOld, ";")))
    px = paste(unique(p1), collapse = ";")
  }
  mut[vx$Gene, vx$Sample] = px
}

geneMap = unique( mutX[, c("Gene", "Entrez")] )
##rownames(geneMap) = geneMap$Gene
saveRDS(list(mut=mut, geneMap=geneMap), file="DATA-raw/Geo_Mut_matrix.rdata")
}

mutLst = readRDS("DATA-raw/Geo_Mut_matrix.rdata")
assaydata = as.matrix(mutLst$mut)
sampleID = colnames(assaydata)
mutMeta = data.frame(sampleID=sampleID)
rownames(mutMeta) = mutMeta$sampleID
phenodata   <- Biobase::AnnotatedDataFrame(data = mutMeta)


featureDF = data.frame(geneName = rownames(assaydata), Entrez=NA)
rownames(featureDF) = featureDF$geneName
#featuredata <- new("AnnotatedDataFrame", data = featureDF)
featuredata <- Biobase::AnnotatedDataFrame(data = featureDF)


mutSeq <- Biobase::ExpressionSet(assayData=assaydata,
                                 phenoData=phenodata,
                                 featureData=featuredata)
geoExp$mutation = mutSeq

##---------------------------------------------------------


saveRDS(geoExp, file = "DATA-raw/Geo_Exp.Rda")

}


creatXevaObject <- function()
{
  processRawData()
  geoExp = readRDS("DATA-raw/Geo_Exp.Rda")
  library(Xeva)
  pdxe = creatXevaSet(name = "PDXE",
               molecularProfiles = list(RNASeq = geoExp$RNASeq, mutation=geoExp$mutation ),
               experiment = geoExp$experiment,
               expDesign  = geoExp$expDesign,
               model = geoExp$model,
               drug  = geoExp$drug)

  setmRECIST(pdxe)<- setmRECIST(pdxe)
  setSlop(pdxe) <- setSlop(pdxe, treatment.only=FALSE)
  setAngle(pdxe) <- setAngle(pdxe)

  save(pdxe, file = "data/pdxe.rda")

  data(pdxe)
}

#creatXevaObject()
