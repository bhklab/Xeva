if(1==2){


source("../PT/src/cxpMeta_mysql_functions.R")
DBcon = getCXPMetaDB_Conn()
qtxt = "SELECT * FROM PDX_Gao2015_NatureMed_RAW;"
rs <- dbSendQuery(DBcon, qtxt)

experimentX <- fetch(rs, n=-1)
closeAllDBconn()

modTr = unique(experimentX[, c("Model", "Treatment")])

mid = lapply( unique(modTr$Model), function(x){
                                   mx = modTr[modTr$Model==x, "Model"]
                                   paste(mx, 1:length(mx), sep=".")})

modTr$model.id = unlist(mid)

experimentX$model.id = apply(experimentX, 1, function(x){ modTr[modTr$Model== x["Model"] &
                                                                modTr$Treatment== x["Treatment"], "model.id"] })

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

experiment$batch  = sapply(strsplit(experiment$model.id, "[.]"), `[[`, 1)

exp.type = experiment$drug.1
exp.type[exp.type=="untreated"] = "control"
exp.type[exp.type!="control"  ] = "treatment"
experiment$exp.type = exp.type

##====================================
## Create design matrix -----
####----- create model matrix ---------------

exp.type = experimentX$Treatment.1
exp.type[exp.type=="untreated"] = "control"
exp.type[exp.type!="control"  ] = "treatment"
edf = experimentX[, c("model.id", "Model")]
colnames(edf) = c("model.id", "batch")
edf$exp.type = exp.type
model = unique(edf)

seqObjId  = sapply(strsplit(model$model.id, "[.]"), `[[`, 1)
model$biobase.id = seqObjId

geoExp = list(experiment=experiment, model = model)


###=========================================================================
##========= creat experiment design list ===================================
#expDesignDf = data.frame(model.id = "",
#                         batch = "",
#                         exp.type = "",
#                         drug = ""
#                         )
##==========================================================================
##====== get the drug datafram =============================================

pasteWithoutNA <- function(L, collapse = " + "){paste(L[!is.na(L)], collapse = collapse)}
drgNames = apply(geoExp$experiment[, c("drug.1","drug.2","drug.3")], 1, pasteWithoutNA)
drgNames = unique(drgNames)
drug = data.frame(drug.id = drgNames,
                 standard.name = drgNames)

geoExp$drug = drug

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


rnaseq <- ExpressionSet(assayData=assaydata,
                      phenoData=phenodata,
                      featureData=featuredata)



geoExp$RNASeq = rnaseq



saveRDS(geoExp, file = "DATA-raw/Geo_Exp.Rda")

#geoExp = readRDS("DATA/Geo_Exp.Rda")

}