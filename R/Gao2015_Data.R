rseq = readRDS("data/Geo_RNAseq_fpkm.Rda")
#crmt = readRDS("data/Geo_curve_metrics.Rda")
#experiment = readRDS("data/Geo_raw_data.Rda")

source("../PT/src/cxpMeta_mysql_functions.R")
DBcon = getCXPMetaDB_Conn()
qtxt = "SELECT * FROM PDX_Gao2015_NatureMed_RAW;"
rs <- dbSendQuery(DBcon, qtxt)

experimentX <- fetch(rs, n=-1)
closeAllDBconn()

experiment = data.frame(model.id = experimentX$Model,
                        drug1 = experimentX$Treatment.1,
                        drug2 = experimentX$Treatment.2,
                        drug3 = experimentX$Treatment.3,
                        time  = experimentX$Days.Post.T0,
                        volume= experimentX$`Volume.(mm3)`,
                        body.weight = experimentX$`body.weight.(g)`)



pdxModels = unique(experiment$Model)


library(Biobase)
assaydata = as.matrix(rseq[1:50,]) ##only 50 genes

sampleID = colnames(assaydata)


models = sapply(sampleID, function(x){ if(x%in%pdxModels) x else NA })

rnaseqMeta = data.frame(sampleID=sampleID, model.id=models )
rownames(rnaseqMeta) = rnaseqMeta$sampleID
phenodata   <- new("AnnotatedDataFrame", data = rnaseqMeta)

featureDF = data.frame(geneName = rownames(assaydata), ensembl.id = NA)
rownames(featureDF) = featureDF$geneName
featuredata <- new("AnnotatedDataFrame", data = featureDF)

eset <- ExpressionSet(assayData=assaydata,
                      phenoData=phenodata,
                      featureData=featuredata)



##---- creat exprement slot --------
##---- creat Model slot ------------
#model = data.frame(model.id=)
