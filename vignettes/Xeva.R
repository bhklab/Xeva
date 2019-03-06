## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle::latex()

## ----options, results='hide', message=FALSE, eval=TRUE, echo=FALSE---------
library(Xeva)

## ----get_lib, results='hide', message=FALSE, eval=FALSE--------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("Xeva", version = "3.8")

## ----githubInst, results='hide', message=FALSE, eval=FALSE-----------------
#  #install devtools if required
#  install.packages("devtools")
#  
#  #install Xeva as:
#  devtools::install_github("bhklab/Xeva")

## ----l, results='hide', message=FALSE, eval=TRUE---------------------------
library(Xeva)

## ----biobase, results='hide', message=FALSE, eval=TRUE, echo=FALSE---------
suppressMessages(library(Biobase))

## ----l2--------------------------------------------------------------------
data(brca)
print(brca)

## ----l3--------------------------------------------------------------------
brca.mod <- modelInfo(brca)
dim(brca.mod)
brca.mod[1:4, ]

## ----expre-----------------------------------------------------------------
model.data <- getExperiment(brca, model.id = "X.1004.BG98")
head(model.data)

## ----batch1----------------------------------------------------------------
batch.name <- batchInfo(brca)
batch.name[1:4]

## ----batch2----------------------------------------------------------------
batchInfo(brca, batch = "X-1004.binimetinib")

## ----plot1, fig.cap="Tumor growth curves for a batch of control and treated PDXs", out.width='4in', fig.wide=TRUE----
plotPDX(brca, batch = "X-4567.BKM120")

## ----pdxplot2, fig.cap="Tumor growth curves for a batch of control and treated PDXs. Here, the volume is normalized and plots are truncated at 40 days", out.width='4in', fig.wide=TRUE----
plotPDX(brca, batch = "X-4567.BKM120", vol.normal = TRUE, control.col = "#a6611a",
        treatment.col = "#018571", major.line.size = 1, max.time = 40)

## ----pdxplot3, fig.cap="Tumor growth curves for a batch of control and treated PDXs generated using patient ID and drug name", out.width='4in', fig.wide=TRUE----
plotPDX(brca, patient.id="X-3078", drug="paclitaxel",control.name = "untreated")

## ----l4--------------------------------------------------------------------
brca.mr <- summarizeResponse(brca, response.measure = "mRECIST")
brca.mr[1:5, 1:4]

## ----mR_BRCA, fig.cap="mRECIST plot for PDXE breast cancer data", fig.width=14.1, fig.height=7.8, fig.wide=TRUE----
plotmRECIST(brca.mr, control.name="untreated", row_fontsize=13, col_fontsize=12)

## ----waterFall1, fig.cap="Waterfall plot for binimetinib drug response in PDXs", fig.width=14.1, fig.height=7.8, fig.wide=TRUE----
waterfall(brca, drug="binimetinib", res.measure="best.average.response")

## ----waterFall2, fig.cap="Waterfall plot for BYL719 drug response in PDXs", fig.width=14.1, fig.height=7.8, fig.wide=TRUE----
mut <- summarizeMolecularProfiles(brca,drug = "BYL719", mDataType="mutation")
model.type <- Biobase::exprs(mut)["CDK13", ]
model.type[grepl("Mut", model.type)] <- "mutation"
model.type[model.type!="mutation"] <- "wild type"
model.color <- list("mutation"="#b2182b", "wild type"="#878787")
waterfall(brca, drug="BYL719", res.measure="best.average.response",
          model.id=names(model.type), model.type= model.type,
          type.color = model.color)


## ----x1--------------------------------------------------------------------
model=read.csv(system.file("extdata", "model.csv", package = "Xeva"))
head(model)

## ----x2--------------------------------------------------------------------
drug = read.csv(system.file("extdata", "drug.csv", package = "Xeva"))
head(drug)

## ----x3--------------------------------------------------------------------
experiment=read.csv(system.file("extdata","experiments.csv",package="Xeva"))
head(experiment)

## ----x4--------------------------------------------------------------------
expDesign=readRDS(system.file("extdata","batch_list.rds",package="Xeva"))
expDesign[[1]]

## ----x5--------------------------------------------------------------------
RNASeq=readRDS(system.file("extdata", "rnaseq.rds", package = "Xeva"))
print(RNASeq)

## ----x6--------------------------------------------------------------------
modToBiobaseMap=read.csv(system.file("extdata","modelToExpressionMap.csv",
                                     package = "Xeva"))
head(modToBiobaseMap)

## ----createXevaSet---------------------------------------------------------
xeva.set=createXevaSet(name="example xevaSet", model=model, drug=drug,
                       experiment=experiment, expDesign=expDesign,
                       molecularProfiles=list(RNASeq = RNASeq),
                       modToBiobaseMap = modToBiobaseMap)
print(xeva.set)

## ----sess------------------------------------------------------------------
sessionInfo()

