## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load, echo=TRUE-----------------------------------------------------
library(Xeva)
data(lpdx)
head(modelInfo(lpdx))

## ---- echo=TRUE----------------------------------------------------------
print(batchNames(lpdx)) 

## ---- echo=TRUE,fig.width = 12, fig.height = 12--------------------------
batchNames <- batchNames(lpdx)
expDesign  <- expDesign(lpdx, batchNames[1])
ang <- calculateAngle(lpdx, expDesign, treatment.only = TRUE, plot=TRUE)
print(ang)

#par(mfrow=c(5,3))
for(I in batchNames)
{
  expDesign  <- expDesign(lpdx, I)
  ang <- calculateAngle(lpdx, expDesign, treatment.only = TRUE, plot=TRUE)
#  print(ang)
}


## ---- echo=TRUE----------------------------------------------------------
lpdx_slop <- summarizeResponse(lpdx, response.measure = "slop", 
                               group.by="patient.id", summary.stat = "mean")



## ---- echo=TRUE----------------------------------------------------------
lpdx_angle <- summarizeResponse(lpdx, response.measure = "angle")


## ---- echo=TRUE----------------------------------------------------------
ldxe_mut <- getMolecularProfiles(lpdx, data.type="mutation")
print(ldxe_mut)

## ---- echo=TRUE----------------------------------------------------------
# get sample names
library(Biobase)
sn <- Biobase::sampleNames(ldxe_mut)
smap <- mapModelSlotIds(lpdx, id=sn, id.name = "biobase.id", map.to = "model.id")
head(smap)

## ---- echo=TRUE----------------------------------------------------------
df = getExperiment(lpdx,"PHLC119_P5.506.B1.3")
print(df[df$time>85 & df$time<109, c("time", "width", "length", "volume", "comment", "dose")])

## ---- echo=TRUE, fig.width = 12, fig.height = 10-------------------------
data(pdxe)
df <- getmRECIST(pdxe)
## add tumor.type information
dfMap <- mapModelSlotIds(object=pdxe, id=df$model.id, id.name="model.id",
                        map.to="tumor.type", unique = FALSE)
if(all(df$model.id==dfMap$model.id)) {df$tumor.type = dfMap$tumor.type}
lungDf = df[df$tumor.type=="NSCLC", ]
#pdf(file="DATA-raw/mRECIST_plot_NSCLC.pdf", width=12, height=10)
plotmRECIST(lungDf, groupBy = "biobase.id", control.name = "untreated")

#pdf(file="DATA-raw/mRECIST_plot_BRCA.pdf", width=12, height=10)
brDF = df[df$tumor.type=="BRCA", ]
plotmRECIST(brDF, groupBy = "biobase.id", control.name = "untreated")
#dev.off()

## ---- echo=TRUE, fig.width = 12, fig.height = 10-------------------------
data(pdxe)

pm = modelInfo(pdxe)
lungPID = unique(pm[pm$tumor.type=="NSCLC", "patient.id"])

pdxe_slop <- summarizeResponse(pdxe, response.measure = "slop", 
                               group.by="patient.id", summary.stat = "mean")

lung_pdxe_slope <- pdxe_slop[, lungPID]

##-------
pdxe_mR <- summarizeResponse(pdxe, response.measure = "mRECIST_recomputed", 
                               group.by="patient.id")

lung_pdxe_mR = pdxe_mR[, lungPID]

slope=c(); mR=c()
for(dn in rownames(lung_pdxe_slope))
{
  for(pi in colnames(lung_pdxe_slope))
  {
    v = c(lung_pdxe_slope[dn,pi], lung_pdxe_mR[dn,pi])
    if(!is.na(v[1]) & !is.na(v[1]))
    { slope = c(slope,v[1]); mR=c(mR,v[2]) }
  }
}

df = data.frame(mR= mR, slope= as.numeric(slope), stringsAsFactors = FALSE)
df$mR= factor(df$mR, c("CR", "PR", "SD", "PD"))

colPalette = c("#377eb8", "#4daf4a", "#fec44f", "#e41a1c")
#pdf(file="DATA-raw/boxplot_lungCancer.pdf", width=12, height=10)
boxplot(slope~mR, data=df, col=colPalette,
  main="mRECIST vs slope", xlab="mRECIST", ylab="slope")
#dev.off()

