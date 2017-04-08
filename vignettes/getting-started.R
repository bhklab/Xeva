## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load, echo=TRUE-----------------------------------------------------
library(Xeva)
data(lpdx)

## ---- echo=TRUE----------------------------------------------------------
lpdx.mod = modelInfo(lpdx)
head(lpdx.mod$model.id)

## ---- echo=TRUE----------------------------------------------------------
modId = lpdx.mod$model.id[82]
df = getExperiment(lpdx, model.id = modId)
head(df)

## ---- echo=TRUE----------------------------------------------------------
df = getExperiment(lpdx, modId, treatment.only = TRUE)
head(df)

## ---- echo=TRUE----------------------------------------------------------
print(batchNames(lpdx))
df = getExperiment(lpdx, batch.name = batchNames(lpdx)[1], treatment.only = TRUE)
head(df)

## ---- echo=TRUE,fig.width = 12, fig.height = 12--------------------------
batchNames <- batchNames(lpdx)
expDesign  <- expDesign(lpdx, batchNames[1])
ang <- calculateAngle(lpdx, expDesign, treatment.only = TRUE, plot=TRUE)
print(ang)

for(I in batchNames)
{
  expDesign  <- expDesign(lpdx, I)
  ang <- calculateAngle(lpdx, expDesign, treatment.only = TRUE, plot=TRUE)
  print(ang)
}


## ---- echo=TRUE----------------------------------------------------------
#lpdx_slop <- summarizeResponse(lpdx, response.measure = "slop", 
#                               group.by="patient.id", summary.stat = "mean")

## ---- echo=TRUE----------------------------------------------------------
#lpdx_angle <- summarizeResponse(lpdx, response.measure = "angle")


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
df = getExperiment(lpdx, "PHLC119_P5.506.B1")
#print(df[df$time>85 & df$time<109, c("time", "width", "length", "volume", "comment", "dose")])

## ---- echo=TRUE, fig.width = 12, fig.height = 10-------------------------
data(pdxe)
#select lung cancer PDXE data
pdxe.lung <- summarizeResponse(pdxe, response.measure = "mRECIST",
                                group.by="patient.id", tumor.type="NSCLC")
## plot matrix
plotmRECIST(pdxe.lung, control.name = "untreated")

## ---- echo=TRUE, fig.width = 12, fig.height = 10-------------------------
data(pdxe)
#select lung cancer PDXE data
pdxe.brca <- summarizeResponse(pdxe, response.measure = "mRECIST",
                                group.by="patient.id", tumor.type="BRCA")

## plot matrix
plotmRECIST(pdxe.brca, control.name = "untreated", control.col = "#238b45")

## ---- echo=TRUE, fig.width = 12, fig.height = 10-------------------------
data(pdxe)

lung_pdxe_slope <- summarizeResponse(pdxe, response.measure = "slope", group.by="patient.id", 
                                     summary.stat = "mean", tumor.type = "NSCLC")

lung_pdxe_mR <- summarizeResponse(pdxe, response.measure = "mRECIST",
                             group.by="patient.id", tumor.type="NSCLC")

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
boxplot(slope~mR, data=df, col=colPalette, main="mRECIST vs slope", 
        xlab="mRECIST", ylab="slope")


