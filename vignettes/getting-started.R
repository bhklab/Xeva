## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load, echo=TRUE-----------------------------------------------------
library(Xeva)
cdf = readRDS("~/CXP/Xeva/DATA-raw/celineData.Rds")
head(ModelInfo(cdf))

## ---- echo=TRUE----------------------------------------------------------
print(batchNames(cdf)) 

## ---- echo=TRUE----------------------------------------------------------
batchNames = batchNames(cdf)
expDesign  = expDesign(cdf, batchNames[4])
ang = calculateAngle(cdf, expDesign, plot=TRUE)
print(ang)

## ---- echo=TRUE, fig.width = 10, fig.height = 10-------------------------
data(pdxe)
df = getmRECIST(pdxe)
## add tumor.type information
dfMap = mapModelSlotIds(object=pdxe, id=df$model.id, id.name="model.id", map.to="tumor.type")
dfx = merge(df, dfMap, by.x = "model.id", by.y = "model.id")
lungDf = dfx[dfx$tumor.type=="NSCLC", ]#c("model.id", "biobase.id", "drug.join.name", "mRECIST", "tumor.type")]
#pdf(file="DATA-raw/mRECIST_plot_NSCLC.pdf", width=12, height=10)
plotmRECIST(lungDf, groupBy = "biobase.id", control.name = "untreated")
#dev.off()

