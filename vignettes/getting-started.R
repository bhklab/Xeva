## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load, echo=TRUE-----------------------------------------------------
library(Xeva)
data(lpdx)
head(modelInfo(lpdx))

## ---- echo=TRUE----------------------------------------------------------
print(batchNames(lpdx)) 

## ---- echo=TRUE----------------------------------------------------------
batchNames <- batchNames(lpdx)
expDesign  <- expDesign(lpdx, batchNames[4])
ang <- calculateAngle(lpdx, expDesign, treatment.only = TRUE, plot=TRUE)
print(ang)

## ---- echo=TRUE----------------------------------------------------------
lpdx_slop <- summarizeResponse(lpdx, response.measure = "slop", 
                               group.by="patient.id", summary.stat = "mean")



## ---- echo=TRUE----------------------------------------------------------
lpdx_angle <- summarizeResponse(lpdx, response.measure = "angle")


## ---- echo=TRUE, fig.width = 12, fig.height = 10-------------------------
data(pdxe)
df <- getmRECIST(pdxe)
## add tumor.type information
dfMap <- mapModelSlotIds(object=pdxe, id=df$model.id, id.name="model.id",
                        map.to="tumor.type", unique = FALSE)
#dfx = merge(df, dfMap, by.x = "model.id", by.y = "model.id")
if(all(df$model.id==dfMap$model.id)) {df$tumor.type = dfMap$tumor.type}
lungDf = df[df$tumor.type=="NSCLC", ]
#pdf(file="DATA-raw/mRECIST_plot_NSCLC.pdf", width=12, height=10)
plotmRECIST(lungDf, groupBy = "biobase.id", control.name = "untreated")

