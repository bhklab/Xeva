## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=TRUE----------------------------------------------------------
library(Xeva)
data(pdxe)

## ---- echo=TRUE----------------------------------------------------------
modelIds = getModelIds(pdxe, drug="paclitaxel", drug.match.exact=TRUE, tumor.type="BRCA")
print(modelIds)

## ---- echo=TRUE----------------------------------------------------------
df = getExperiment(pdxe, model.id="X.1655.LE11.biib")
head(df)

## ----cars----------------------------------------------------------------
summary(cars)

## ----pressureX, , include=TRUE, echo=FALSE-------------------------------
plot(pressure)

