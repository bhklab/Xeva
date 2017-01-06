## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=TRUE----------------------------------------------------------
library(Xeva)
data(pdxe)

## ---- echo=TRUE----------------------------------------------------------
#modelIds = getModelIds(pdxe, drug="paclitaxel", drug.match.exact=TRUE, tumor.type="BRCA")
#print(modelIds)
df = getExperiment(pdxe, experiment.id="X.1004.pael.paclitaxel")
head(df)

## ---- echo=TRUE----------------------------------------------------------
df = getExperiment(pdxe, experiment.id="X.1004.pael.paclitaxel")
head(df)

## ---- echo=TRUE----------------------------------------------------------
library(Xeva)
data("celineData")
head(ModelInfo(lpdx))

## ----cars----------------------------------------------------------------
summary(cars)

## ----pressureX, , include=TRUE, echo=FALSE-------------------------------
plot(pressure)

