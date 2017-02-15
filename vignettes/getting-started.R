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

