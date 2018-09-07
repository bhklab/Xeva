## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load, echo=TRUE-----------------------------------------------------
library(Xeva)
data(brca)
brca

## ---- echo=TRUE, tidy=TRUE-----------------------------------------------
brca.mod <- modelInfo(brca)
head(brca.mod)

## ---- echo=TRUE, tidy=TRUE-----------------------------------------------
brca.mr <- summarizeResponse(brca, response.measure = "mRECIST")
head(brca.mr)
plotmRECIST(brca.mr, control.name = "untreated")

## ---- echo=TRUE, tidy=TRUE-----------------------------------------------
plotBatch(brca, batch = "X-5249.BKM120")

plotBatch(brca, patient.id="X-5249", drug="BKM120", control.name="untreated")

