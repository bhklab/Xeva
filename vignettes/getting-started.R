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

