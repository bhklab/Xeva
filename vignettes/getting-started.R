## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load, echo=TRUE-----------------------------------------------------
library(Xeva)
cdf = readRDS("~/CXP/Xeva/DATA-raw/celineData.Rds")
head(ModelInfo(cdf))

## ---- echo=TRUE----------------------------------------------------------
print(cdf@expDesign[[1]]) 

## ---- echo=TRUE----------------------------------------------------------
batch = cdf@expDesign[[1]]
ang = calculateAngle(cdf, batch, plot=TRUE)
print(ang)

