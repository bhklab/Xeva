if(1==2){ print("this is running....")

geoExp = readRDS("data/Geo_Exp.Rda")
experiment = geoExp$experiment
model = geoExp$model

expSlot = experimentSlotfromDf(experiment)

model = checkModel(model, expSlot)

expDesign = creatExperimentDesign(model, expSlot)





















}