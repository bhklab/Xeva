##to for datafram to Exp slot
#library(BBmisc)
creatListFromDF <- function(exp.mod.dg)
{
  rtx = list()
  exp.mod.dg <- data.frame(lapply(exp.mod.dg, as.character), stringsAsFactors=FALSE)

  rtx$model.id = unique(exp.mod.dg$model.id)

  drgColName.No = colnames(exp.mod.dg)[grep("drug\\.",colnames(exp.mod.dg))]

  drug = list("join.name" = unique(exp.mod.dg$drug))
  if(length(drgColName.No)>1)
  {
    drug.N = sapply(drgColName.No, function(x){unique(exp.mod.dg[,x])  })
    drug.N = drug.N[!is.na(drug.N)]
    drug[["names"]] = drug.N
  }

  rtx$drug = drug

  ##--- elements with unique function
  #for (u in c("batch", "exp.type", "formula", "tumor.type"))
  for (u in c("formula", "tumor.type"))
  {
    if(is.element(u, colnames(exp.mod.dg)))
    {
      rtx[[u]] = unique(exp.mod.dg[,u])
      if(length(rtx[[u]])>1)
      {
        msg = sprintf("For model.id==%s and drug==%s, %s is not unique",
                      rtx$model.id, rtx$drug$join.name, u)
        stop(msg)
      }
    }
  }

  #if(!is.null(rtx$exp.type) & !is.element(rtx$exp.type, c("control", "treatment")))
  #{
  #  msg = sprintf("For model.id %s and drug %s, exp.type should be control OR treatment
  #                it is %s\n", rtx$model.id, rtx$drug$join.name, rtx$exp.type)
  #  stop(msg)
  #}

  doseColsNames = c("dose", gsub("drug", "dose", names(rtx$drug$names)))
  dataColName = c("time", "volume", "width","length", doseColsNames, "weight", "date")
  for (w in dataColName)
  {
      if(is.element(w, colnames(exp.mod.dg))==FALSE)
      {
        exp.mod.dg[,w] = NA
      }
  }

  ##---- should we add dose.1 + dose.2 .... to dose

  rtxData = data.frame(lapply(exp.mod.dg[,dataColName], as.character), stringsAsFactors=FALSE)
  ##------ change column type for each column -----------------------------------------------
  rtxData$time  = as.numeric(rtxData$time)
  rtxData$volume= as.numeric(rtxData$volume)
  rtxData$width = as.numeric(rtxData$width)
  rtxData$length= as.numeric(rtxData$length)
  rtxData$weight= as.numeric(rtxData$weight)
  rtxData$date  = as.Date(rtxData$date)

  rtxData[ ,doseColsNames]= sapply(doseColsNames, function(x){as.numeric(rtxData[ ,x])} )

  rtxData = BBmisc::sortByCol(rtxData , dataColName, asc = rep(TRUE, length(dataColName)))

  rtx$data= rtxData
  return(rtx)
}


#' @export
experimentSlotfromDf <- function(experiment)
{
  #experiment = readRDS("data/DC_pdxDf_test_object.Rda")

  response = "volume" ## Default response colomn
  drugColsName = colnames(experiment)[grep("drug",colnames(experiment))]

  requredCols = c("model.id", "time", response, drugColsName)
  colAbsent = setdiff(requredCols, colnames(experiment))
  if(length(colAbsent)>0)
  {
    msg = sprintf("These colums are required\n%s", paste(colAbsent, collapse = ', '))
    stop(msg)
  }


  if(length(drugColsName)==0)
  {
    msg = sprintf("Column with drug information requred\nDrug infromation column should be named drug, drug.1 ...\n")
    stop(msg)
  } else{
    msg = sprintf("Drug columns are\n%s\n", paste(drugColsName, collapse = ', '))
    cat(msg)
  }

  doseColsName = colnames(experiment)[grep("dose",colnames(experiment))]
  standardCols = c(requredCols, doseColsName, "width","length", "date", "time",
                   "formula", "body.weight", "tumor.type") #"batch", "exp.type"
  extraCol = setdiff(colnames(experiment), standardCols)
  if(length(extraCol)>0)
  {
    msg = sprintf("These colums are not part of standard information,
                  therefor will be stored but not processed\n%s\n", paste(extraCol, collapse = ', '))
    warning(msg)
  }

  ##---- reformat drug column -----------
  drgColName.No = colnames(experiment)[grep("drug\\.",colnames(experiment))]
  if(length(drgColName.No)>0)
  {
    msg = sprintf("drug column will be replaced by %s\n", paste(drgColName.No, collapse = " + "))
    cat(msg)
    pasteWithoutNA <- function(L, collapse = " + "){paste(L[!is.na(L)], collapse = collapse)}
    experiment[, "drug"] = apply(experiment[,drgColName.No], 1, pasteWithoutNA)
  }

  ##------- if drug names are already in drug1 + drug2 split them ----------

  expSlot = list()
  u.modDrg.id = unique(experiment[, c("model.id", "drug")])
  for (i in 1:dim(u.modDrg.id)[1])
  {
    exp.mod.dg = subset(experiment,
                     experiment$model.id== u.modDrg.id[i, "model.id"] &
                     experiment$drug == u.modDrg.id[i, "drug"] )

    expSlot[[i]] = creatListFromDF(exp.mod.dg)
  }

  mod.ids = sapply(expSlot , "[[" , "model.id" )
  if(any(table(mod.ids)!=1))
  {
    msg = sprintf("These model.id are repeated\n%s", paste(mod.ids[table(mod.ids)!=1], collapse = ', '))
    stop(msg)
  }

  expNames = make.unique(sapply(expSlot, function(x){ sprintf("%s.%s", x$model.id, x$drug$join.name)} ), sep="_")
  names(expSlot) = mod.ids
  for(i in 1:length(expSlot))
  { expSlot[[i]][["experiment.id"]] = expNames[i] }

  return(expSlot)
}





###=================================================================
if(1==2){
checkModel <- function(model, expSlot)
{
  reqColName = c("model.id", "batch", "exp.type", "biobase.id")#, "drug.join.name")
  if(all(reqColName %in% colnames(model))==FALSE)
  {
    msg = sprintf("The required colmns for model are\n%s", paste(reqColName, collapse = ', '))
    stop(msg)
  }

  ##---- add drug column ---------------
  model$drug.join.name = sapply(model$model.id, function(x){ expSlot[[x]]$drug$join.name})
  return(model)
}



creatExperimentDesign <- function(model, expSlot)
{
  drgNames  = unique(sapply(expSlot, "[[", c("drug", "join.name")))
  tretBatch = unique(model[,  "batch"])

  rtx=list()
  for(drgI in drgNames)
  {
    for(batI in tretBatch)
    {
      Lx = list(drug.join.name = drgI, batch = batI)
      Lx$treatment = unique(model[model$drug.join.name == drgI &
                                    model$batch == batI &
                                    model$exp.type == "treatment", "model.id"] )

      Lx$control = unique(model[model$batch == batI & model$exp.type == "control", "model.id"] )

      if(length(Lx$treatment) > 0)
        {
        namx = sprintf("%s.%s", drgI, batI)
        rtx[[namx]]= Lx
        }
    }
  }

  trG = sapply(rtx, "[[", "treatment", simplify = TRUE)
  cnG = sapply(rtx, "[[", "control", simplify = TRUE)
  modInExpDes = unique(unlist(list(trG,cnG), recursive=TRUE))
  modInExpSlot= sapply(expSlot, "[[", "model.id")
  modInExpSlotNotPre = modInExpSlot[!modInExpSlot %in% modInExpDes]

  if(length(modInExpSlotNotPre)>0)
  {
    msg = sprintf("These model ids are not present in experiment design\n%s", paste(modInExpSlotNotPre, collapse = ', '))
    stop(msg)
  }

  return(rtx)
}

#experiment = readRDS("data/DC_pdxDf_test_object.Rda")
#DC_pdxDf_exp_object = experimentSlotfromDf(readRDS("data/DC_pdxDf_test_object.Rda"))
#saveRDS(DC_pdxDf_exp_object, file= "DC_pdxDf_ExpSlot_object.Rda")

Geo_Exp = readRDS("data/Geo_Exp.Rda")
Gao_PharPx_obj = list()
Gao_PharPx_obj[["experiment"]] = experimentSlotfromDf(Geo_Exp$experiment)
Gao_PharPx_obj[["model"]] = checkModel(Geo_Exp$model, Gao_PharPx_obj[["experiment"]])

Gao_PharPx_obj[["expDesign"]] = creatExperimentDesign(Gao_PharPx_obj[["model"]], Gao_PharPx_obj[["experiment"]])


saveRDS(Gao_PharPx_obj, file= "data/Gao_PharPx_obj.Rda")
}

















