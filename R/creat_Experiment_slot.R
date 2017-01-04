##to for datafram to Exp slot

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
    drug[["names"]] = as.list(drug.N)
  }

  rtx$drug = drug

  ##--- elements with unique function -----------------
  eleUnq = c("formula", "tumor.type", "batch", "exp.type")
  for (u in eleUnq)
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
  dataColName = c("time", "volume", "width","length", doseColsNames, "body.weight", "date")
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
  rtxData$body.weight= as.numeric(rtxData$body.weight)
  rtxData$date  = as.Date(rtxData$date)

  rtxData[ ,doseColsNames]= sapply(doseColsNames, function(x){as.numeric(rtxData[ ,x])} )

  rtxData = BBmisc::sortByCol(rtxData , dataColName, asc = rep(TRUE, length(dataColName)))

  rtx$data= rtxData
  return(rtx)
}


#' @export
experimentSlotfromDf <- function(experiment)
{
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
                   "formula", "body.weight" ) #, "tumor.type", "batch", "exp.type")
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
  if(any(is.na(u.modDrg.id$model.id)))
  { stop("model.id is NA") }

  if(any(is.na(u.modDrg.id$drug)))
  { stop("drug is NA") }

  for (i in 1:dim(u.modDrg.id)[1])
  {
    exp.mod.dg = subset(experiment,
                     experiment$model.id== u.modDrg.id[i, "model.id"] &
                     experiment$drug == u.modDrg.id[i, "drug"] )

    expSlot[[i]] = creatListFromDF(exp.mod.dg)
  }

  mod.ids = unlist(sapply(expSlot , "[[" , "model.id" ))
  if(any(table(mod.ids)!=1))
  {
    msg = sprintf("These model.id are repeated\n%s", paste(mod.ids[table(mod.ids)!=1], collapse = ', '))
    stop(msg)
  }

  expNames = make.unique(unlist(sapply(expSlot, function(x){ sprintf("%s.%s", x$model.id, x$drug$join.name)} )), sep="_")
  for(i in 1:length(expSlot))
  { expSlot[[i]][["experiment.id"]] = expNames[i] }
  names(expSlot) = expNames

  return(expSlot)
}

