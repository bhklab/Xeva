##to for datafram to Exp slot

creatListFromDF <- function(exp.mod.dg, extraCol=NULL)
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

  ##------------ set extra col -----------------------------------------------------
  if(!is.null(extraCol))
  {
    for(ec in c(extraCol))
    {
      vx = exp.mod.dg[, ec]
      if(length(unique(vx))==1)
      { vx = unique(vx) }
      rtx[[ec]] = vx
    }
  }

  ##-------------------------------------------------------------------------
  doseColsNames = c("dose", gsub("drug", "dose", names(rtx$drug$names)))
  dataColName = c("time", "volume", "width","length",
                  doseColsNames, "body.weight", "date", "comment")
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

  ##------ remove all cols which are NA -----------
  rtxData = rtxData[, colSums(is.na(rtxData)) != nrow(rtxData)]

  rtx$data= rtxData

  return(rtx)
}


###----- define standard column names -----------
.getColumnsDF <- function()
{
  standCol = c("model.id", "drug", "time", "volume", "width","length",
               "date", "body.weight","formula")

  requredCols = c("model.id", "time", "volume", "drug")
  rtz = list(standCol=standCol,requredCols=requredCols)
  return(rtz)
}


## @export
experimentSlotfromDf <- function(experiment)
{
  clnm = .getColumnsDF() ## change this Define all columns in function getColumnsDF and use it

  drugColsName = colnames(experiment)[grep("drug",colnames(experiment))]

  requredCols = c("model.id", "time", "volume", drugColsName)
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
  if(length(doseColsName)==0)
  {
    msg = sprintf("No dose column found\n")
    warning(msg)
  }

  standardCols = c(requredCols, doseColsName, "width","length", "date", "time",
                   "body.weight", "comment" )
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


  u.modDrg.id = unique(experiment[, c("model.id", "drug")])
  if(any(is.na(u.modDrg.id$model.id)))
  { stop("model.id is NA") }

  mdup = u.modDrg.id$model.id[duplicated(u.modDrg.id$model.id)]
  if(length(mdup)>0)
  {
    msg = sprintf("Duplicated model.id\n%s\nuse different model.id for different drugs\n", paste(mdup, collapse = "\n"))
    stop(msg)
  }


  expSlot = list()
  for (i in 1:dim(u.modDrg.id)[1])
  {
    exp.mod.dg = subset(experiment,
                     experiment$model.id== u.modDrg.id[i, "model.id"] &
                     experiment$drug == u.modDrg.id[i, "drug"] )

    expSlot[[i]] = creatListFromDF(exp.mod.dg, extraCol=extraCol)
  }

  mod.ids = unlist(sapply(expSlot , "[[" , "model.id" ))
  if(any(table(mod.ids)!=1))
  {
    msg = sprintf("These model.id are repeated\n%s", paste(mod.ids[table(mod.ids)!=1], collapse = ', '))
    stop(msg)
  }

  #expNames = make.unique(unlist(sapply(expSlot, function(x){ sprintf("%s.%s", x$model.id, x$drug$join.name)} )), sep="_")
  #for(i in 1:length(expSlot))
  #{ expSlot[[i]][["experiment.id"]] = expNames[i] }

  expNames = unlist(sapply(expSlot, "[[", "model.id"))
  names(expSlot) = expNames

  return(expSlot)
}

