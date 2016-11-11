##to for datafram to Exp slot
library(BBmisc)
creatListFromDF <- function(exp.mod.dg)
{
  rtx = list()
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
  for (u in c("batch", "type", "formula"))
  {
    if(is.element(u, colnames(exp.mod.dg)))
    {
      rtx[[u]] = unique(exp.mod.dg[,u])
      if(length(rtx[[u]])>1)
      {
        msg = sprintf("For model.id==%s and drug==%s, %s is not unique", rtx$model.id, rtx$drug, u)
        stop(msg)
      }
    }
  }

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

  rtxData = data.frame(exp.mod.dg[,dataColName], stringsAsFactors = FALSE)
  rtxData = sortByCol(rtxData , dataColName, asc = rep(TRUE, length(dataColName)))
  #rtxData$width = as.numeric(rtxData$width)
  rtx$data= rtxData
  return(rtx)
}


#' @export
experimentSlotfromDf <- function(experiment)
{
  experiment = readRDS("data/DC_pdxDf_test_object.Rda")

  response = "volume" ## Default response colomn
  requredCols = c("model.id", "time", response, drugColsName)
  colAbsent = setdiff(requredCols, colnames(experiment))
  if(length(colAbsent)>0)
  {
    msg = sprintf("These colums are required\n%s", paste(colAbsent, collapse = ', '))
    stop(msg)
  }

  drugColsName = colnames(experiment)[grep("drug",colnames(experiment))]
  if(length(drugColsName)==0)
  {
    msg = sprintf("Column with drug information requred\nDrug infromation column should be named drug, drug.1 ...")
    stop(msg)
  } else{
    msg = sprintf("Drug columns are\n%s", paste(drugColsName, collapse = ', '))
    cat(msg)
  }

  doseColsName = colnames(experiment)[grep("dose",colnames(experiment))]
  standardCols = c(requredCols, doseColsName, "width","length", "date", "time", "batch", "type", "formula", "body.weight")
  extraCol = setdiff(colnames(experiment), standardCols)
  if(length(extraCol)>0)
  {
    msg = sprintf("These colums are not part of standard information,
                  therefor will be stored but not processed\n%s", paste(extraCol, collapse = ', '))
    warning(msg)
  }

  ##---- reformat drug column -----------
  drgColName.No = colnames(experiment)[grep("drug\\.",colnames(experiment))]
  if(length(drgColName.No)>0)
  {
    msg = sprintf("drug column will be replaced by %s", paste(drgColName.No, collapse = " + "))
    cat(msg)
    pasteWithoutNA <- function(L,collapse = " + "){paste(L[!is.na(L)], collapse = collapse)}
    experiment[, "drug"] = apply(dinfo, 1, pasteWithoutNA)
  }

  expSlot = list()
  u.modDrg.id = unique(experiment[, c("model.id", "drug")])
  for (i in 1:dim(u.modDrg.id)[1])
  {
    exp.mod.dg = subset(experiment,
                     experiment$model.id== u.modDrg.id[i, "model.id"] &
                     experiment$drug == u.modDrg.id[i, "drug"] )

    expSlot[[i]] = creatListFromDF(exp.mod.dg)
  }

  return(expSlot)
}

