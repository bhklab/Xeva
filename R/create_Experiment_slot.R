.createListFromDF <- function(exp.mod.dg, extraCol=NULL)
{
  rtx <- list()
  exp.mod.dg <- data.frame(lapply(exp.mod.dg, as.character), stringsAsFactors=FALSE)

  rtx$model.id = unique(exp.mod.dg$model.id)

  drgColName.No = colnames(exp.mod.dg)[grep("drug",colnames(exp.mod.dg))]

  drug = list("join.name" = unique(exp.mod.dg$drug))
  if(length(drgColName.No)>1)
  {
    drug.N <- vapply(drgColName.No, function(x){unique(exp.mod.dg[,x])}, FUN.VALUE = character(1))
    drug.N <- drug.N[!is.na(drug.N)]
    drug[["names"]] = as.list(drug.N)
  }

  rtx$drug = drug

  ##------------ set extra col ------------------------------------
  if(!is.null(extraCol))
  {
    for(ec in c(extraCol))
    {
      vx = exp.mod.dg[, ec]
      if(length(unique(vx))==1)
      { vx <- unique(vx) }
      rtx[[ec]] <- vx
    }
  }

  doseColsNames <- c("dose", gsub("drug", "dose", names(rtx$drug$names)))
  dataColName <- c("time", "volume", "width","length",
                  doseColsNames, "body.weight", "date", "comment")
  dataColName <- unique(c(dataColName, colnames(exp.mod.dg)))

  naCols <- setdiff(dataColName, colnames(exp.mod.dg))
  exp.mod.dg[, naCols] <- NA

  ##---- add dose.1 + dose.2 .... to dose
  rtxData <- data.frame(lapply(exp.mod.dg[,dataColName], as.character),
                        stringsAsFactors=FALSE)
  ##------ change column type for each column ---------------------------
  rtxData$time  <- as.numeric(rtxData$time)
  rtxData$volume<- as.numeric(rtxData$volume)
  rtxData$width <- as.numeric(rtxData$width)
  rtxData$length<- as.numeric(rtxData$length)
  rtxData$body.weight<- as.numeric(rtxData$body.weight)
  rtxData$date  <- as.Date(rtxData$date)

  for(doseCi in doseColsNames)
  {
    rtxData[ ,doseCi] <- as.numeric(rtxData[ ,doseCi])
  }

  rtxData <- BBmisc::sortByCol(rtxData , "time", asc = TRUE)
  rtx$data<- rtxData
  return(rtx)
}

## list of PDXMI variables
## Source
## Meehan, Terrence F., et al. "PDX-MI: minimal information for patient-derived
## tumor xenograft models." Cancer research 77.21 (2017): e62-e66.
## http://cancerres.aacrjournals.org/lookup/doi/10.1158/0008-5472.CAN-17-0582
modelClassS4Vars <- function()
{
  return(
    c("model.id", "drug", "data", "treatment.type", "treatment.target",
      "patient.id", "patient.sex", "patient.age", "patient.diagnosis",
      "patient.consent.to.share.data", "patient.ethnicity",
      "patient.current.treatment.drug",
      "patient.current.treatment.protocol", "patient.prior.treatment.protocol",
      "patient.response.to.prior.treatment", "patient.virology.status",

      "tumor.id", "tumor.tissue.of.origin",
      "tumor.primary.metastasis.recurrence",
      "tumor.specimen.tissue", "tumor.tissue.histology", "tumor.tumor.grade",
      "tumor.disease.stage", "tumor.specific.markers",
      "tumor.fom.untreated.patient",
      "tumor.original.sample.type", "tumor.from.existing.pdx.model",

      "model.submitter.pdx.id", "model.mouse.strain.source",
      "model.strain.immune.system.humanized",
      "model.type.of.humanization", "model.tumor.preparation",
      "model.injection.type.and.site",
      "model.mouse.treatment.for.engraftment", "model.engraftment.rate",
      "model.engraftment.time",

      "model.tumor.characterization.technology",
      "model.tumor.confirmed.not.to.be.of.mouse.origin",
      "model.response.to.standard.of.care",
      "model.animal.health.status", "model.passage.qa.performed",

      "model.treatment.passage", "model.treatment.protocol",
      "model.treatment.response", "model.tumor.omics",
      "model.development.of.metastases.in.strain",
      "model.doubling.time.of.tumor",

      "pdx.model.availability", "governance.restriction.for.distribution",
      "id.publication.data")
  )
}

makePDXModClassS4 <- function(exp.mod.dg, extraCol)
{
  pdxS3 <- .createListFromDF(exp.mod.dg, extraCol=extraCol)
  extraInfo <- list()
  if(!is.null(extraCol))
  { extraInfo <- pdxS3[extraCol] }

  pdxS4 <- PDXmodClass(model.id = pdxS3$model.id, drug = pdxS3$drug,
                       data=pdxS3$data, meta=extraInfo)

  pS4SlN<- modelClassS4Vars()


  for(s in pS4SlN)
  {
    if(!is.null(pdxS3[[s]]))
    { slot(pdxS4, s) <- pdxS3[[s]] }
  }

  return(pdxS4)
}

###----- define standard column names -----------
.getColumnsDF <- function()
{
  standCol <- c("model.id", "drug", "time", "volume", "width","length",
               "date", "body.weight","formula")

  requredCols <- c("model.id", "time", "volume", "drug")
  rtz <- list(standCol=standCol,requredCols=requredCols)
  return(rtz)
}


experimentSlotfromDf <- function(experiment)
{
  clnm <- .getColumnsDF()
  drugColsName <- colnames(experiment)[grep("drug",colnames(experiment))]

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

  doseColsName <- colnames(experiment)[grep("dose",colnames(experiment))]
  if(length(doseColsName)==0)
  {
    msg = sprintf("No dose column found\n")
    #warning(msg)
  }

  standardCols <- unique(unlist(c(requredCols, doseColsName, "width","length",
                                  "date", "time", "body.weight", "comment",
                                  modelClassS4Vars())))

  extraCol <- setdiff(colnames(experiment), standardCols)
  if(length(extraCol)>0)
  {
    msg <- sprintf("These colums are not part of standard information, therefor will be stored but not processed\n%s\n", paste(extraCol, collapse = ', '))
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
  u.modDrg.id <- unique(experiment[, c("model.id", "drug")])
  if(any(is.na(u.modDrg.id$model.id)))
  { stop("model.id is NA") }

  mdup = u.modDrg.id$model.id[duplicated(u.modDrg.id$model.id)]
  if(length(mdup)>0)
  {
    msg = sprintf("Duplicated model.id\n%s\nuse different model.id for different drugs\n", paste(mdup, collapse = "\n"))
    stop(msg)
  }


  expSlot = list()
  for (i in seq_len(dim(u.modDrg.id)[1]))
  {
    exp.mod.dg <- subset(experiment,
                     experiment$model.id== u.modDrg.id[i, "model.id"] &
                     experiment$drug == u.modDrg.id[i, "drug"] )

    expSlot[[i]] <- makePDXModClassS4(exp.mod.dg, extraCol=extraCol)
  }

  mod.ids <- vapply(expSlot, function(mod){slot(mod, "model.id")}, FUN.VALUE = character(1))

  if(length(mod.ids) != length(unique(mod.ids)))
  {
    msg <- sprintf("These model.id are repeated\n%s",
                   paste(mod.ids[table(mod.ids)!=1], collapse = ', '))
    stop(msg)
  }
  names(expSlot) <- mod.ids

  return(expSlot)
}
