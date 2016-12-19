
creatExperimentDesign <- function(expSlot)
{

  expModDrg = lapply(expSlot, function(x){ c(x$experiment.id, x$model.id, x$drug$join.name, x$batch, x$exp.type)})
  expModDrg = as.data.frame(do.call(rbind, expModDrg), stringsAsFactors = FALSE)
  colnames(expModDrg) = c("experiment.id", "model.id", "drug.join.name", "batch", "exp.type")

  dxBat = as.data.frame(table(expModDrg[expModDrg$exp.type=="treatment", c("drug.join.name", "batch") ]))
  dxBat = dxBat[dxBat$Freq>0, ]

  rtx=list()

  for(I in 1:dim(dxBat)[1])
  {
    #print(I)
    Lx = list(batch = as.character(dxBat[I, "batch"]))
    modBatchx = expModDrg[expModDrg$batch == Lx$batch, ]
    Lx$control = as.character(modBatchx[modBatchx$exp.type=="control",  "model.id"] )

    Lx$drug.join.name = as.character(dxBat[I, "drug.join.name"])
    Lx$treatment = as.character(modBatchx[modBatchx$exp.type=="treatment" &
                                            modBatchx$drug.join.name==Lx$drug.join.name, "model.id"])


    #namx = sprintf("%s.%s", Lx$drug.join.name, Lx$batch) ##-- no need for name here
    namx = length(rtx)+1
    rtx[[namx]]= Lx
  }
  rtx = .checkExperimentDesign(rtx)
  return(rtx)
}


.checkExperimentDesign<- function(expDesign)
{
  modNoControl = c()
  modNoTreatme = c()
  for(I in expDesign)
  {
    if(length(I$control)==0 & length(I$treatment)==0)
    {stop("Error Treatmetn and Control are missing in expDesign!")}

    if(length(I$control)==0)  { modNoControl = c(modNoControl, I$treatment)}
    if(length(I$treatment)==0){ modNoTreatme = c(modNoTreatme, I$control)}
  }


  if(!is.null(modNoControl))
  {
    txt = sprintf("These models have no Controls\n%s", paste(unique(modNoControl), collapse = "\n"))
    cat(txt)
  }
  return(expDesign)
}








