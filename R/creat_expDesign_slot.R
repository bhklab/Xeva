


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

  ##------- setting name -----------------
  bnam = sapply(expDesign, "[[", "batch.name")
  bnamDup = bnam[duplicated(bnam)]
  if(length(bnamDup)>0)
  {
    txt = sprintf("These batch names are duplicated\n%s", paste(bnamDup, collapse = "\n"))
    stop(txt)
  }
  names(expDesign) = bnam

  return(expDesign)
}








