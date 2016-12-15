
summaryData <- function(pdxe)
{
  ##----------------------------
  totalExp = length(pdxe@experiment)
  controlName = sapply(pdxe@experiment, function(x){if(x$exp.type=="control"){return(x$experiment.id)}else{return(NA)}})
  controlName = controlName[!is.na(controlName)]

  treatmentName = sapply(pdxe@experiment, function(x){if(x$exp.type=="treatment"){return(x$experiment.id)}else{return(NA)}})
  treatmentName = treatmentName[!is.na(treatmentName)]


  ##--- count total drugs ------
  dn = lapply(pdxe@experiment, "[[", c("drug", "names"))

}
