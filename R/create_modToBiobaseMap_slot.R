.createModToBiobaseMap <- function(modToBiobaseMap, molecularProfiles, model)
{
  w <- names(molecularProfiles)[!(names(molecularProfiles) %in% colnames(modToBiobaseMap))]
  if(length(w)>0)
  {
    msg <- sprintf("Molecular data type %s not present in modToBiobaseMap", paste(w, collapse = "\n"))
    stop(msg)
  }
  return(modToBiobaseMap)
}

