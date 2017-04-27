<<<<<<< HEAD
.checkmodToBiobaseMapSlot <- function(modToBiobaseMap, molecularProfiles)
{
  rqdCol <- c("model.id", "biobase.id", "mDataType")
  for(cx in rqdCol)
  {
    if(is.element(cx, colnames(modToBiobaseMap))==FALSE)
    {
      msg <- sprintf("column %s not present is modToBiobaseMap\nmodToBiobaseMap must have the columns\n%s\n",
                     cx, paste(rqdCol, collapse = "\n"))
      stop(msg)
    }
  }

  mbDataTypes <- unique(as.character(modToBiobaseMap$mDataType))
  w <- names(molecularProfiles)[!(names(molecularProfiles) %in% mbDataTypes)]
  if(length(w)>0)
  {
    msg <- sprintf("Id mapping for molecular data type %s not present in modToBiobaseMap", paste(w, collapse = "\n"))
    warning(msg)
  }

  return(modToBiobaseMap)
}
=======
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

>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
