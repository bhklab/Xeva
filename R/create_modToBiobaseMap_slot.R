.checkmodToBiobaseMapSlot <- function(modToBiobaseMap, molecularProfiles)
{
  if(nrow(modToBiobaseMap) > 0 & length(molecularProfiles)>0)
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
    #return(modToBiobaseMap)
  } else
  {
    if(nrow(modToBiobaseMap) > 0)
    { warning("modToBiobaseMap not present")}

    if(length(molecularProfiles)>0)
    { warning("molecularProfiles not present")}
  }

  return(modToBiobaseMap)
}
