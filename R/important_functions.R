convertListToDataFram <- function(inLst)
{
  ##-- this will convert list to data.fram
  mxInd = max(sapply(inLst, length))
  rlst  = as.data.frame(do.call(rbind,lapply(inLst, `length<-`, mxInd )), stringsAsFactors = FALSE)
  #rlst  = data.frame(as.matrix(rlst), stringsAsFactors = FALSE)
  #colnames(rlst) <- names(lst[[which.max(indx)]])
  return(rlst)
}


##---- gives index of element in vector
##---- also works for NA
getIndex <- function(inVec, indxOf)
{
  if(is.na(indxOf)){return(which(is.na(inVec)))}
  return(which(inVec==indxOf))
}