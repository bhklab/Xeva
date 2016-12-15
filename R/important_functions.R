#' @export
.convertListToDataFram <- function(inLst)
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


.castDataFram <- function(df, row.var, col.var, value)
{
  uqR = unique(as.character(df[,row.var]))
  uqC = unique(as.character(df[,col.var]))
  dfx = data.frame( matrix(NA, nrow = length(uqR), ncol = length(uqC) ))
  rownames(dfx) = uqR ; colnames(dfx) = uqC
  for(I in 1:dim(df)[1])
  {
    dfx[as.character(df[I,row.var]), as.character(df[I,col.var])] = df[I, value]
  }
  return(dfx)
}


.appendToList <- function(in.list, value)
{
  in.list[[length(in.list)+1]]=value
  return(in.list)
}
