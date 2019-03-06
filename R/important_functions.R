## given a list of data frames it will give a new data frame
## x <- data.frame(a=1:3, b= LETTERS[1:3])
## y <- data.frame(a=(1:3)*10, b= letters[1:3])
## inLst <- list(x=x, y=y)
## .rbindListOfDataframs(inLst)
.rbindListOfDataframs<- function(inList)
{
  allColNames = lapply(inList, colnames)
  allColNames = unique(unlist(allColNames))
  rtx = data.frame(); ncolX=0;
  for(dfi in inList)
  {
    if( length(colnames(dfi)) > ncolX)
    { ncolX= length(colnames(dfi))
    maxCols =colnames(dfi)}

    for(cx in allColNames)
    {
      if(is.element(cx, colnames(dfi))==FALSE)
      {dfi[, cx]=NA}
    }
    rtx = rbind(rtx, dfi[, allColNames])
  }

  rtx1 = rtx[, maxCols]
  rtx2 = rtx[, setdiff(colnames(rtx), maxCols)]
  RTz = cbind(rtx1, rtx2)
  return(RTz)
}


##---- gives index of element in vector
##---- also works for NA
getIndex <- function(inVec, indxOf)
{
  if(is.na(indxOf)){return(which(is.na(inVec)))}
  return(which(inVec==indxOf))
}


.castDataFram <- function(df, row.var, col.var, value, collapse = ";")
{
  uqR <- unique(as.character(df[,row.var]))
  uqC <- unique(as.character(df[,col.var]))
  dfx <- matrix(NA, nrow = length(uqR), ncol = length(uqC) )
  rownames(dfx) <- uqR ; colnames(dfx) <- uqC

  ml <- list()
  for(r in uqR)
  {
    cl <- list()
    for(c in uqC)
    { cl[[c]] <- c() }
    ml[[r]] <- cl
  }

  for(I in seq_len(nrow(df)))
  {
    v <- ml[[df[I, row.var]]][[df[I, col.var]]]
    ml[[df[I, row.var]]][[df[I, col.var]]] <- c(v, df[I, value]) #vx
  }

  for(r in rownames(dfx))
  {
    for(c in colnames(dfx))
    {
      v <- ml[[r]][[c]]
      if(length(v)==0){v <- NA}
      vz <- NULL
      if(collapse =="mean")
      {
        vz <- mean(as.numeric(v[!is.na(v)]))
      } else if (collapse =="median")
      {
        vz <- median(as.numeric(v[!is.na(v)]))
      } else
      {
        if(length(v)>1)
        {
          vz <- pasteWithoutNA(v, collapse = collapse)
        } else
        { vz <- v }
        if(!is.na(vz) & vz==""){vz <- NA}
      }
      if(is.nan(vz)){vz <- NA}
      dfx[r,c] <- vz
    }
  }
  dfx <- data.frame(dfx, stringsAsFactors = FALSE, check.names = FALSE)
  return(dfx)
}


.appendToList <- function(in.list, value)
{
  in.list[[length(in.list)+1]]=value
  return(in.list)
}

##---------------------------------------------------------------------------
##---------------remove NA col-----------------------------------------------
.removeNAcol <- function(df)
{ return(df[, !apply(is.na(df), 2, all)]) }

##------------------------------------------------------------------------
##---------------reorder column ------------------------------------------

.reorderCol <- function(df, columnName, newIndx)
{
  OtherCN = colnames(df)[colnames(df)!=columnName]
  newCN = append(OtherCN, columnName, after= (newIndx-1) )
  return(df[,newCN])
}

##-------------------------------------------------------------------------
##-------------------------------------------------------------------------
## paste vector elements together while removing NAs
## \code{pasteWithoutNA} paste vector elements together while removing NAs
## @param L A vector with values and NAs
## @param collapse Collapse string, default "+"
## @return  Returns a string with vector values pasted together
## @examples
## L = c("A", NA, "B", NA, NA, "C")
## pasteWithoutNA(L, collapse = " + ")
pasteWithoutNA <- function(L, collapse = " + "){paste(L[!is.na(L)],
                                                      collapse = collapse)}



##------------------------------------------------------------------------
##--- this will creat empty theme for ggplot -----------------------------
.ggplotEmptyTheme <- function(plt){
  plt +  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        panel.background = ggplot2::element_blank(),
                        legend.key=element_blank(), ##removes legend background
                        axis.line = ggplot2::element_line(colour = "black"))
}

##-------------------
# Function to print data.frame in message
#
# \code{printAndCapture} prints data.frame in stop or warning functions
# @examples
# df <- data.frame(a=1:5, b=11:15)
# msg <- sprintf("data frame is:\n%s", printAndCapture(df))
# warning(msg)
printAndCapture <- function(x)
{
  paste(capture.output(print(x)), collapse = "\n")
}

###----------------------------
##Normalize a vector between 0 and 1
.normalize01 <- function(x) { (x-min(x))/(max(x)-min(x)) }

###------------------------------
##-------------------
## Function to remove low variance features
## Function to remove low variance features
## @examples
## data(cars)
## removeZeroVar(cars, varCutoff=0)
removeZeroVar <- function(df, varCutoff=0, sort=TRUE)
{
  dfR <- apply(df,2, stats::var)
  dfR <- dfR[dfR>varCutoff]

  if(sort==TRUE)
  {
    dfR <- sort(dfR, decreasing = TRUE)
  }
  return(df[, names(dfR), drop=FALSE])
}




extractBetweenTags <- function(inVec, start.tag=0, end.tag=0)
{
  inVIndx= seq_along(inVec)
  stIndx = min(inVIndx[inVec!=start.tag])

  V2 = inVec[stIndx:length(inVec)]
  v2end = which(V2==end.tag)
  if(length(v2end)>0)
  {
    enIndx = min(v2end) -1
    enIndxR= enIndx + stIndx -1
  } else
  {enIndxR = length(inVec)}

  Vi = stIndx:enIndxR
  return(Vi)
}

