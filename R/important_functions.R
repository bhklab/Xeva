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


.castDataFram <- function(df, row.var, col.var, value, collapse = ";")
{
  uqR = unique(as.character(df[,row.var]))
  uqC = unique(as.character(df[,col.var]))
  dfx = data.frame( matrix(NA, nrow = length(uqR), ncol = length(uqC) ))
  rownames(dfx) = uqR ; colnames(dfx) = uqC

  for(r in rownames(dfx))
  {
    for(c in colnames(dfx))
    {
      v = df[df[,row.var]==r & df[,col.var]==c, value]

      if(collapse =="mean") {
         vz = mean(as.numeric(v[!is.na(v)]))
      } else if (collapse =="median"){
        vz = median(as.numeric(v[!is.na(v)]))
      } else {
        vz= pasteWithoutNA(v, collapse = collapse)
        if(vz==""){vz=NA}
      }
      if(is.nan(vz)){vz = NA}
      dfx[r,c]=vz
    }
  }
  return(dfx)
}


.appendToList <- function(in.list, value)
{
  in.list[[length(in.list)+1]]=value
  return(in.list)
}


#' @export
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

##------------------------------------------------------------------------------------------
##---------------remove NA col--------------------------------------------------------------
.removeNAcol <- function(df)
{ return(df[, !apply(is.na(df), 2, all)]) }

##------------------------------------------------------------------------------------------
##---------------reorder column ------------------------------------------------------------
.reorderCol <- function(df, columnName, newIndx)
{
  OtherCN = colnames(df)[colnames(df)!=columnName]
  newCN = append(OtherCN, columnName, after= (newIndx-1) )
  return(df[,newCN])
}

##------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------
#' paste a vector elements togather while removing NA
#'
#' \code{pasteWithoutNA} paste a vector elements togather while removing NA
#'
#' @param L A vector with values and NA
#' @param collapse Collapse string default " + "
#'
#' @return  Returns an string with vector values paste togather
#'
#' @examples
#' L = c("A", NA, "B", NA, NA, "C")
#' pasteWithoutNA(L, collapse = " + ")
#' @export
pasteWithoutNA <- function(L, collapse = " + "){paste(L[!is.na(L)], collapse = collapse)}

##------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------
#' paste a data.frame columns togather while removing NA
#'
#' \code{pasteColTogather} paste a data.frame columns togather while removing NA
#'
#' @param df A data.frame
#' @param collapse Collapse string default " + "
#'
#' @return  Returns an vector of strings where column values paste togather
#'
#' @examples
#' df = data.frame(x= 1:6, y = c("A", NA, "B", NA, NA, "C"))
#' pasteColTogather(df, collapse = " + ")
#' @export
pasteColTogather <- function(df, collapse = " + ")
{
  apply(df, 1, pasteWithoutNA, collapse =collapse)
}



##------------------------------------------------------------------------
##--- this will creat empty theme for ggplot -----------------------------
.ggplotEmptyTheme <- function(plt){
  plt +  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        panel.background = ggplot2::element_blank(),
                        legend.key=element_blank(),   ##removes legend background
                        axis.line = ggplot2::element_line(colour = "black"))
}

