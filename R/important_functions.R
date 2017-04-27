##-- this will convert list to data.fram -------
<<<<<<< HEAD
# v <- 1:3; names(v) <- LETTERS[1:3]
# inLst <- list(x=v, y=v*10)
# .convertListToDataFram(inLst)
=======
#' v <- 1:3; names(v) <- LETTERS[1:3]
#' inLst <- list(x=v, y=v*10)
#' .convertListToDataFram(inLst)
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
## @export
.convertListToDataFram <- function(inLst)
{
  name.Rows = names(inLst)
  if(is.null(name.Rows))
  { name.Rows= 1:length(inLst)}

  name.Cols = unique(unlist(lapply(inLst, names)))
  if(is.null(name.Cols))
  { name.Cols = 1:max(unlist(lapply(inLst, length))) }

  rd = data.frame(matrix(NA, nrow = length(name.Rows), ncol = length(name.Cols)))
  rownames(rd) = name.Rows
  colnames(rd) = name.Cols
  for(rx in name.Rows)
  {
    v = inLst[[rx]]
    if(is.null(names(v))){names(v) = 1:length(v)}
    rd[rx, names(v)] = unlist(v, use.names = TRUE)[names(v)]
  }
  return(rd)
}

###-------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------
##------ given a list of datafram it will give a new datafram --------------------------
<<<<<<< HEAD
# x <- data.frame(a=1:3, b= LETTERS[1:3])
# y <- data.frame(a=(1:3)*10, b= letters[1:3])
# inLst <- list(x=x, y=y)
# .rbindListOfDataframs(inLst)
=======
#' x <- data.frame(a=1:3, b= LETTERS[1:3])
#' y <- data.frame(a=(1:3)*10, b= letters[1:3])
#' inLst <- list(x=x, y=y)
#' .rbindListOfDataframs(inLst)
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
## @export
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

###-------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------

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
      if(length(v)==0){v=NA}
      vz=NULL
      if(collapse =="mean")
      {
        vz = mean(as.numeric(v[!is.na(v)]))
      } else if (collapse =="median")
      {
        vz = median(as.numeric(v[!is.na(v)]))
      } else
      {
        if(length(v)>1)
        {
          vz= pasteWithoutNA(v, collapse = collapse)
        } else
        { vz=v }

        if(!is.na(vz) & vz==""){vz=NA}
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
##---------------Convert data.frame columns from factors to characters ---------------------
## This will convert all factor type columns to character columns
.factor2char <- function(df)
{
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  return(df)
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



#' Symmetric set differences
#'
#' \code{symmetricSetDiff} give symmetric set differences. Symmetric set difference (disjunctive union), of two sets is the set of elements which are not in their intersection.
#'
#' @param a vector
#' @param b vector
#'
#' @return  Returns a vector
#'
#' @examples
#' a <- c(1, 2, 3)
#' b <- c(2, 3, 4)
#' symmetricSetDiff(a, b)
#' @export
symmetricSetDiff <- function(a,b){ unique(c(setdiff(a,b), setdiff(b,a))) }



##-------------------
#' Function to print
#'
#' \code{printf} is shortcut for cat(sprintf())
#' @examples
#' printf("Hello %s\n", "world" )
#' @export
printf <- function(...) cat(sprintf(...))




##-------------------
#' Function to print data.frame in massage
#'
#' \code{printAndCapture} prints data.frame in stop or warning functions
#' @examples
#' df <- data.frame(a=1:5, b=11:15)
#' msg <- sprintf("data frame is:\n%s", printAndCapture(df))
#' warning(msg)
#' @export
printAndCapture <- function(x)
{
  paste(capture.output(print(x)), collapse = "\n")
}


###----------------------------
##Normalize a vector between 0 and 1
.normalize01 <- function(x) { (x-min(x))/(max(x)-min(x)) }

###------------------------------
##-------------------
#' Function to remove low variance features
#'
#' Function to remove low variance features
#' @examples
#' data(cars)
#' removeZeroVar(cars, varCutoff=0)
#' @export
removeZeroVar <- function(df, varCutoff=0, sort=TRUE)
{
  dfR <- apply(df,2, var)
  dfR <- dfR[dfR>varCutoff]

  if(sort==TRUE)
  {
    dfR <- sort(dfR, decreasing = TRUE)
  }
  return(df[, names(dfR), drop=FALSE])
}


