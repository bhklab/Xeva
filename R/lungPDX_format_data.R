if(1==2){

slideFunct <- function(data, window, step)
{
  windowX = window-1
  maxIndx = length(data)-windowX
  spots <- seq(from=1, to=maxIndx, by=step)
  rtx = list()
  for(i in 1:length(spots))
  {
    rtx[[i]] = data[spots[i]:(spots[i]+windowX)]
  }
  return(rtx)
}

.appendToList <- function(in.list, value)
{
  in.list[[length(in.list)+1]]=value
  return(in.list)
}

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

##-------------------------------------------------------------------------
##--------------------for one Row -----------------------------------------
getDFfor1Row <- function(rw, dataHeader, expDate)
{
  #rw = dataX[1,]
  weightIndx = grep("Weight",dataHeader)
  spltRowIndx = slideFunct(c(weightIndx, (length(rw)+1)), window=2, step=1)
  dfl= list()
  for(I in spltRowIndx)
  {
    di = rw[I[1]:(I[2]-1)]
    diHead = dataHeader[I[1]:(I[2]-1)]
    colnames(di) = diHead

    dix = di[, !apply(is.na(di), 2, all)]
    if(ncol(dix)>0)
    {
      diDate = expDate[I[1]]
      dix[, "Date"] = diDate
      #print(dix); cat("\n")
      dfl = .appendToList(dfl, dix)
    }
  }
  dfx = .rbindListOfDataframs(dfl)
  return(dfx)
}
##-------------------------------------------------------------------------
##-------------------------------------------------------------------------

library(XLConnect)
##----- read a file -----
formatOnePDXFile <- function(flx)
{
  #flx = fileLst[1]
  cat(sprintf("##=== processing  %s", flx))
  lx <- readWorksheetFromFile(flx, sheet="Treatment", header = FALSE )

  dataStrIndx = grep("Animal", lx[,1])
  dataEndIndx = grep("Extra Information", lx[,1])
  if(length(dataStrIndx)!=1 | length(dataEndIndx)!=1)
  {
    msg= sprintf("Error in data start or end index of file\n%s", flx)
    stop(msg)
  }
  dataHeader = lx[dataStrIndx,]
  data = lx[(dataStrIndx[1] + 1): (dataEndIndx[1]-1),]

  ##--- remove AVERAGE rows ------------------------------
  AVERAGEcols = grep("AVERAGE", data[,1])
  ##----- remove rows which have first cell NA -----------
  rwWithNa = which(is.na(data[,1]))
  rw2Remove = c(AVERAGEcols, rwWithNa)
  newRow = setdiff(1:nrow(data), rw2Remove)
  dataX = data[newRow, ]

  ##---------get other information ---------------------
  expDate = lx[1,]

  expID = lx[2,1]
  expID = gsub(" \\(", "_", expID)
  expID = gsub("\\)", "", expID)

  Donor = lx[2,2]
  Donor = gsub("Donor:", "", Donor)

  DFlist = list()
  for(rwi in 1:nrow(dataX))
  {
    rw = dataX[rwi,]
    dfi = getDFfor1Row(rw, dataHeader, expDate)

    model.id = sprintf("%s.%s.%s", expID, rw[,1], rw[,2])
    dfi$model.id = model.id

    dfi$drug = rw[,3]
    dfi$batch= expID
    dfi$Donor = Donor

    DFlist=.appendToList(DFlist, dfi)
  }

  DF = .rbindListOfDataframs(DFlist)
  return(DF)
}


dirPath = "/home/arvind/CXP/XG/Data/Celine_Data/Nhu-An_21Dec2016"
fileLst = list.files(path = dirPath, pattern = NULL, full.names = TRUE)





}











