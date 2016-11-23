if(1==2){

library(XLConnect)
processFile <- function(xlFile)
{
bc <- readWorksheetFromFile(xlFile, sheet="Data", header = FALSE )


getBlocks = which(bc[,1] == "Treatment")
getBlocks = c(getBlocks, length(bc[,1])+1)

infoBlock = list()

i=1; strtIndx = 1
for (endIndx in getBlocks)
{
  infoBlock[[i]] = bc[strtIndx:(endIndx-1),]
  i = i+1
  strtIndx = endIndx
}



##This is the model (REF001BL) and the passage P1. Could separate these into 2 fields?
mpx = strsplit(infoBlock[[1]][1,1], " ")[[1]]
MODPSG = list(model= paste(mpx[1:length(mpx)-1], collapse = "_"),
              passage = mpx[length(mpx)])


if(infoBlock[[1]][2,2] == "Date")
  { expDates = infoBlock[[1]][2,]
  } else{
    msg = sprintf("cell %s%d should be Date", LETTERS[2], 2)
    stop(msg)}


##==== process a block ============

processOneExperment <- function(expVal, dateVec)
{
  expVal= as.vector(t(expVal))
  dateVec=as.vector(t(dateVec))
  widIndx = seq(1,length(expVal),3)
  lenIndx = seq(2,length(expVal),3)
  volIndx = seq(3,length(expVal),3)
  rtx = data.frame(width = expVal[widIndx],
                   length= expVal[lenIndx],
                   volume= expVal[volIndx],
                   stringsAsFactors = FALSE)

  rtx$date = dateVec[volIndx]
  rtx$date = as.Date(rtx$date)
  rtx$time = rtx$date - rtx$date[1]
  return(rtx)
}

##-------------------------------------------------------
##--- creat a model Id for each mouse used in experement
creatModelId <- function(mouseID, RowName)
{
  mid = sprintf("%s_%s_%s_%s", MODPSG$model,MODPSG$passage, mouseID, RowName)
  return(mid)
}

##==================================
allExp = list(); indxAllExp = 1

expNdf = data.frame()
for(blokID in 2:length(infoBlock))
{
  #blok = infoBlock[[4]]  ### for test
  blok = infoBlock[[blokID]]
  if(blok[1,1]=="Treatment")
    { drug = blok[2,1]
    } else {
    msg = sprintf("cell %s%d should be Treatment", LETTERS[4],1)
    stop(msg)}

  for (expri in 2:6) ## max 5 experements in each block
  {
    if(!is.na(blok[expri,2]))
    {
      expVal = blok[expri,3:dim(blok)[2]]
      dateVec = expDates[3:length(expDates)]
      expN = processOneExperment(expVal, dateVec)
      expN$formula = "width*width*length/2"

      expN$model.id = creatModelId( mouseID=blok[expri,2], RowName=rownames(blok[expri,]) )
      expN$drug  = drug
      expN$batch = xlFile
      if(drug=="H20")
      {expN$type  = "Control"} else { expN$type  = "Treatment"}

      expNdf = rbind(expNdf, expN)
      allExp[[indxAllExp]] = expN
      indxAllExp = indxAllExp + 1
    }
  }
}

##--- remove empty rows from df ---------------------------
chkCol = c("width", "length", "volume")
valuesToAvoid = c(NA,0)

rwTest = sapply(chkCol, function(x){ expNdf[,x] %in% valuesToAvoid})
rw2Remove = apply(rwTest, 1, all)
expNdf = expNdf[!rw2Remove, ]

##----- split drug ---------------
source("R/important_functions.R")
#dg  = unique(expNdf$drug)
dg  = expNdf$drug
dgx = strsplit(dg, "\\+")
dgy = convertListToDataFram(dgx)
dgy = as.data.frame(apply(dgy, 2, trimws), stringsAsFactors = FALSE)
class(dgy)

colnames(dgy) = paste0("drug.", 1:dim(dgy)[2])

dosx = data.frame(matrix(NA, nrow = dim(dgy)[1], ncol = dim(dgy)[2]))
colnames(dosx) = paste0("dose.", 1:dim(dosx)[2])

##---- manually put dose values -------------------------
dosx[dgy=="Intermittent TTK 6mg"] = 6
dgy[ dgy=="Intermittent TTK 6mg"] = "Intermittent TTK"

dosx[dgy=="Taxol 40mg"] = 40
dgy[ dgy=="Taxol 40mg"] = "Taxol"

dosx[dgy=="Taxol 40mg/kg"] = 40
dgy[ dgy=="Taxol 40mg/kg"] = "Taxol"
##-------------------------------------------------------

expNdf = cbind(expNdf, dgy, dosx)
return(expNdf)

}

#xlFile = "data/REF001BLP1 (H2O_945_Intermittent TTK60mg+Taxol40mg_Taxol_Carbo) June 10_16.xlsx"
xlFile = "data/NOTCHB01P5-(PBSctrl_DPctrl,Fluva_DP,PBSctrl_DP,Fluva_DPctrl)-Nov 2_16.xlsx"

expNdf = processFile(xlFile)
##save datafram for test
saveRDS(expNdf, file="data/DC_pdxDf_test_object2.Rda")





if(1==2){

##==== Trim each experiment =====================
getIndex <- function(inVec, indxOf)
{
  if(is.na(indxOf)){return(which(is.na(inVec)))}
  return(which(inVec==indxOf))
}

trimExperimentData <- function(IAllExp)
{
  ##--- remove all rows for which all ("width", "length", "volume") == NA or 0
  desiredCols = c("width", "length", "volume")
  valuesToAvoid = c(NA,0)
  toKeep = c()
  for(I in 1:max(sapply(desiredCols, function(x)length(IAllExp[[x]]))))
  {
    valChk = sapply(desiredCols, function(x){IAllExp[[x]][I] %in% valuesToAvoid})
    if(all(valChk==TRUE)){next}
    toKeep = c(toKeep, I)
  }

  IAllExpOut = IAllExp
  for(j in c(desiredCols, "date", "time"))
  { IAllExpOut[[j]] = IAllExp[[j]][toKeep]}


  lenChk = sapply(desiredCols, function(x){length(IAllExpOut[[x]])})

  if(all(lenChk==0)){IAllExpOut =NULL}
  return(IAllExpOut)
}


allExp = lapply(allExp, function(IAllExp){trimExperimentData(IAllExp)})
##-- remove null elements
allExp = allExp[lapply(allExp,length)>0]

saveRDS(allExp, file="data/PDX_test_object.Rda")

colors <- rainbow(length(allExp));
plot(c(1,100), c(1,1000),type="n", xlab="Time", ylab="Tumor volume" )
colIndx = 1
for(i in allExp)
{
  lines(i$time, i$volume, type = "o", pch=19, col=colors[colIndx])
  colIndx = colIndx+1
  print(i$drug)
}

drgName = sapply(allExp, function(x)x$drug)

legend("topleft", legend=drgName, cex=0.8, col=colors, pch=19, lty=1, title="Drug")


}







}
