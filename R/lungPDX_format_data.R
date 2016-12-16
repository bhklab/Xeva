
if(1==2){
xlFile = "data/PHLC110 Ludo (Every Mon) May01.15 final.xlsx"

library(XLConnect)
lx <- readWorksheetFromFile(xlFile, sheet="Treatment", header = FALSE )
date = lx[1,]

}
