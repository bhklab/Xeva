
xlFile = "data/PHLC110 Ludo (Every Mon) May01.15 final.xlsx"
##-- read the Lung xls
#library(openxlsx)
#sheetLG = getSheetNames(xlFile)
#lc1 <- read.xlsx(xlFile, sheet = "Treatment", colNames = FALSE,
#                 rowNames = FALSE,detectDates = TRUE)

library(XLConnect)
lc1 <- readWorksheetFromFile(xlFile, sheet="Treatment", header = FALSE )


