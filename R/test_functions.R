
test_mResest <- function()
{

  library(XLConnect)
  xlFile = "DATA-raw/Geo_PCT_curve_matrics.xlsx"
  gr = readWorksheetFromFile(xlFile, sheet="PCT curve metrics", header = TRUE )
  gr$Treatment = gsub('figitumumab\"', 'figitumumab', gr$Treatment)
  gr = BBmisc::sortByCol(gr , "Day_BestResponse", asc = rep(TRUE, 1))


  data("pdxe")
  setmRECIST(pdxe)<- setmRECIST(pdxe)

  dfx = lapply(pdxe@experiment, function(x)
    {
    rx = c(x$batch, x$drug$join.name, x$mRECIST, NA)
    names(rx) = c("model", "drug","xeno.mR", "org.mR")
    gmr = gr[gr$Model==x$batch & gr$Treatment==x$drug$join.name, "ResponseCategory"]
    if(length(gmr)>0){ rx["org.mR"] = gmr}
    return(rx)
    } )

  z= .convertListToDataFram(dfx)
  z[is.na(z$xeno.mR),]
  z[is.na(z$org.mR),]

  org.mR1st = sapply(strsplit(z$org.mR,"-"), function(x) x[1])
  z[,"org.mR1st"] = org.mR1st
  k = z[z$xeno.mR!=z$org.mR1st,]
  dim(k)

  pdxe@experiment[["X-6047.18.paclitaxel"]]

}
