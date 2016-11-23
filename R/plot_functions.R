library(ggplot2)
#' Given a PDX object this will plot drug response curve
#' i.e. Time vs Volume
#' @param PDX a PDX object
#' @return ggplot object
#' @examples
#' plotDrugResponse(PDX)
#' @export



getExperimentDF <- function(expSlot, modelID, value)
{
  modelID = c(modelID)
  rtx = lapply(modelID, function(x){expSlot[[x]][[value]]})

  if(length(rtx)>1)
  {
    ##merge it, they are replacte
  }
  ## for now select 1
  rtx = rtx[[1]]

  return(rtx[,c("time", "volume")])
}


creatPlotDF <- function(expSlot, expList)
{

  col = rainbow(length(expDesX))
  px = list()

  for(I in 1:length(expDesX))
  {
    TrmodelID= expDesX[[I]][["treatment"]]
    tdf = getExperimentDF(expSlot, modelID = TrmodelID, value="data")

    CnmodelID= expDesX[[I]][["control"]]
    cdf = getExperimentDF(expSlot, modelID=CnmodelID, value="data")

    ##--- add batch ID -----------------
    tdf$batch = expDesX[[I]][["batch"]]
    cdf$batch = expDesX[[I]][["batch"]]

    ##---add drug name -----------------
    tdf$drug = expDesX[[I]][["drug.join.name"]]
    cdf$drug = expDesX[[I]][["drug.join.name"]]

    ## same color for one treatment and control
    tdf$color = col[I]
    cdf$color = col[I]

    tdf$pch = 16
    cdf$pch = 18

    px[[I]] = rbind(tdf, cdf)
  }

  plotData = do.call(rbind, px)
  return(plotData)
}



plotDrugResponse <- function(expSlot, expList)
{
  ## read the PDX object
  toPlt = readRDS("data/toPlot.Rda")
  expSlot=toPlt$expSlot
  expDesX=toPlt$expDesX

  plotData = creatPlotDF(expSlot, expList)

  ##---- plot the df in ggplot and return plot --------

}





