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

    ##---add drug name -----------------
    tdf$drug = expDesX[[I]][["drug.join.name"]]
    cdf$drug = expDesX[[I]][["drug.join.name"]]

    ##--- add batch ID -----------------
    tdf$batch = expDesX[[I]][["batch"]]
    cdf$batch = expDesX[[I]][["batch"]]

    ##--- add type ID -----------------
    tdf$exp.type = "treatment"
    cdf$exp.type = "control"

    ## same color for all treatment and control
    tdf$color = col[I]
    cdf$color = col[I]

    tdf$pch = 16
    cdf$pch = 18

    px[[I]] = rbind(tdf, cdf)
  }

  plotData = do.call(rbind, px)
  return(plotData)
}


.get_Data <- function()
{
  if(1==2){
  ## read the PDX object
  toPlt = readRDS("DATA-raw/toPlot.rdata")
  expSlot=toPlt$expSlot
  expDesX=toPlt$expDesX
  plotData = creatPlotDF(expSlot, expList)
  plotData2 = plotData[plotData$batch=="X-2094",]
  SE = (plotData2$volume/sum(plotData2$volume))*1000
  plotData2$SE.lower = plotData2$volume - SE
  plotData2$SE.upper = plotData2$volume + SE
  saveRDS(list(plotData, plotData2), file = "DATA-raw/toPlot_DF.Rda")
  }
  plotData = readRDS("DATA-raw/toPlot_DF.Rda")
  return(plotData)
}

plotDrugResponse <- function(expSlot, expList)
{

  dataX = .get_Data()
  plotData = dataX[[1]] ## data without error bars
  plotData2= dataX[[2]] ## data with error bars

  ##this seems a bug ----------------------------
  plotData = plotData[plotData$batch=="X-2094",] ##select data for one patient
  ### The plot should give 2 lines: one control and one treatment
  ### but this gives only one line

  title <- paste(length(unique(plotData$batch)), 'Experiments for', unique(plotData$drug))

  ##---- plot the df in ggplot and return plot --------
  p1 <- ggplot(plotData, aes(time, volume, group = exp.type)) + xlab('Time') + ylab('Volume') + ggtitle(title)

  p_line <- p1 + geom_line(aes(time, volume, colour = color), data = plotData, size = 0.7, alpha = 0.6) + guides(color=FALSE)

  p_point <- p_line + geom_point(aes(shape = factor(pch)))

  p_legend <- p_point + scale_shape_discrete(name = "", breaks=c(16, 18), labels=c("Treatment", "Control"))

  plotObject <- p_legend + theme(panel.background=element_rect(fill="#F0F0F0"),
                                 plot.background=element_rect(fill="#F0F0F0"),
                                 panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
                                 legend.position="bottom",  plot.title=element_text(face="bold",colour="Black",size=10),
                                 axis.text.x=element_text(size=11,colour="#535353",face="bold"),
                                 axis.text.y=element_text(size=11,colour="#535353",face="bold"),
                                 axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
                                 axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5))

  return(plotObject)
}






