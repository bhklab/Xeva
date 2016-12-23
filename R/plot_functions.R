library(ggplot2)

#' To plot drug response curve
#' @examples
#' data(pdxe)
#' plotDrugResponse(pdxe, drug="LEE011 + binimetinib")
#' @param object The \code{XevaSet} to replace drug info in
#' @param drug Name of the drug
#' @return Updated \code{XevaSet}
setGeneric(name= "plotDrugResponse", def = function(object, drug) {standardGeneric("plotDrugResponse")} )

#' @export
setMethod( f=plotDrugResponse,
           signature=c(object = "XevaSet"),
           definition=function(object, drug)
           {
             #object@expDesign = value
             return(object)
           } )



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

# function to assign a seperate color
getColor <- function(plotData)
{
  plotData$color <- ifelse(plotData$exp.type == 'treatment', 'blue', 'grey')
  return(plotData)
}


plotDrugResponse_old <- function(expSlot, expList, err.bars, colors = 'different')
{

  ##---- plot the df in ggplot and return plot --------
  dataX = .get_Data()
  plotData = dataX[[1]] ## data without error bars
  plotData2= dataX[[2]] ## data with error bars

  ##this seems a bug ----------------------------
  plotData = plotData[plotData$batch=="X-2094",] ##select data for one patient
  ### The plot should give 2 lines: one control and one treatment
  ### but this gives only one line
  title <- paste(length(unique(plotData$batch)), 'Experiments for', unique(plotData$drug))

  # recode color for plotData (different color for treatment and control)
  if (colors == 'different') {
    plotData <- getColor(plotData)
    plotData2 <- getColor(plotData2)

    p1 <- ggplot(plotData2, aes(time, volume, group = exp.type)) + xlab('Time') + ylab('Volume') + ggtitle(title)

    if (err.bars) {
      # scale_fill_manual(name="Bar",values=cols, guide="none")
      p1 <- p1 + geom_errorbar(aes(ymin = SE.lower, ymax = SE.upper))

    }

    p_line <- p1 + geom_line(aes(time, volume, colour = color), data = plotData, size = 0.7, alpha = 0.6)
    # + guides(color=FALSE)
    p_legend <- p_line + scale_colour_manual(name = "", values=c("blue", "black"), labels=c("Treatment", "Control"))

    p_point <- p_line + geom_point(aes(shape = factor(pch)))


  } else {


    p1 <- ggplot(plotData2, aes(time, volume, group = exp.type)) + xlab('Time') + ylab('Volume') + ggtitle(title)

    if (err.bars) {

      p1 <- ggplot(plotData2, aes(time, volume, group = exp.type )) + xlab('Time') + ylab('Volume') + ggtitle(title)
      p1 <- p1 + geom_errorbar(aes(ymin = SE.lower, ymax = SE.upper))

    }

    p_line <- p1 + geom_line(aes(time, volume, colour = color), data = plotData, size = 0.7, alpha = 0.6) + guides(color=FALSE)

    p_point <- p_line + geom_point(aes(shape = factor(pch)))

    p_legend <- p_point + scale_shape_discrete(name = "", breaks=c(16, 18), labels=c("Treatment", "Control"))

  }

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






