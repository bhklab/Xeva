library(ggplot2)

#' To plot drug response curve
#' @examples
#' data(pdxe)
#' plotDrugResponse(pdxe, drug="LEE011 + binimetinib")
#' plotDrugResponse(pdxe, drug="paclitaxel", drug.match.exact=TRUE, tumor.type="BRCA")
#' @param object The \code{XevaSet} to replace drug info in
#' @param drug Name of the drug
#' @return Updated \code{XevaSet}
setGeneric(name= "plotDrugResponse",
           def = function(object,
                          drug, drug.match.exact=TRUE,
                          tumor.type=NULL,
                          control=TRUE)
             {standardGeneric("plotDrugResponse")} )

#' @export
setMethod( f=plotDrugResponse,
           signature=c(object = "XevaSet"),
           definition=function(object,
                               drug, drug.match.exact=TRUE,
                               tumor.type=NULL,
                               control=TRUE)
           {

             allExpIds = selectModelIds(object, drug=drug, drug.match.exact=drug.match.exact,
                                     tumor.type=tumor.type)

             if(length(allExpIds)==0)
             {stop("No experiment present")}

             #allExpIds = allExpIds[1:3]

             #allExpModIds = mapModelSlotIds(object, id=allExpIds, id.name="experiment.id", map.to="model.id")
             #allBatchIds = c(unlist(sapply(allExpModIds$model.id, function(x) getBatchName(object, x))))

             allBatchIds = unlist(sapply(allExpIds, function(x) getBatchName(object, x)))

             expDesign2plt = lapply(allBatchIds, function(x) expDesignInfo(object)[[x]])

             batchNames = names(expDesign2plt)
             ##------ get color for each list element -----------
             colVec = grDevices::rainbow(length(batchNames))
             names(colVec) = batchNames

             DFlstx = list()
             for(en in batchNames)
             {

               dx = getTimeVarData(object, expDesign2plt[[en]], var = "volume")

               DFlstx[[en]] = .addColorPchLty(dx, colVec[en], treatment.lty="solid", control.lty="dashed")
             }

             DF = do.call(rbind, DFlstx)
             return(DF)
           } )



.addColorPchLty <- function(DFx, col, treatment.lty="solid", control.lty="dashed")
{
  lty = as.character(DFx$exp.type)
  lty[lty=="treatment"] = treatment.lty
  lty[lty=="control"] = control.lty
  DFx$lty = lty
  DFx$color = col
  return(DFx)
}

##----------------------------------------------------------------

#' @export
#' @import ggplot2
NewPlotFunction <- function(DF, drug.join.name)
{

  DF = readRDS("DATA-raw/toPlot_DF.Rda")
  drug.join.name = "paclitaxel"


  ##--- get range ------------------
  xr = range(DF$time, na.rm = TRUE)
  yr = range(DF$mean, DF$upper, DF$lower, na.rm = TRUE)

  ##--------
  #plt = ggplot() + geom_line(data = DF, aes(x=time, y=c2))
  plt = ggplot(DF, aes(x = time, y = mean))+ xlim(xr) + ylim(yr) + geom_blank()

  ## add one line
  dt = DF[DF$patient.id=="X-1004" & DF$exp.type=="treatment",]
  dc = DF[DF$patient.id=="X-1004" & DF$exp.type=="control",]

  addLine2Plot <- function(plt, dx, colour="blue")
  {
    plt = plt + geom_line(colour= colour, aes_string(x="time", y="mean"), dx)
    return(plt)
  }

  plt = addLine2Plot(plt, dt, colour="blue")
  plt = addLine2Plot(plt, dc, colour="red")


  plt = plt +  geom_line(colour="blue", aes_string(x="time", y="mean"), dt)
  plt = plt +  geom_line(colour="red", aes_string(x="time", y="mean"), dc)






  d<-data.frame(x=1:5, y1=1:5, y2=2:6)
  ggplot(d, aes(x)) +
    geom_line(aes(y=y1, colour="1")) +
    geom_line(aes(y=y2, colour="2")) +
    scale_colour_manual(values=c("red", "blue"))










  # create title objects
  drug_name <- unique(DF$drug.join.name[DF$drug.join.name != 'untreated'])
  title <- paste(length(DF$drug.join.name), 'Experiments for', drug_name)



  # create a unique variable combining patient.id and exp.type
  DF$patient.exp <- paste(DF$patient.id, DF$exp.type, sep = '_')

  # create basic plot object
  plot.1 <- ggplot(DF, aes(time, mean, group = patient.exp)) +
    xlab('Time') + ylab('Volume') + ggtitle(title)


    # if (!all(is.na(err.bars))) {
    #   # add error bars if not all of the error bar columns are NA
    #   plot.1 <- plot.1 + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper))
    #
    # }

  # plot add legend to plot.1 and points
  plot.final <- plot.1 + geom_line(aes(time, mean, colour = exp.type), data = DF, size = 0.7, alpha = 0.6) +
  scale_colour_manual(name = "", values=c("blue", "grey"), labels=c("Control", "Treatment"))

  # add final theme to plot.final
  plotObject <- plot.final + theme(panel.background=element_rect(fill="#F0F0F0"),
                                   plot.background=element_rect(fill="#F0F0F0"),
                                   panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
                                   legend.position="bottom",  plot.title=element_text(face="bold",colour="Black",size=10),
                                   axis.text.x=element_text(size=11,colour="#535353",face="bold"),
                                   axis.text.y=element_text(size=11,colour="#535353",face="bold"),
                                   axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
                                   axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5))

  return(plotObject)

}



##===============================================
.creatPlotDF <- function(object, DFx)
{
  rtx = list()
  for(I in 1:nrow(DFx))
  {
    expData = getExperiment(object, model.id= DFx[I, "model.id"])
    expData$batch.name= DFx[I, "batch.name"]
    expData$exp.type  = DFx[I, "exp.type"]
    rtx = .appendToList(rtx, expData)
  }
  rtz = do.call(rbind, rtx)

  return(rtz)
}

##================================================

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


plotDrugResponse_old <- function(expSlot, expList, diff.colors, patient.id = "X-1004")
{

  ##---- plot the df in ggplot and return plot --------
  dataX = .get_Data()

  ##---- get error bar information -------
  err.bars <- append(dataX$upper, dataX$lower)

  ##this seems a bug ----------------------------
  dataX = dataX[dataX$patient.id==patient.id,] ##select data for one patient
  ### The plot should give 2 lines: one control and one treatment
  ### but this gives only one line
  drug_name <- unique(dataX$drug.join.name[dataX$drug.join.name != 'untreated'])
  title <- paste(length(dataX$drug.join.name), 'Experiments for', drug_name)

  # create basic plot object
  plot.1 <- ggplot2::ggplot(dataX, ggplot2::aes(time, mean, group = exp.type)) +
  xlab('Time') +
  ylab('Volume') +
  ggtitle(title)
  # recode color for dataX (different color for treatment and control)
  if (diff.colors) {

    dataX <- getColor(dataX)

    if (!all(is.na(err.bars))) {
      # add error bars if not all of the error bar columns are NA
      plot.1 <- plot.1 + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper))

    }

    # plot add legend to plot.1 and points
    plot.final <- plot.1 + ggplot2::geom_line(ggplot2::aes(time, mean, colour = color), data = dataX, size = 0.7, alpha = 0.6) +
      ggplot2::scale_colour_manual(name = "", values=c("blue", "grey"), labels=c("Treatment", "Control")) +
      ggplot2::geom_point()


  } else {

    if (!all(is.na(err.bars))) {

      # add error bars if not all of the error bar columns are NA
      plot.1 <- plot.1 + ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper))

    }

    # plot add legend to plot.1 and points - on this one the no color distinction and
    # legend reflects pch, not color.
    plot.final <- plot.1 + ggplot2::geom_line(aes(time, mean, colour = color), data = dataX, size = 0.7, alpha = 0.6) +
      guides(color=FALSE) + ggplot2::geom_point(aes(shape = factor(lty))) +
      scale_shape_discrete(name = "", breaks=c("solid", "dashed"), labels=c("Treatment", "Control"))

  }

  # add final theme to plot.final
  plotObject <- plot.final + theme(panel.background=element_rect(fill="#F0F0F0"),
                                 plot.background=element_rect(fill="#F0F0F0"),
                                 panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
                                 legend.position="bottom",  plot.title=element_text(face="bold",colour="Black",size=10),
                                 axis.text.x=element_text(size=11,colour="#535353",face="bold"),
                                 axis.text.y=element_text(size=11,colour="#535353",face="bold"),
                                 axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
                                 axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5))

  return(plotObject)
}

# plot once with different colors, once with same colors, using pch to differentiate.
plotDrugResponse_old(diff.colors = T)
plotDrugResponse_old(diff.colors = F)






