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



  ## loop through experiement data and find all model ids with specified drugname
  slotModelID <- list()
  for (i in 1:length(expSlot)) {
    slotModelID[[i]] <- expSlot[[i]]$model.id[drugName == expSlot[[i]]$drug$join.name]
  }
  ## get vector of model ids in data frame
  slotModelID <- as.data.frame(unlist(slotModelID))


  ## loop through experimental design information and get ids for treatment and control groups with
  ## specified drug name
  desXTreatment <- list()
  desXControl <- list()
  for (i in 1:length(expDesX)) {
    desXTreatment[[i]] <- expDesX[[i]]$treatment[drugName == expDesX[[i]]$drug.join.name]
    desXControl[[i]] <- expDesX[[i]]$control[drugName == expDesX[[i]]$drug.join.name]
  }

  ## get vector of treatment and control group ids in data frame
  desXTreatment <- as.data.frame(unlist(desXTreatment))
  desXControl <- as.data.frame(unlist(desXControl))


  ## find the intersection of both treatment and control with model ids for each drugname
  treatment <- intersect(slotModelID$`unlist(slotModelID)`, desXTreatment$`unlist(desXTreatment)`)

  #control <- intersect(slotModelID$`unlist(slotModelID)`, desXControl$`unlist(desXTreatment)`)
  control <- desXControl$`unlist(desXControl)`


  treatmentCol = "darkred"
  controlCol = "green"

  ## get Experimental data associated with treatment ids. add columns to that data to identify
  ## the data frame as treatment with a corresponding color (data in long format for ggplot)
  plotTreatmentData <- list()
  for (id in treatment) {
    data <- expSlot[[id]]$data
    data$id <- id
    data$type <- 'treatment'
    data$col <- treatmentCol
    plotTreatmentData[[id]] <- data
  }

  ## collapse list to create data frame
  plotTreatmentData <- do.call(rbind, plotTreatmentData)

  ## get Experimental data associated with control ids. add columns to that data to identify
  ## the data frame as control with a corresponding color (data in long format for ggplot)
  plotControlData <- list()
  for (id in control) {
    data <- expSlot[[id]]$data
    data$id <- id
    data$type <- 'control'
    data$col <- controlCol
    plotControlData[[id]] <- data
  }

  ## collapse list to create data frame
  plotControlData <- do.call(rbind, plotControlData)

  ## combine treatement and control data
  plotData <- rbind(plotControlData, plotTreatmentData)

  ## get number of patients in graph
  n_patients <- length(unique(plotData$id))

  ## assign a color vector for ggplot2 - to distinguish treatment vs control
  col <- plotData$col


  ## plot data
  plotObject <- ggplot(plotData, aes(time, volume, group = id)) +
    geom_line(aes(x = time, y = volume, colour = type), data = plotData, size = 0.7, alpha = 0.6) +
    scale_color_manual(name = '',labels = c('Control','Treatment' ),
                       values=c("treatment"=treatmentCol, "control" = controlCol))  +
    geom_point(aes(x = time, y = volume), size = 0.7, colour = plotData$col, alpha = 1) +
    xlab('Time') +
    ylab('Volume')
  plotObject = plotObject+ ggtitle(paste("Showing", n_patients, "Patients", sep = ' '))

  plotObject = plotObject + theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.background = element_blank(),
                                  axis.line = element_line(colour = "black"),
                                  legend.background = element_rect(colour = "black"),
                                  plot.title = element_text(hjust = 0.5))

  #+
    # theme(panel.background=element_rect(fill="white"),
    #       plot.background=element_blank(), #element_rect(fill="#F0F0F0"),
    #       #panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
    #       legend.background = element_blank(),#element_rect(fill = 'lightgrey'),
    #       legend.position = 'bottom',
    #       axis.text.x=element_text(size=11,colour="#535353",face="bold"),
    #       axis.text.y=element_text(size=11,colour="#535353",face="bold"),
    #       axis.title.y=element_text(size=11,colour="#535353",face="bold"),
    #       axis.title.x=element_text(size=11,colour="#535353",face="bold"),
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank())

  return(plotObject)


}

plotDrugResponse(drugName = 'LCL161 + paclitaxel')

