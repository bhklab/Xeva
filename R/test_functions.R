#'
#' test_mResest <- function()
#' {
#'
#'   if(1==2)
#'   {
#'
#'     if(1==2)
#'     {
#'
#'       ########----------------------------------------------------------------------------
#'       #' @import ggplot2
#'       ##-------------------------------------------------------
#'       #library(ggplot2)
#'       #' To plot drug response curve
#'       #' @examples
#'       #' data(pdxe)
#'       #' plotDrugResponse(pdxe, drug="LEE011 + binimetinib")
#'       #' plotDrugResponse(pdxe, drug="paclitaxel", drug.match.exact=TRUE, tumor.type="BRCA")
#'       #' @param object The \code{XevaSet} to replace drug info in
#'       #' @param drug Name of the drug
#'       #' @return Updated \code{XevaSet}
#'       setGeneric(name= "plotDrugResponse",
#'                  def = function(object,
#'                                 drug, drug.match.exact=TRUE,
#'                                 tumor.type=NULL,
#'                                 control=TRUE)
#'                  {standardGeneric("plotDrugResponse")} )
#'
#'       #' @export
#'       setMethod( f=plotDrugResponse,
#'                  signature=c(object = "XevaSet"),
#'                  definition=function(object,
#'                                      drug, drug.match.exact=TRUE,
#'                                      tumor.type=NULL,
#'                                      control=TRUE)
#'                  {
#'
#'                    allExpIds = selectModelIds(object, drug=drug, drug.match.exact=drug.match.exact,
#'                                               tumor.type=tumor.type)
#'
#'                    if(length(allExpIds)==0)
#'                    {stop("No experiment present")}
#'
#'                    #allExpIds = allExpIds[1:3]
#'
#'                    #allExpModIds = mapModelSlotIds(object, id=allExpIds, id.name="experiment.id", map.to="model.id")
#'                    #allBatchIds = c(unlist(sapply(allExpModIds$model.id, function(x) getBatchName(object, x))))
#'
#'                    allBatchIds = unlist(sapply(allExpIds, function(x) getBatchName(object, x)))
#'
#'                    expDesign2plt = lapply(allBatchIds, function(x) expDesignInfo(object)[[x]])
#'
#'                    batchNames = names(expDesign2plt)
#'                    ##------ get color for each list element -----------
#'                    colVec = grDevices::rainbow(length(batchNames))
#'                    names(colVec) = batchNames
#'
#'                    DFlstx = list()
#'                    for(en in batchNames)
#'                    {
#'
#'                      dx = getTimeVarData(object, expDesign2plt[[en]], var = "volume")
#'
#'                      DFlstx[[en]] = .addColorPchLty(dx, colVec[en], treatment.lty="solid", control.lty="dashed")
#'                    }
#'
#'                    DF = do.call(rbind, DFlstx)
#'                    return(DF)
#'                  } )
#'
#'
#'
#'       .addColorPchLty <- function(DFx, col, treatment.lty="solid", control.lty="dashed")
#'       {
#'         lty = as.character(DFx$exp.type)
#'         lty[lty=="treatment"] = treatment.lty
#'         lty[lty=="control"] = control.lty
#'         DFx$lty = lty
#'         DFx$color = col
#'         return(DFx)
#'       }
#'
#'       ##----------------------------------------------------------------
#'       .takeLogOfDF <- function(DF, colN, removeNan=TRUE)
#'       {
#'         DF[, colN] = log(DF[, colN])
#'         if(removeNan==TRUE){
#'           DF = DF[is.finite(DF[, colN]), ] }
#'         return(DF)
#'       }
#'
#'
#'
#'
#'
#'
#'       ##===============================================
#'       .creatPlotDF <- function(object, DFx)
#'       {
#'         rtx = list()
#'         for(I in 1:nrow(DFx))
#'         {
#'           expData = getExperiment(object, model.id= DFx[I, "model.id"])
#'           expData$batch.name= DFx[I, "batch.name"]
#'           expData$exp.type  = DFx[I, "exp.type"]
#'           rtx = .appendToList(rtx, expData)
#'         }
#'         rtz = do.call(rbind, rtx)
#'
#'         return(rtz)
#'       }
#'
#'       ##================================================
#'       plotDrugResponse_old <- function(expSlot, expList, diff.colors, patient.id = "X-1004")
#'       {
#'
#'         ##---- plot the df in ggplot and return plot --------
#'         dataX = .get_Data()
#'
#'         ##---- get error bar information -------
#'         err.bars <- append(dataX$upper, dataX$lower)
#'
#'         ##this seems a bug ----------------------------
#'         dataX = dataX[dataX$patient.id==patient.id,] ##select data for one patient
#'         ### The plot should give 2 lines: one control and one treatment
#'         ### but this gives only one line
#'         drug_name <- unique(dataX$drug.join.name[dataX$drug.join.name != 'untreated'])
#'         title <- paste(length(dataX$drug.join.name), 'Experiments for', drug_name)
#'
#'         # create basic plot object
#'         plot.1 <- ggplot2::ggplot(dataX, ggplot2::aes(time, mean, group = exp.type)) +
#'           xlab('Time') +
#'           ylab('Volume') +
#'           ggtitle(title)
#'         # recode color for dataX (different color for treatment and control)
#'         if (diff.colors) {
#'
#'           dataX <- getColor(dataX)
#'
#'           if (!all(is.na(err.bars))) {
#'             # add error bars if not all of the error bar columns are NA
#'             plot.1 <- plot.1 + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper))
#'
#'           }
#'
#'           # plot add legend to plot.1 and points
#'           plot.final <- plot.1 + ggplot2::geom_line(ggplot2::aes(time, mean, colour = color), data = dataX, size = 0.7, alpha = 0.6) +
#'             ggplot2::scale_colour_manual(name = "", values=c("blue", "grey"), labels=c("Treatment", "Control")) +
#'             ggplot2::geom_point()
#'
#'
#'         } else {
#'
#'           if (!all(is.na(err.bars))) {
#'
#'             # add error bars if not all of the error bar columns are NA
#'             plot.1 <- plot.1 + ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper))
#'
#'           }
#'
#'           # plot add legend to plot.1 and points - on this one the no color distinction and
#'           # legend reflects pch, not color.
#'           plot.final <- plot.1 + ggplot2::geom_line(aes(time, mean, colour = color), data = dataX, size = 0.7, alpha = 0.6) +
#'             guides(color=FALSE) + ggplot2::geom_point(aes(shape = factor(lty))) +
#'             scale_shape_discrete(name = "", breaks=c("solid", "dashed"), labels=c("Treatment", "Control"))
#'
#'         }
#'
#'         # add final theme to plot.final
#'         plotObject <- plot.final + theme(panel.background=element_rect(fill="#F0F0F0"),
#'                                          plot.background=element_rect(fill="#F0F0F0"),
#'                                          panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
#'                                          legend.position="bottom",  plot.title=element_text(face="bold",colour="Black",size=10),
#'                                          axis.text.x=element_text(size=11,colour="#535353",face="bold"),
#'                                          axis.text.y=element_text(size=11,colour="#535353",face="bold"),
#'                                          axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
#'                                          axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5))
#'
#'         return(plotObject)
#'       }
#'
#'     }
#'
#'   }
#'
#' }
