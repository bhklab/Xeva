
#' Extract controls for a given model.id
#'
#' @examples
#' data(pdxe)
#' # extract controls for a given model.id
#' getControls(object=pdxe, model.id="X-6047.16")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return a \code{vector} with control model.id
setGeneric(name = "getControls", def = function(object, model.id) {standardGeneric("getControls")} )


#' @export
setMethod( f=getControls,
           signature="XenoSet",
           definition= function(object,model.id)
           {
             rtx = list()
             for(ed in object@expDesign)
             {
               if(is.element(model.id, ed$treatment))
               { rtx = .appendToList(rtx, ed$control) }
             }
            return(unlist(rtx))
            })




#' Extract treatment for a given model.id
#'
#' @examples
#' data(pdxe)
#' # extract treatment for a given model.id
#' getTreatment(object=pdxe, model.id="X-6047.21")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return a \code{vector} with treatment model.id
setGeneric(name = "getTreatment", def = function(object, model.id) {standardGeneric("getTreatment")} )


#' @export
setMethod( f=getTreatment,
           signature="XenoSet",
           definition= function(object,model.id)
           {
             rtx = list()
             for(ed in object@expDesign)
             {
               if(is.element(model.id, ed$control))
               { rtx = .appendToList(rtx, ed$treatment) }
             }
             return(unlist(rtx))
           })



.getModelIdIndexInExpDesign<- function(object, model.id)
{
  cntrL = c(); tretL = c()
  for(i in 1:length(object@expDesign))
  {
    ed = object@expDesign[[i]]
    if(is.element(model.id, ed$control))
    { cntrL = c(cntrL, i) }

    if(is.element(model.id, ed$treatment))
    { tretL = c(cntrL, i) }
  }
  return(list(control.indx=cntrL, treat.indx=tretL))
}

#' Get experiment type (treatment or control) for a given model.id
#'
#' @examples
#' data(pdxe)
#' # get experiment type for model.id
#' getTreatment(object=pdxe, model.id="X-6047.21")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return returns \code{treatment} or \code{control}
setGeneric(name = "expType", def = function(object, model.id) {standardGeneric("expType")} )

###------- don't look in pdxe@experiment slot -----------
###------- use always pdxe@expDesign slot ---------------
#' @export
setMethod( f=expType,
           signature="XenoSet",
           definition= function(object,model.id)
           {
             ct.indx = .getModelIdIndexInExpDesign(object, model.id)
             if(length(ct.indx$control.indx)>0 & length(ct.indx$treat.indx)==0)
             {return("control")}

             if(length(ct.indx$control.indx)==0 & length(ct.indx$treat.indx)>0)
             {return("treatment")}

             if(length(ct.indx$control.indx)>0 & length(ct.indx$treat.indx)>0)
             {return("control and treatment")}
             return(NA)
           })






###-------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------

plotTreatmentControl <- function(object, type="compact")
{

  ##----------------- GET STAT ---------------------------------------------------------------------
  groupBy = "biobase.id"
  tret = unlist(sapply(object@expDesign, "[[", "treatment"))
  cont = unlist(sapply(object@expDesign, "[[", "control"))
  exp.type = c(rep("treatment",length(tret)), rep("control", length(cont)))
  df = data.frame(model.id = c(tret, cont),
                  exp.type = exp.type, stringsAsFactors = FALSE)

  dfz = merge(df, object@model[, c("model.id", groupBy)], by.x = "model.id", by.y = "model.id")
  if(type=="compact"){dfz = unique(dfz)}
  dfx = as.data.frame.matrix(table(dfz[, c(groupBy,"exp.type")]))
  dfx[,groupBy] = rownames(dfx)

  modelsWithoutCntr =  dfx[dfx$control==0,groupBy]
  ##----------------------------------------------------------------------------------------------

  dfx$sortOrd = dfx$control*dfx$treatment

  dfx = BBmisc::sortByCol(dfx, c("sortOrd","treatment", "control"), asc = c(FALSE, FALSE, FALSE))
  tx = dfx[, c(groupBy, "treatment")]
  colnames(tx) = c(groupBy, "value")
  tx$variable = "treatment"

  cx = dfx[, c(groupBy, "control")]
  colnames(cx) = c(groupBy, "value")
  cx$variable = "control"

  dfp = rbind(tx,cx)
  dfp[, groupBy] = factor(dfp[, groupBy], levels = tx[, groupBy])
  dfp[dfp$variable=="control", "value"] = -(dfp[dfp$variable=="control", "value"])

  plt = ggplot(dfp, aes_string(x=groupBy, y="value", fill="variable")) +
    geom_bar(stat="identity", position="identity")
  plt = plt + scale_fill_manual(values=c('#91cf60','#fc8d59'))
  plt = plt +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))

  plt = plt + theme( axis.text.x = element_blank())
  #pdf("DATA-raw/treatmentContr.pdf", width = 11, height = 5)
  print(plt)
  #dev.off()
}




