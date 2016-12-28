##----- get expDesignInfo -------------
#' expDesignInfo Generic
#' Generic for expDesignInfo method
#'
#' @examples
#' data(pdxe)
#' expDesignInfo(pdxe)
#' @param object The \code{XevaSet} to retrieve drug info from
#' @return a \code{list} with the all experiment designs
setGeneric(name = "expDesignInfo", def = function(object) {standardGeneric("expDesignInfo")} )

#' @export
setMethod( f=expDesignInfo, signature="XevaSet",
           definition=function(object)
           { object@expDesign} )


#' expDesignInfo<- Generic
#' Generic for expDesignInfo replace method
#' @examples
#' data(pdxe)
#' expDesignInfo(pdxe) <- expDesignInfo(pdxe)
#' @param object The \code{XevaSet} to replace drug info in
#' @param value A \code{list} with the experiment designs
#' @return Updated \code{XevaSet}
setGeneric(name= "expDesignInfo<-", def = function(object, value) {standardGeneric("expDesignInfo<-")} )

#' @export
setMethod( f="expDesignInfo<-",
           signature=c(object = "XevaSet", value="list"),
           definition=function(object, value)
           {
             object@expDesign = value
             return(object)
           } )



##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
## get ExpDesign Slot DF
#' Given a model.id it will return a data.fram of experiemt design
#'
#' @examples
#' data(pdxe)
#' # extract controls for a given model.id
#' getExpDesignDF(object=pdxe, model.id="X.1655.LE11.biib")
#' getExpDesignDF(object=pdxe, model.id="X-6047.16")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return a \code{data.fram} with treatment, control and batch.name
setGeneric(name = "getExpDesignDF", def = function(object, model.id) {standardGeneric("getExpDesignDF")} )

#' @export
setMethod( f=getExpDesignDF,
           signature="XevaSet",
           definition= function(object,model.id)
           {
             rtx = data.frame("treatment" = character(), "control" = character(),
                              "batch.name" = character(), stringsAsFactors=FALSE)

             for(ed in expDesignInfo(object))
             {
               if(is.element(model.id, ed$treatment) | is.element(model.id, ed$control))
               {
                 for(i in ed$treatment)
                 {
                   for(j in ed$control)
                   {
                     w = data.frame(treatment = i, control=j, batch.name=ed$batch.name,stringsAsFactors=FALSE)
                     rtx = rbind(rtx,w)
                   }
                 }
               }
             }
             return(rtx)
           })


##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
## collapse time-vol data based on expDesign

.collapseRplicate <- function(inLst, var = "volume")
{
  timeAll = sort(unique(unlist(lapply(inLst, "[[", "time" ))))
  rd = data.frame()
  for(t in timeAll)
  {
    vx = unlist(sapply(inLst, function(x){ x[x$time==t, var]}))
    vz = as.list(Rmisc::STDERR(vx))
    rd = rbind(rd, data.frame(time=t, mean=vz$mean, upper= vz$upper, lower= vz$lower))
  }
  return(rd)
}


.checkDrugNameSame <- function(inLst)
{
  ##----- add drug.join.name -------------
  Tdfx = do.call(rbind, inLst)
  drgNames = c(unique(Tdfx$drug.join.name))
  if(length(drgNames)>1)
  {
    txt = sprintf("You are collapsing experiment with different drug. Drug used are:\n%s\n",
                  paste(drgNames, collapse = "\n"))
    warning(txt)
    drgNamesX = paste(drgNames, collapse = " OR ")
  } else{drgNamesX =drgNames[1]}

  return(drgNamesX)
}

.checkPatientIDSame <- function(inLst)
{
  Tdfx = do.call(rbind, inLst)
  pid = c(unique(Tdfx$patient.id))
  if(length(pid)>1)
  {
    txt = sprintf("You are collapsing togather experiment from different patients. Patient ids are:\n%s\n",
                  paste(pid, collapse = "\n"))
    warning(txt)
    pidX = paste(pid, collapse = " OR ")
  } else{pidX = pid[1]}

  return(pidX)
}

##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
## get ExpDesign Slot DF
#' Given a model.id it will return a data.fram of experiemt design
#'
#' @examples
#' data(pdxe)
#' # extract controls for a given model.id
#' ExpDesign = list(batch.name="myBatch", treatment=c("X.010.BG98"), control=c("X.010.uned"))
#' getTimeVarData(object=pdxe, ExpDesign, var = "volume", collapse=TRUE)
#' @param object The \code{Xeva} dataset
#' @param ExpDesign A list with batch.name, treatment and control
#' @return a \code{data.fram} with treatment, control and batch.name
setGeneric(name = "getTimeVarData", def = function(object, ExpDesign, var, collapse) {standardGeneric("getTimeVarData")} )

#' @export
setMethod( f=getTimeVarData,
           signature=c(object="XevaSet", ExpDesign="list"),
           definition= function(object, ExpDesign, var = "volume", collapse=TRUE)
           {
             rtxTret = list()
             for(i in ExpDesign$treatment)
             {
               expID = mapModelSlotIds(object, id=i, id.name="model.id", map.to="experiment.id")
               expData = getExperiment(object, experiment.id= expID[1, "experiment.id"])
               expData = expData[, c("model.id", "experiment.id", "drug.join.name", "time", var)]
               patient.idx = mapModelSlotIds(object, id=expData$model.id, id.name="model.id", map.to="patient.id")
               expData$patient.id = patient.idx[1, "patient.id"]
               rtxTret = .appendToList(rtxTret, expData)
             }
             ###-------------------------------------------------------------------------------------
             rtxCont = list()
             for(j in ExpDesign$control)
             {
               conID = mapModelSlotIds(object, id=j, id.name="model.id", map.to="experiment.id")
               conData = getExperiment(object, experiment.id= conID[1, "experiment.id"])
               conData = conData[, c("model.id", "experiment.id", "drug.join.name", "time", var)]
               patient.idx = mapModelSlotIds(object, id=conData$model.id, id.name="model.id", map.to="patient.id")
               conData$patient.id = patient.idx[1, "patient.id"]

               rtxCont = .appendToList(rtxCont, conData)
             }

             if(collapse==TRUE)
             {
               trD = .checkDrugNameSame(rtxTret)
               tpi = .checkPatientIDSame(rtxTret)

               rtxTretX =.collapseRplicate(rtxTret, var)
               rtxTretX$drug.join.name = trD
               rtxTretX$patient.id = tpi

               ##----------for control ------------------------
               cnD = .checkDrugNameSame(rtxCont)
               cpi = .checkPatientIDSame(rtxCont)

               rtxContX =.collapseRplicate(rtxCont, var)
               rtxContX$drug.join.name = cnD
               rtxContX$patient.id = cpi

             } else
             {
               rtxTretX = do.call(rbind, rtxTret)
               rtxContX = do.call(rbind, rtxCont)
             }

             rtxTretX$exp.type="treatment"
             rtxContX$exp.type="control"
             rtX = rbind(rtxTretX, rtxContX)
             rtX$batch.name = ExpDesign$batch.name
             return(rtX)
           })



##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
#' Extract controls for a given model.id
#'
#' @examples
#' data(pdxe)
#' # extract controls for a given model.id
#' getControls(object=pdxe, model.id="X.1655.LE11.biib")
#' getControls(object=pdxe, model.id="X-6047.16")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return a \code{vector} with control model.id
setGeneric(name = "getControls", def = function(object, model.id) {standardGeneric("getControls")} )

#' @export
setMethod( f=getControls,
           signature="XevaSet",
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
#' getTreatment(object=pdxe, model.id="X.1655.uned")
#' getTreatment(object=pdxe, model.id="X-6047.21")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return a \code{vector} with treatment model.id
setGeneric(name = "getTreatment", def = function(object, model.id) {standardGeneric("getTreatment")} )

#' @export
setMethod( f=getTreatment,
           signature="XevaSet",
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


#' Extract batch.name for a given model.id
#'
#' @examples
#' data(pdxe)
#' # extract treatment for a given model.id
#' getBatchName(object=pdxe, model.id="X.1655.uned")
#' getBatchName(object=pdxe, model.id="X.010.fiab")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return a \code{vector} with treatment model.id
setGeneric(name = "getBatchName", def = function(object, model.id){standardGeneric("getBatchName")} )

#' @export
setMethod( f=getBatchName,
           signature="XevaSet",
           definition= function(object,model.id)
           {
             rtx = list()
             for(ed in object@expDesign)
             {
               if(is.element(model.id, ed$treatment) | is.element(model.id, ed$control))
               { rtx = .appendToList(rtx, ed$batch.name) }
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
#' experimentType(object=pdxe, model.id="X.1655.LE11.biib")
#' getTreatment(object=pdxe, model.id="X-6047.21")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return returns \code{treatment} or \code{control}
setGeneric(name = "experimentType", def = function(object, model.id) {standardGeneric("experimentType")} )

###------- don't look in pdxe@experiment slot -----------
###------- use always pdxe@expDesign slot ---------------
#' @export
setMethod( f=experimentType,
           signature="XevaSet",
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
#' @export
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

  plt = ggplot2::ggplot(dfp, ggplot2::aes_string(x=groupBy, y="value", fill="variable")) +
    ggplot2::geom_bar(stat="identity", position="identity")
  plt = plt + ggplot2::scale_fill_manual(values=c('#91cf60','#fc8d59'))
  plt = .ggplotEmptyTheme(plt)# +  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = element_blank(),
               #panel.background = element_blank(), axis.line = ggplot2::element_line(colour = "black"))

  plt = plt + ggplot2::theme( axis.text.x = ggplot2::element_blank())
  #pdf("DATA-raw/treatmentContr.pdf", width = 11, height = 5)
  print(plt)
  #dev.off()
}

#' @export
experimentDesignSummary <- function(object)
{
  cat(sprintf("number of experiment designs = %d\n", length(object@expDesign)))

  tretAll = sapply(object@expDesign, "[[", "treatment")
  contAll = sapply(object@expDesign, "[[", "control")

  tret = unlist(tretAll)
  cont = unlist(contAll)
  exp.type = c(rep("treatment",length(tret)), rep("control", length(cont)))
  df = data.frame(model.id = c(tret, cont),
                  exp.type = exp.type, stringsAsFactors = FALSE)
  df = unique(df)
  cat(sprintf("number of experiments (in experiment designs) = %d\n", dim(df)[1]))
  cat(sprintf("number of experiments decleared as control= %d\n", dim(df[df$exp.type=="control",])[1]))
  cat(sprintf("number of experiments decleared as treatment= %d\n", dim(df[df$exp.type=="treatment",])[1]))

  tlx = sapply(tretAll, length)
  cat(sprintf("experiment designs without treatment= %d\n", length(tlx[tlx==0])))

  ctx = sapply(contAll, length)
  cat(sprintf("experiment designs without control = %d\n", length(ctx[ctx==0])))

  tretWithoutCont =c()
  contWithoutTret =c()
  for(I in object@expDesign)
  {
    if(length(I$control)==0)  { tretWithoutCont = c(tretWithoutCont, I$treatment) }
    if(length(I$treatment)==0){ contWithoutTret = c(contWithoutTret, I$control) }
  }
  tretWithoutCont =unique(tretWithoutCont)
  contWithoutTret =unique(contWithoutTret)

  cat(sprintf("number of experiments without control = %d\n", length(tretWithoutCont)))
  cat(sprintf("number of experiments without treatment = %d\n", length(contWithoutTret)))

  mdfTr = accessModel(object, tretWithoutCont)
  cat(sprintf("Patients without control = %d\n", length(unique(mdfTr$patient.id))))

  mdfCn = accessModel(object, contWithoutTret)
  cat(sprintf("Patients without treatment = %d\n", length(unique(mdfCn$patient.id))))


}






