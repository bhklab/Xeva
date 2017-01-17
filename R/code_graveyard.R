
##=======================================================================================
##=======================================================================================
.checkDrugNameSame <- function(inLst)
{
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



#' ## get time vs volume data with standard error
#' #'
#' #' Given a treatment and control model ids
#' #' it will return a data.fram with time vs volume (or any other variable)
#' #'
#' #' @examples
#' #' data(pdxe)
#' #' # creat a experiment desing
#' #' ExpDesign = list(batch.name="myBatch", treatment=c("X.010.BG98"), control=c("X.010.uned"))
#' #' df = getTimeVarData(object=pdxe, ExpDesign, var = "volume", collapse=TRUE)
#' #' ## if collapse=FALSE it will not calculate standard error
#' #' df2= getTimeVarData(object=pdxe, ExpDesign, var = "volume", collapse=FALSE)
#' #'
#' #' @param object The \code{Xeva} dataset
#' #' @param ExpDesign A list with batch.name, treatment and control
#' #' @param var Name of the variable, default \code{volume}
#' #' @param collapse Default \code{TRUE}. It will summerize all models which belongs to same treatment or control group
#' #' @return a \code{data.fram} with treatment, control and batch.name
#' setGeneric(name = "Old_getTimeVarData", def = function(object, ExpDesign, var, collapse) {standardGeneric("Old_getTimeVarData")} )
#'
#' ########## @export
#' setMethod( f=Old_getTimeVarData,
#'            signature=c(object="XevaSet", ExpDesign="list"),
#'            definition= function(object, ExpDesign, var = "volume", collapse=TRUE)
#'            {
#'
#'              rtxTret = list()
#'              for(i in ExpDesign$treatment)
#'              {
#'                expData = getExperiment(object, model.id= i)
#'                expData = expData[, c("model.id", "drug.join.name", "time", var)]
#'                patient.idx = mapModelSlotIds(object, id=expData$model.id, id.name="model.id", map.to="patient.id")
#'                expData$patient.id = patient.idx[1, "patient.id"]
#'                rtxTret = .appendToList(rtxTret, expData)
#'              }
#'              ###-------------------------------------------------------------------------------------
#'              rtxCont = list()
#'              for(j in ExpDesign$control)
#'              {
#'                conData = getExperiment(object, model.id= i)
#'                conData = conData[, c("model.id", "drug.join.name", "time", var)]
#'                patient.idx = mapModelSlotIds(object, id=conData$model.id, id.name="model.id", map.to="patient.id")
#'                conData$patient.id = patient.idx[1, "patient.id"]
#'
#'                rtxCont = .appendToList(rtxCont, conData)
#'              }
#'
#'              if(collapse==TRUE)
#'              {
#'                trD = .checkDrugNameSame(rtxTret)
#'                tpi = .checkPatientIDSame(rtxTret)
#'
#'                rtxTretX =.collapseRplicate(rtxTret, var)
#'                rtxTretX$drug.join.name = trD
#'                rtxTretX$patient.id = tpi
#'
#'                ##----------for control ------------------------
#'                cnD = .checkDrugNameSame(rtxCont)
#'                cpi = .checkPatientIDSame(rtxCont)
#'
#'                rtxContX =.collapseRplicate(rtxCont, var)
#'                rtxContX$drug.join.name = cnD
#'                rtxContX$patient.id = cpi
#'
#'              } else
#'              {
#'                rtxTretX = do.call(rbind, rtxTret)
#'                rtxContX = do.call(rbind, rtxCont)
#'              }
#'
#'              rtxTretX$exp.type="treatment"
#'              rtxContX$exp.type="control"
#'              rtX = rbind(rtxTretX, rtxContX)
#'              rtX$batch.name = ExpDesign$batch.name
#'              return(rtX)
#'            })
#'
#'
#'
#'
#' ##-------------------------------------------------------------------------------------
#' ##-------------------------------------------------------------------------------------







computAngelFor1ExpDesign_old <- function(object, expDegI, var="volume", treatment.only=TRUE, plot=FALSE)
{
  DFx = getTimeVarData(object, expDegI, var=var, treatment.only=treatment.only)


  dfC = DFx[DFx$exp.type=="control",]
  dfT = DFx[DFx$exp.type=="treatment",]

  if(nrow(dfC)==0)
  {
    warning("Control have no data!")
    return(NA)
  }

  if(nrow(dfT)==0)
  {
    warning("treatment have no data!")
    return(NA)
  }
  fitC = .computSlopFun(dfC$time, dfC$mean)
  fitT = .computSlopFun(dfT$time, dfT$mean)

  angDiff = fitC$angel - fitT$angel

  if(plot==TRUE)
  { .plotAngelAndFit(dfT, dfC, fitT, fitC) }

  return(angDiff)
}
