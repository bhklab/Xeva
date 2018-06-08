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
#' Get all batch names
#'
#' Get all batch.name from a Xeva dataset
#' @examples
#' data(pdxe)
#' batchNames(pdxe)
#' @param object The \code{XevaSet} to replace drug info in
#' @return A \code{Vector} with all batch.name
setGeneric(name= "batchNames", def = function(object) {standardGeneric("batchNames")} )

#' @export
setMethod( f="batchNames",
           signature=c(object = "XevaSet"),
           definition=function(object)
           {
             rtx = names(expDesignInfo(object))
             return(rtx)
           } )

#' Given a batch.name get batch
#'
#' Given a batch.name get batch from a Xeva dataset
#' @examples
#' data(pdxe)
#' expDesign(pdxe, batch.name = "X-6047.paclitaxel")
#' @param object The \code{XevaSet}
#' @param object The \code{batch.name}
#' @return A \code{Vector} with all batch.name
setGeneric(name= "expDesign", def = function(object, batch.name) {standardGeneric("expDesign")} )

#' @export
setMethod( f="expDesign",
           signature=c(object = "XevaSet"),
           definition=function(object, batch.name)
           {
             #sapply(object@expDesign, function(x){x$batch.name == batch.name})
             #for(x in object@expDesign)
             #{
             # if(x$batch.name == batch.name)return(x)
             #}
             bt <- expDesignInfo(object)[[batch.name]]
             if(is.null(bt))
             {
               msg = sprintf("batch name %s not present\nuse batchNames(object) to see all batch names", batch.name)
               stop(msg)
             }
             if(!is.null(bt) & bt$batch.name!= batch.name)
             {
               msg = sprintf("Batch slot name are different then batch.name")
               stop(msg)
             }
             return(bt)
           } )

##-------------------------------------------------------------------------------------
##-------------------------------expDesign sanity check -------------------------------
.sanityCheckExpDesign <- function(object, expDesign)
{
  if(length(expDesign$control)==0 & length(expDesign$treatment)==0)
  {stop("Error Treatmetn and Control are missing in expDesign!")}


  if(length(expDesign$treatment)>0)
  {
    tr = lapply(expDesign$treatment, function(x){getExperiment(object, model.id=x)})
    trx= .rbindListOfDataframs(tr)

    drgNamesTr = unique(trx$drug.join.name)
    if(length(drgNamesTr)>1)
    {

      msg = sprintf("For batch.name = %s drug names in treatment is not same. Drug names are:\n%s\n",
                    expDesign$batch.name, paste(drgNamesTr, collapse = "\n"))
      warning(msg)
    }


    patient.idx = mapModelSlotIds(object, id=trx$model.id, id.name="model.id", map.to="patient.id")
    trx = merge(trx, patient.idx, by.x = "model.id", by.y = "model.id")
    patName= unique(trx$patient.id)
    if(length(patName)>1)
    {
      msg = sprintf("For batch.name = %s patient.id in treatment are not same. patient.id are:\n%s\n",
                    expDesign$batch.name, paste(patName, collapse = "\n"))
      warning(msg)
    }
  }

  if(length(expDesign$control)>0)
  {
    cr = lapply(expDesign$control, function(x){getExperiment(object, model.id=x)})
    crx= .rbindListOfDataframs(cr)

    drgNamesCr = unique(crx$drug.join.name)
    if(length(drgNamesCr)>1)
    {
      msg = sprintf("For batch.name = %s drug name in control is not same. Drug names are:\n%s\n",
                    expDesign$batch.name, paste(drgNamesTr, collapse = "\n"))
      warning(msg)
    }

    patient.idx = mapModelSlotIds(object, id=crx$model.id, id.name="model.id", map.to="patient.id")
    trx = merge(crx, patient.idx, by.x = "model.id", by.y = "model.id")
    patName= unique(crx$patient.id)
    if(length(patName)>1)
    {
      msg = sprintf("For batch.name = %s patient.id in control are not same. patient.id are:\n%s\n",
                    expDesign$batch.name, paste(patName, collapse = "\n"))
      warning(msg)
    }
  }

}

##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
#' Get controls for a given model.id
#'
#' Get controls for a given model.id.
#' If no control found it will return NULL
#'
#' @examples
#' data(pdxe)
#' # extract controls for a given model.id
#' getControls(object=pdxe, model.id="X.1655.LE11.biib")
#' # if no control found it will return NULL
#' getControls(object=pdxe, model.id="X.1655.uned")
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


#' Get treatment for a given model.id
#'
#' Get treatment for a given model.id.
#' If no treatment found it will return NULL
#' @examples
#' data(pdxe)
#' # extract treatment model.id for a given model.id
#' getTreatment(object=pdxe, model.id="X.1655.uned")
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


#' Get batch.name for a given model.id
#'
#'
#' Get batch.name for a given model.id.
#' If no batch.name found it will return NULL
#' @examples
#' data(pdxe)
#' # extract batch.name for a given model.id
#' getBatchName(object=pdxe, model.id="X.1655.uned")
#' getBatchName(object=pdxe, model.id="X.010.fiab")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id} for which batch name required
#' @return a \code{vector} with all batch names
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
#' experimentType(object=pdxe, model.id="X.1655.uned")
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#' @return returns \code{treatment} or \code{control}
setGeneric(name = "experimentType", def = function(object, model.id) {standardGeneric("experimentType")} )

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



##-------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------
## get ExpDesign Slot DF
#' Given a model.id it will return a data.fram of experiemt design
#' with columns as "treatment", "control", "batch.name"
#'
#' @examples
#' data(pdxe)
#' # This will give a data.fram with columns as "treatment", "control", "batch.name"
#' getExpDesignDF(object=pdxe, model.id="X.1655.LE11.biib")
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


##------------------------------------------------------------------------------
.smoothCurveOld <- function(dt, x, y)
{
  #lo <- loess(y~x)
  naIndx <- which(is.na(dt[,y]))
  notNaIndx <- which(!is.na(dt[,y]))

  naDF <- dt[naIndx, ]
  notNaDF <- dt[notNaIndx, ]

  lo <- loess(dt[,y]~dt[,x])
  #lo <- loess(notNaDF[,y]~notNaDF[,x])

  for(i in 1:nrow(dt))
  {
    if(is.na(dt[i, y]))
    {
      dt[i, y] <- predict(lo, newdata = dt[i, x])
    }
  }
  return(dt)
}

.smoothCurve <- function(trDF, tsDF, x, y)
{
  ##lo <- loess(y~x)
  #lo <- loess(trDF[,y]~trDF[,x])
  #tsDF[,y] <- predict(lo, newdata = tsDF[,x])
  #

  lmdf <- data.frame(X=trDF[,x], Y=trDF[,y])
  lmod <- lm(Y~X, data = lmdf)
  prdf <- data.frame(X=tsDF[,x])
  tsDF[,y] <- predict(lmod, prdf, se.fit = F)
  return(tsDF)
}



.smoothModel <- function(dw, timeVec, var = "volume")
{
  dw$impute.value <- "NO"
  texra <- timeVec[timeVec<max(dw$time)]
  texra <- setdiff(texra, dw$time)
  texra <- texra[is.numeric(texra)]

  if(length(texra)>0)
  {
    nwRws <- data.frame(matrix(data = NA, nrow = length(texra), ncol = ncol(dw)),
                        stringsAsFactors = F)
    colnames(nwRws) <- colnames(dw)
    nwRws$time <- texra
    nwRws$model.id <- dw$model.id[1]
    nwRws$drug.join.name <- dw$drug.join.name[1]
    nwRws$impute.value <- "YES"

    nwRws <- .smoothCurve(dw, nwRws, x="time", y=var)

    dt <- rbind(dw, nwRws)
    dt <- BBmisc::sortByCol(dt, c("time"))

    #ds <- .smoothCurve(dt, x="time", y=var)
    dt$volume.normal <- .normalizeVolume(dt$volume)
    return(dt)
  }

  return(dw)
}

##------------------------------------------------------------------------------
## collapse time-vol data based on expDesign
.collapseRplicate <- function(inLst, var = "volume", impute.value=TRUE)
{
  if(is.null(names(inLst))){names(inLst) <- sapply(inLst, function(x){ x$model.id[1]})}


  if(impute.value==TRUE)
  {
    inLst2 <- list()
    timeVec <- sort(unique(unlist(lapply(inLst, "[[", "time" ))))
    for(mid in names(inLst))
    {
      dw <- inLst[[mid]]
      inLst2[[mid]] <- .smoothModel(dw, timeVec=timeVec, var =var)
    }
    inLst <- inLst2
  }

  #dw<- reshape2::dcast(dq, as.formula(paste("time", "model.id", sep="~")),
  #                     value.var = var)

  timeAll <- sort(unique(unlist(lapply(inLst, "[[", "time" ))))
  rd = data.frame()
  for(t in timeAll)
  {
    vx <- unlist(sapply(inLst, function(x){ x[x$time==t, var]}))
    vx <- vx[!is.na(vx)]
    vz <- as.list(Rmisc::STDERR(vx))
    rd <- rbind(rd, data.frame(time=t, mean=vz$mean, upper= vz$upper, lower= vz$lower))
  }
  return(rd)
}


## get time vs volume data with standard error
#'
#' Get time vs volume data with standard error
#'
#'
#' Given a batch (treatment and control model ids)
#' it will return a data.fram with time vs volume (or any other variable)
#' with standard error calculated. Note that this function do not check if
#' model.id in given batch belongs to same patient
#' Note: Write a function to check integrity of a batch
#'
#' @examples
#' data(pdxe)
#' # creat a experiment desing
#' ExpDesign = list(batch.name="myBatch", treatment=c("X.010.BG98", "X.010.BG98"), control=c("X.010.uned"))
#' df = getTimeVarData(object=pdxe, ExpDesign, var = "volume")
#'
#' @param object The \code{Xeva} dataset
#' @param ExpDesign A list with batch.name, treatment and control
#' @param var Name of the variable, default \code{volume}
#' @param drug.name \code{FALSE}. If \code{TRUE} will return drug name also
#' @return a \code{data.fram} with treatment, control and batch.name
setGeneric(name = "getTimeVarData",
           def = function(object, ExpDesign=NULL, batchName=NULL, var = "volume",
                          treatment.only=FALSE, drug.name=FALSE,
                          vol.normal=FALSE, impute.value=TRUE)
  {standardGeneric("getTimeVarData")} )

#' @export
setMethod( f=getTimeVarData,
           signature=c(object="XevaSet"),
           definition= function(object, ExpDesign=NULL, batchName=NULL, var = "volume",
                                treatment.only=FALSE, drug.name=FALSE,
                                vol.normal=FALSE, impute.value=TRUE)
           {

             if(is.null(batchName) & is.null(ExpDesign))
             {
               stop("please provide 'batchName' or 'ExpDesign'")
             }

             if(!is.null(batchName) & is.null(ExpDesign))
             {
               ExpDesign <- expDesign(object, batchName)
             }

             if(is.null(ExpDesign$treatment) & is.null(ExpDesign$control))
             {stop("treatment and control both NULL")}
             if(length(ExpDesign$treatment)==0 & length(ExpDesign$control)==0)
             {stop("treatment and control both NULL")}

             trDF = data.frame()
             cnDF = data.frame()
             if(!is.null(ExpDesign$treatment) & length(ExpDesign$treatment)>0 )
             {
               trLs <- .getExperimentMultipalIDs(object, ExpDesign$treatment,
                                                treatment.only=treatment.only
                                                ,vol.normal=vol.normal
                                                )
               trDF <- .collapseRplicate(trLs, var = var, impute.value=impute.value)
               trDF$exp.type = "treatment"
               trDF$batch.name = ExpDesign$batch.name
               if(drug.name==TRUE)
               {
                 drugAll = sort(unique(unlist(lapply(trLs, "[[", "drug.join.name" ))))
                 if(length(drugAll)>1)
                 {
                   msg <- sprintf("multipal drugs for batch, will colleps by ;")
                   warning(msg)
                 }
                 trDF$drug.name <- paste(drugAll, collapse = ";")
               }
             }

             if(!is.null(ExpDesign$control) & length(ExpDesign$control)>0)
              {
               cnLs = .getExperimentMultipalIDs(object, ExpDesign$control,
                                                treatment.only=treatment.only
                                                ,vol.normal=vol.normal
                                                )
               cnDF = .collapseRplicate(cnLs, var = var, impute.value=impute.value)
               cnDF$exp.type = "control"
               cnDF$batch.name = ExpDesign$batch.name

               if(drug.name==TRUE)
               {
                 drugAll = sort(unique(unlist(lapply(cnLs, "[[", "drug.join.name" ))))
                 if(length(drugAll)>1)
                 {
                   msg <- sprintf("multipal drugs for batch, will colleps by ;")
                   warning(msg)
                 }
                 cnDF$drug.name <- paste(drugAll, collapse = ";")
               }

             }

             # if(vol.normal==TRUE)
             # {
             #   if(nrow(trDF)>0)
             #   {
             #     trDF$mean.raw <- trDF$mean
             #     trDF$mean <- .normalizeVolume(trDF$mean.raw)
             #   }
             #
             #   if(nrow(cnDF)>0)
             #   {
             #     cnDF$mean.raw <- cnDF$mean
             #     cnDF$mean <- .normalizeVolume(cnDF$mean.raw)
             #   }
             # }

             rdf <- rbind(trDF, cnDF)
             return(rdf)
           })

.getExperimentMultipalIDs <- function(object, mids, treatment.only=TRUE, vol.normal=FALSE
                                      )
{
  rtx = list()
  for(i in mids)
  {
    miD = getExperiment(object, model.id= i,treatment.only=treatment.only)
    if(vol.normal==TRUE)
    {
      miD$volume.raw <- miD$volume
      miD$volume <- miD$volume.normal
    }
    rtx = .appendToList(rtx, miD)
  }
  return(rtx)
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

