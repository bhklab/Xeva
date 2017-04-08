
.getDrugsForABatch <- function(object, batch)
{
  getDrugForMod <- function(x){ object@experiment[[x]]$drug$join.name}
  treatment = unlist(sapply(batch$treatment, getDrugForMod))
  control   = unlist(sapply(batch$control  , getDrugForMod))
  return(list(treatment=treatment, control=control))
}


.getAngleSummary <- function(object, tumor.type)
{
  rdf = data.frame()
  for(batch in expDesignInfo(object))
  {
    if(!is.null(tumor.type))
    {
      tt <- modelInfo(object)[modelInfo(object)$model.id %in% c(batch$treatment,batch$control),"tumor.type"]

      if(is.element(tumor.type, tt)==FALSE)
      { next }

      if(length(unique(tt))>1)
      {
        msg <- sprintf("batch %s contains models with different tumor type", batch$batch.name)
        warning(msg)
      }
    }

    drugX = .getDrugsForABatch(object, batch)
    drug = unique(drugX$treatment)
    rdf = rbind(rdf, data.frame(batch.name = batch$batch.name,
                                drug.join.name=drug))
  }
  rdf <- .factor2char(rdf)
  rdf$angle <- slot(object, "sensitivity")[["batch"]][rdf$batch.name, ]
  return(rdf)
}

.summarizePerBatchResponse <- function(object, response.measure = "angle", group.by=NULL,
                                       summary.stat, tumor.type)
{

  if(response.measure == "angle")
  {
    df = .getAngleSummary(object, tumor.type=tumor.type)
  }

  if(!is.null(group.by) & group.by != "batch.name")
  {
    mapId = mapModelSlotIds(object, id=df$batch.name, id.name="batch.name",
                            map.to=group.by, unique = TRUE)

    for(I in unique(mapId$batch.name))
    {
      di = unique(mapId[mapId$batch.name==I, group.by])
      if(length(di)>1)
      {
        msg1 = sprintf("batch.name mapped to multipal %s therefore all such %s will contain same information", group.by, group.by)
        msg2 = sprintf("\nbatch.name = %s %s = %s\n", I, group.by, paste(di, collapse =","))
        warning(msg1, msg2)
      }
    }

    df = merge(df, mapId, by.x = "batch.name", by.y="batch.name")

    for(I in unique(df[, group.by]))
    {
      di = df[df[, group.by]==I, ]
      if(nrow(di)>1)
      {
        drg <- di$drug.join.name
        if(length(drg) > length(unique(drg)))
        {
          msg1 <- sprintf("'%s' mapped to multipal batch.name, values will be collapsed using '%s'\n", group.by, summary.stat)
          df2Pr<- unique(di[, c("batch.name", group.by)]); rownames(df2Pr)<- NULL
          msg2 <- paste(capture.output(print(df2Pr)), collapse = "\n")
          warning(msg1, msg2)
        }
      }
    }
  }

  mat = .castDataFram(df, row.var="drug.join.name", col.var = group.by,
                      value=response.measure, collapse = summary.stat)
  return(mat)
}



.getValueFrom1Model <- function(object, model.id, values)
{
  mod = slot(object, "experiment")[[model.id]]
  vz = list() #c()
  for(v in c(values))
  {
    if(v=="drug.join.name")
    {
      vl = mod$drug$join.name
    } else
    {
      vl = mod[[v]]
    }

    if(is.null(vl)) {vl = NA}
    vz[[v]] = vl
  }
  return(vz)
}

.getValueFromAllModel <- function(object, values, tumor.type)
{
  mInfo <- modelInfo(object)
  allModNames = mInfo[, "model.id"]
  if(!is.null(tumor.type))
  {
    allModNames = mInfo[mInfo$tumor.type==tumor.type, "model.id"]
    if(length(allModNames)<1)
    {
      tx <- paste(unique(mInfo$tumor.type), collapse = "\n")
      msg <- sprintf("No model found for tumor.type %s\ntumor.type are %s", tumor.type, tx)
      stop(msg)
    }
  }

  #allModNames = names(slot(object, "experiment"))
  dfL = lapply(allModNames, function(model.id)
              { .getValueFrom1Model(object, model.id, values)})
  df = .convertListToDataFram(dfL)
  return(df)
}

.mapAndAttachColumn <- function(object, df, id.name, map.to)
{
  dfMap = mapModelSlotIds(object, id=df[, id.name], id.name=id.name,
                          map.to= map.to, unique = FALSE)

  if(all(df[, id.name] ==dfMap[, id.name]))
  {df[,map.to] = dfMap[,map.to]} else{stop("column do not match")}
  return(df)
}

.summarizePerModelResponse <- function(object, response.measure, group.by, summary.stat, tumor.type)
{
  if(is.element(response.measure,  colnames(slot(object, "sensitivity")[["model"]]) )==FALSE)
  {
    msg <- sprintf("'%s' is not present in sensitivity slot\n", response.measure)
    stop(msg)
  }

  dfVal <- slot(object, "sensitivity")[["model"]] [,c("model.id", response.measure)]

  df <- modelInfo(object)
  df[, response.measure] <- dfVal[df$model.id, response.measure]

  if(!is.null(tumor.type))
  {
    df <- df[df$tumor.type==tumor.type,]
  }

  if(is.element(group.by, colnames(df))==FALSE)
  {
    msg <- sprintf("'group.by' %s not present in model\n", group.by)
    stop(msg)
  }

  mat = .castDataFram(df, row.var="drug", col.var = group.by,
                      value=response.measure, collapse = summary.stat)

  ###------------------------------------------------
  #values <- c("model.id","drug.join.name")
  #valueColName <- NULL

  #if(response.measure =="mRECIST_recomputed")
  #{
  #  valueColName =="mRECIST"
  #  summary.stat=";"
  #}
  #
  # if(response.measure =="mRECIST_recomputed" | response.measure =="mRECIST")
  # {
  #   valueColName <- "mRECIST"
  #   summary.stat=";"
  # }
  #
  # if(response.measure =="mRECIST_published")
  # {
  #   valueColName = "mRECIST_published"
  #   summary.stat=";"
  # }
  #
  # if(response.measure =="slope")
  # {
  #   valueColName <- "slope"
  # }
  #
  # ##------------------------------------------------------
  # if(is.null(valueColName))
  # {
  #   varPresTest <- checkExperimentSlotVariable(object, response.measure)
  #   if(varPresTest==TRUE)
  #   { valueColName <- response.measure }
  # }
  #
  # values <- c("model.id","drug.join.name", valueColName)
  # df <- .getValueFromAllModel(object, values, tumor.type)
  #
  # if(all(is.na(df[, valueColName]))==TRUE)
  # {
  #   msg = sprintf("all values for %s are NA", valueColName)
  #   warning(msg)
  # }
  #
  # if(group.by!="model.id")
  # {
  #   df = .mapAndAttachColumn(object, df, id.name="model.id", map.to=group.by)
  # }

  #mat = .castDataFram(df, row.var="drug.join.name", col.var = group.by,
  #                    value=valueColName, collapse = summary.stat)

  return(mat)
}



#####================= summarizeResponse ==================
#' Summarize Response of PDXs
#'
#' Summarize Response of PDXs.
#'
#' @param object The \code{XevaSet}
#' @param response.measure \code{character} . Which response measure to use? Use the responseMeasures function to find out what measures are available for each Xeva set.
#' @param group.by default \code{patient.id}. How the models should be grouped togather. See details
#' @param summary.stat which summary method to use if multipal ids were found
#' @return a \code{matrix} with rows as drug names, coulmn as \code{group.by} and each cell contains \code{response.measure} for the pair.
#'
#' @details
#' There can be two types of response measure
#' \itemize{
#' \item{per model response : One response value for each Model, e.g. mRECIST_recomputed for each model}
#' \item{per batch response : One response value for each Batch, e.g. angle between treatment and control groups}
#' }
#' In case of \code{per model response} output columns will be \code{model.id} (or group.by).
#' For \code{per batch response} \code{group.by} value can be \code{"batch.name"} .
#'
#' @examples
#' data(pdxe)
#' pdxe_mR <- summarizeResponse(pdxe, response.measure = "mRECIST_recomputed", group.by="patient.id")
#' #to get only lung PDXE
#' pdxe_mR <- summarizeResponse(pdxe, response.measure = "mRECIST_recomputed",
#'                              group.by="patient.id", tumor.type="NSCLC")
#' @export
summarizeResponse <- function(object, response.measure = "mRECIST_recomputed",
                              group.by="patient.id", summary.stat=c(";", "mean", "median"),
                              tumor.type=NULL)
{

  summary.stat = c(summary.stat)[1]

  if(group.by=="batch.name")
  {
    if(is.element(response.measure, c("angle"))==FALSE)
    {
      msg1 = sprintf("group.by== 'batch.name' is only allowed in per batch response" )
      stop(msg1)
    }
  }

  if(response.measure =="angle")
  {
    mat <- .summarizePerBatchResponse(object, response.measure = "angle",
                                      group.by=group.by, summary.stat=summary.stat,
                                      tumor.type=tumor.type)

  } else
  {
    mat <- .summarizePerModelResponse(object, response.measure=response.measure,
                                      group.by=group.by, summary.stat=summary.stat,
                                      tumor.type=tumor.type)
  }


  return(mat)
}






