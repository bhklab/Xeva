

###----sensitivity/response ----------
#
# make_response_slot <- function(object)
# {
#   allDrg = unique(unlist(sapply(object@experiment, "[[", c("drug", "join.name"))))
#   allMod = unique(unlist(sapply(object@experiment, "[[", "model.id")))
#   modRes = matrix(NA, nrow = length(allDrg), ncol = length(allMod))
#
#   for(I in object@experiment)
#   {
#
#
#   }
#
# }

.getDrugsForABatch <- function(object, batch)
{
  getDrugForMod <- function(x){ object@experiment[[x]]$drug$join.name}
  treatment = unlist(sapply(batch$treatment, getDrugForMod))
  control   = unlist(sapply(batch$control  , getDrugForMod))
  return(list(treatment=treatment, control=control))
}


.getAngleSummary <- function(object)
{
  rdf = data.frame()
  for(batch in expDesignInfo(object))
  {
    drugX = .getDrugsForABatch(object, batch)
    drug = unique(drugX$treatment)
    rdf = rbind(rdf, data.frame(batch.name = batch$batch.name,
                                drug.join.name=drug,
                                angle= batch$angle))
  }

  rdf <- .factor2char(rdf)
  return(rdf)
}

.summarizePerBatchResponse <- function(object, response.measure = "angle", group.by=NULL, summary.stat ="mean")
{

  if(response.measure == "angle")
  {
    df = .getAngleSummary(object)
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



getValueFromModel <- function(object, model.id, values)
{
  mod = object@experiment[[model.id]]
  vz = c()
  for(v in c(values))
  {
    if(v=="drug.join.name")
    {
      vz[v] = mod$drug$join.name
    } else{ vz[v] = mod[[v]] }
  }
  return(vz)
}

.mapAndAttachColumn <- function(object, df, id.name, map.to)
{
  dfMap = mapModelSlotIds(object, id=df[, id.name], id.name=id.name,
                          map.to= map.to, unique = FALSE)

  if(all(df[, id.name] ==dfMap[, id.name]))
  {df[,map.to] = dfMap[,map.to]} else{stop("column do not match")}
  return(df)
}

.summarizePerModelResponse <- function(object, response.measure = "mRECIST_recomputed",
                                       group.by="patient.id", summary.stat="mean")
{
  if(response.measure =="mRECIST_recomputed")
  {
    df = getmRECIST(object)
    valueColName = "mRECIST"
    summary.stat=";"
  }

  if(response.measure =="slop")
  {
    allModNames = names(object@experiment)
    df = lapply(allModNames, function(model.id){
                values=c("model.id","drug.join.name", "slop")
                getValueFromModel(object, model.id, values)} )
    df = .convertListToDataFram(df)
    valueColName = "slop"
  }

  if(group.by!="model.id")
  {
    df = .mapAndAttachColumn(object, df, id.name="model.id", map.to=group.by)
  }

  mat = .castDataFram(df, row.var="drug.join.name", col.var = group.by,
                      value=valueColName, collapse = summary.stat)

  return(mat)
}



#.checkGroupByValue <- function(object, group.by, response.measure)
#{
#  if(group.by=="batch.name" & )
#}

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
#' @export
summarizeResponse <- function(object, response.measure = "mRECIST_recomputed",
                              group.by="patient.id", summary.stat=c("mean", "median", ";"))
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

  if(response.measure =="mRECIST_recomputed")
  {
    mat <- .summarizePerModelResponse(object, response.measure=response.measure,
                               group.by=group.by)
  }

  if(response.measure =="slop")
  {
    mat <- .summarizePerModelResponse(object, response.measure=response.measure,
                               group.by=group.by, summary.stat=summary.stat)
  }

  if(response.measure =="angle")
  {
    mat <- .summarizePerBatchResponse(object, response.measure = "angle",
                               group.by=group.by, summary.stat=summary.stat)

  }

  return(mat)
}






