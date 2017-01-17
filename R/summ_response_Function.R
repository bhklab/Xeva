

###----sensitivity/response ----------

make_response_slot <- function(object)
{
  allDrg = unique(unlist(sapply(object@experiment, "[[", c("drug", "join.name"))))
  allMod = unique(unlist(sapply(object@experiment, "[[", "model.id")))
  modRes = matrix(NA, nrow = length(allDrg), ncol = length(allMod))

  for(I in object@experiment)
  {


  }

}

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
  return(rdf)
}

.summarizePerBatchResponse <- function(object, response.measure = "angle")
{
  if(response.measure == "angle")
  {
    df = .getAngleSummary(object)
  }

  mat = .castDataFram(df, row.var="drug.join.name", col.var = "batch.name",
                      value=response.measure, collapse = ";")
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
    #if(any(is.na(df$slop)))
    #{
    #  df = df[!is.na(df$slop), ]
    #  warning("removing slope with NA")
    #}
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




#####================= summarizeResponse ==================
#' Summarize Response of PDXs
#'
#' Summarize Response of PDXs.
#'
#' @param object The \code{XevaSet}
#' @param response.measure \code{character} . Which response measure to use? Use the responseMeasures function to find out what measures are available for each Xeva set.
#' @param group.by default \code{patient.id}. How the models should be grouped togather.
#' @param summary.stat which summary method to use if more then multipal model.id response is found
#' @return a \code{matrix} with rows as drug names, coulmn as \code{group.by} and each cell contains \code{response.measure} for the pair.
#'
#' @details
#' There can be two types of response measure
#' \itemize{
#' \item{per model response : One response value for each Model, e.g. mRECIST_recomputed for each model}
#' \item{per batch response : One response value for each Batch, e.g. angle between treatment and control groups}
#' }
#' In case of \code{per model response} output columns will be \code{model.id} (or group.by).
#' For \code{per batch response} output columns will be \code{batch.name}.
#'
#' @examples
#' data(pdxe)
#' pdxe_mR <- summarizeResponse(pdxe, response.measure = "mRECIST_recomputed", group.by="patient.id")
#' @export
summarizeResponse <- function(object, response.measure = "mRECIST_recomputed",
                              group.by="patient.id", summary.stat=c("mean", "median", ";"))
{

  summary.stat = c(summary.stat)[1]

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
    mat <- .summarizePerBatchResponse(object, response.measure = "angle")
  }

  return(mat)
}






