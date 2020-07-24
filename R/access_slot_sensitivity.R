## get sensitivity for an Xeva object
#'
#' Get sensitivity for an Xeva object
#' @description
#' Given a Xeva object, it will return a \code{data.frame} detailing sensitivity information.
#'
#' @examples
#' data(brca)
#' head(sensitivity(brca, type="batch"))
#' head(sensitivity(brca, type="model"))
#' @param object The \code{Xeva} dataset.
#' @param type Sensitivity type (either model or batch).
#' @param sensitivity.measure Name of the \code{sensitivity.measure}. Default \code{NULL} will return all sensitivity measures.
#'
#' @return A \code{data.frame} with model or batch ID and sensitivity values.
#' @export
sensitivity <- function(object, type=c("model", "batch"), sensitivity.measure=NULL)
{
  type <- match.arg(type)
  senNames <- names(slot(object, "sensitivity"))

  if(is.element(type, senNames)==FALSE)
  {
    msg <- sprintf("sensitivity 'type' can be:\n%s", paste(senNames, collapse = "\n"))
    stop(msg)
  }

  sm <- slot(object, "sensitivity")[[type]]
  if(!is.null(sensitivity.measure))
  {
    sm2take <- intersect(sensitivity.measure, colnames(sm))
    if(length(sm2take)>0)
    {
      if(length(sensitivity.measure)!=length(sm2take))
      {
        notPresent <- setdiff(sensitivity.measure, sm2take)
        msg <- sprintf("some sensitivity.measure are not present in sensitivity slot and will be ignored\n%s\n",
                       paste(notPresent, collapse="\n") )
        warning(msg)
      }
      if(type=="model"){ sm <- sm[,c("model.id", sm2take)] }
      if(type=="batch"){ sm <- sm[,c("batch.name", sm2take)]}

    } else
    {
      msg <- sprintf("sensitivity.measure are not present in sensitivity slot\n%s\n",
                     paste(sensitivity.measure,collapse="\n") )
      stop(msg)
    }
  }
  return(sm)
}



## add new sensitivity in a Xeva object
#'
#' add new sensitivity in a Xeva object
#' @description
#' This will add a add new sensitivity in a Xeva object
#'
#' @examples
#' data(brca)
#' s <- sensitivity(brca, "model")
#' brca <- setSensitivity(brca, "model", "mR", s$mRECIST)
#' @param object The \code{Xeva} dataset
#' @param type sensitivity type (either model or batch)
#' @param name name of new sensitivity column
#' @param value a vector of values. If vector is named, values will be filled by name. Missing values will be NA
#'
#' @return a \code{Xeva} object with updated with updated sensitivity
#' @keywords internal
#' @noRd
setGeneric(name= "setSensitivity", def = function(object, type, name, value)
  {standardGeneric("setSensitivity")} )

#' @keywords internal
#' @noRd
####@export
setMethod( f="setSensitivity",
           signature= signature(object="XevaSet"),
           definition=function(object, type, name, value)
           {
             if(is.element(type, names(slot(object, "sensitivity")))==TRUE )
             {
               if(is.null(names(value)))
               {
                 if(length(value)!= nrow(slot(object, "sensitivity")[[type]]))
                 {
                   msg <- sprintf("vector 'value' must be of same length as other sensitivity")
                   stop(msg)
                 }
                 slot(object, "sensitivity")[[type]][, name] <- NA
                 slot(object, "sensitivity")[[type]][, name] <- value
               } else
               {
                 slot(object, "sensitivity")[[type]][, name] <- NA
                 slot(object, "sensitivity")[[type]][, name] <-
                   value[rownames(slot(object, "sensitivity")[[type]])]
               }
             } else
             {
               msg <- sprintf("sensitivity 'type' can be:\n%s",
                              paste(names(slot(object, "sensitivity")), collapse = "\n"))
               stop(msg)
             }
             return(object)
           } )

