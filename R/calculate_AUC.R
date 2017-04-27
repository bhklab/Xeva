##copied from 'pracma' package trapz function
.trapz <- function (x, y)
{
  if (missing(y)) {
    if (length(x) == 0)
      return(0)
    y <- x
    x <- seq(along = x)
  }
  if (length(x) == 0 && length(y) == 0)
    return(0)
  if (!(is.numeric(x) || is.complex(x)) || !(is.numeric(y) || is.complex(y)))
    stop("Arguments 'x' and 'y' must be real or complex vectors.")
  m <- length(x)
  if (length(y) != m)
    stop("Arguments 'x', 'y' must be vectors of the same length.")
  if (m <= 1)
    return(0)
  xp <- c(x, x[m:1])
  yp <- c(numeric(m), y[m:1])
  n <- 2 * m
  p1 <- sum(xp[1:(n - 1)] * yp[2:n]) + xp[n] * yp[1]
  p2 <- sum(xp[2:n] * yp[1:(n - 1)]) + xp[1] * yp[n]
  return(0.5 * (p1 - p2))
}

##.normalize01 <- function(x) { (x-min(x))/(max(x)-min(x)) }

##-------- aac for 1 model ------------------------------
.getAAC <- function(x, y)
{
  x <- .normalize01(x)
  y <- .normalize01(y)
  #aac <- 1- ( .trapz(x, y) )

  #x <- x-x[1]
  aac <- .trapz(x, y)
  return(aac)
}

##-----------------------------

## compute area above the curve (1-auc)
#'
#' compute area above the curve (1-auc)
#' @description
#' Given a model.id or batch it will return
#' aac (area above the curve)
#'
#' @examples
#' data(pdxe)
#' aac(pdxe, model.id = "X.015.BY19")
#' # creat a experiment desing
#' myDesign = list(batch.name="myBatch", treatment=c("X.015.BY19"), control=c("X.015.uned"))
#' aac(pdxe, batch = myDesign)
#'
#' @param object The \code{Xeva} dataset
#' @param model.id A model id
#' @param batch A list with treatment and control
#' @param average For batch. Either "curve" or "area". See details
#' @param treatment.only default {TRUE}. If TRUE only treatment periode will be considered
#'
#' @details If average = "curve", aac of average curve will be calculated. If average = "area", aac value for each curve will be calculated seperately and values will be average after.
#'
#' @return aac values for model.id or treatment and control aac for batch
setGeneric(name = "aac",
           def = function(object, model.id=NULL, batch=NULL, average = c("curve", "area"), treatment.only=TRUE)
           {standardGeneric("aac")} )

#' @export
setMethod( f=aac,
           signature=c(object="XevaSet"),
           definition= function(object, model.id=NULL, batch=NULL,
                                average = c("curve", "area"),
                                treatment.only=TRUE)
           {
             if(is.null(model.id) & is.null(batch))
             {stop("please specify 'model.id' or 'batch'")}

             average = average[1]

             ###---------------------------------------------------------------------------------------
             .aacFor1ModelId <- function(object, model.id, treatment.only)
             {
               DFx <- getExperiment(object, model.id=model.id, treatment.only = treatment.only)
               aac <- .getAAC(DFx$time, DFx$volume)
               return(aac)
             }
             ###---------------------------------------------------------------------------------------

             if(!is.null(model.id))
             {
               if(!is.null(batch))
               { warning("ignoring 'batch'")}
               aac <- .aacFor1ModelId(object, model.id, treatment.only)
               return(aac)
             }

             if(!is.null(batch))
             {
               rtx <- list(treatment=NA, control=NA)

               if(average=="curve")
               {
                 bx <- getTimeVarData(object, batch, treatment.only = treatment.only)
                 tr <- bx[bx$exp.type=="treatment",]
                 cr <- bx[bx$exp.type=="control",  ]
                 if(nrow(tr)>0){ rtx$treatment <- .getAAC(tr[ ,"time"], tr[ ,"mean"])}
                 if(nrow(cr)>0){ rtx$control   <- .getAAC(cr[ ,"time"], cr[ ,"mean"])}
                 return(rtx)
               }

               if(average=="area")
               {
                 tr <- sapply(batch$treatment, .aacFor1ModelId, object=object, treatment.only=treatment.only)
                 cn <- sapply(batch$control, .aacFor1ModelId, object=object, treatment.only=treatment.only)
                 if(length(tr)>0) {rtx$treatment <- mean(tr)}
                 if(length(cr)>0) {rtx$control   <- mean(cr)}
                 return(rtx)
               }
             }
           })





###########################################

#####================= setAAC ==================
#' setAAC generic
#'
#' Generic for setAAC method
#'
#' @examples
#' data(cm.pdxe)
#' setAAC(cm.pdxe) <- setAAC(object = cm.pdxe)
#' @param object The \code{XevaSet}
#' @return a \code{list} with model and batch AAC
setGeneric(name = "setAAC", def = function(object) {standardGeneric("setAAC")} )

#' @export
setMethod( f=setAAC, signature="XevaSet",
           definition=function(object)
           {
             value <- list()
             value$model <- sapply(modelInfo(object)$model.id,function(mid) {aac(object, model.id = mid)})

             value$batch <- lapply(names(object@expDesign),
                                   function(bid) {aac(object, batch = expDesign(object,bid))})
             names(value$batch) <- names(object@expDesign)

             return(value)
           } )


#' setAAC<- Generic
#'
#' Generic for setAAC replace method
#'
#' @examples
#' data(cm.pdxe)
#' setAAC(cm.pdxe) <- setAAC(object = cm.pdxe)
#' @param object The \code{XevaSet} to replace drug info in
#' @return Updated \code{XevaSet}
setGeneric(name= "setAAC<-", def = function(object, value) {standardGeneric("setAAC<-")} )

#' @export
setMethod( f="setAAC<-",
           signature=c(object = "XevaSet"),
           definition=function(object, value)
           {
             object@sensitivity$model[ , "AAC"] <- value$model[object@sensitivity$model$model.id]

             object@sensitivity$batch[, "AAC.treatment"]<- sapply(object@sensitivity$batch$batch.name,
                                                                  function(bid){value$batch[[bid]][["treatment"]]})
             object@sensitivity$batch[, "AAC.control"]  <- sapply(object@sensitivity$batch$batch.name,
                                                                  function(bid){value$batch[[bid]][["control"]]})
             return(object)
           } )
