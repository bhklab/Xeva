.checkExperimentDesign<- function(expDesign)
{
  modNoControl = c()
  modNoTreatme = c()
  for(I in expDesign)
  {
    if(length(I$control)==0 & length(I$treatment)==0)
    {stop("Error Treatmetn and Control are missing in expDesign!")}

    if(length(I$control)==0)  { modNoControl = c(modNoControl, I$treatment)}
    if(length(I$treatment)==0){ modNoTreatme = c(modNoTreatme, I$control)}
  }


  if(!is.null(modNoControl))
  {
    txt = sprintf("These models have no Controls\n%s\n", paste(unique(modNoControl), collapse = "\n"))
    cat(txt)
  }

  ##------- setting name -----------------
  ###bnam <- xapply(expDesign, "[[", "batch.name")
  bnam <- c()
  for(i in expDesign){bnam <- c(bnam, i[["batch.name"]])}

  bnamDup <- bnam[duplicated(bnam)]
  if(length(bnamDup)>0)
  {
    txt <- sprintf("These batch names are duplicated\n%s\n", paste(bnamDup, collapse = "\n"))
    stop(txt)
  }
  names(expDesign) <- bnam

  for(i in names(expDesign))
  {
    class(expDesign[[i]]) <- append(class(expDesign[[i]]),"pdxBatch")
  }

  return(expDesign)
}

#' Print the pdx batch
#' @param x pdxBatch object
#' @param ... Other arguments
#' @return prints pdxBatch
#' @export
print.pdxBatch <- function(x, ...)
{
  if(is.null(x$control))  { x$control  <- "NA"}
  if(is.null(x$treatment)){ x$treatment<- "NA"}
  txt <- sprintf("name = %s\ncontrol = %s\ntreatment = %s\n", x$batch.name,
                 paste0(x$control, collapse = ", "),
                 paste0(x$treatment, collapse = ", "))
  for(i in names(x))
  {
    if(!i %in% c("batch.name", "control", "treatment"))
    { txt <- sprintf("%s%s = %s\n", txt, i, paste0(x[[i]], collapse = ", ")) }
  }
  cat(txt)
}
