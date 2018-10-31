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
  bnam <- sapply(expDesign, "[[", "batch.name")
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

#' @export
print.pdxBatch <- function(b)
{
  if(is.null(b$control))  { b$control  <- "NA"}
  if(is.null(b$treatment)){ b$treatment<- "NA"}
  txt <- sprintf("name = %s\ncontrol = %s\ntreatment = %s\n", b$batch.name,
                 paste0(b$control, collapse = ", "),
                 paste0(b$treatment, collapse = ", "))
  cat(txt)
}






