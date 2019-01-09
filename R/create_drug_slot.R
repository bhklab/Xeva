.checkDrugSlot <- function(drf)
{
  #drug.id   standard.name
  drf <- data.frame(apply(drf, 2, as.character), stringsAsFactors = FALSE)
  return(drf)
}
