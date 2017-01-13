#####================= getMolecularProfiles ==================
#' Get Molecular Profiles
#'
#' Get Molecular Profiles
#'
#' @param object The \code{XevaSet}
#' @param data.type \code{character}, which one of the molecular data types is needed
#' @return a \code{ExpressionSet} where sample names are \code{biobase.id} of model
#' @examples
#' data(pdxe)
#' pdxe_RNA <- getMolecularProfiles(pdxe, data.type="RNASeq")
#' @export
getMolecularProfiles <- function(object, data.type)
{
  if(is.element(data.type, names(object@molecularProfiles))==FALSE)
  {
    msg = sprintf("avalibale expression sets are\n%s\n",
                  paste(names(object@molecularProfiles), collapse ="\n"))
    stop(msg)
  }
  expset = object@molecularProfiles[[data.type]]
  return(expset)
}
