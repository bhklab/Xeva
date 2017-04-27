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
<<<<<<< HEAD
  if(is.element(data.type, names(slot(object, "molecularProfiles")))==FALSE)
=======
  if(is.element(data.type, names(object@molecularProfiles))==FALSE)
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
  {
    msg = sprintf("available molecular data are\n%s\n",
                  paste(names(object@molecularProfiles), collapse ="\n"))
    stop(msg)
  }
<<<<<<< HEAD
  expset <- slot(object, "molecularProfiles")[[data.type]]
=======
  expset = object@molecularProfiles[[data.type]]
>>>>>>> 9f9947748d00443b9546698266dd7eb78c636ce4
  return(expset)
}
