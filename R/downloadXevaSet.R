availableXevaSet <- function()
{
  path <- "https://www.pmgenomics.ca/bhklab/sites/default/files/downloads/XevaSets/XevaSet.csv"
  x <- read.csv(path, sep = ",", stringsAsFactors = FALSE)
  return(x)
}

#' @import downloader
getXevaSet <- function(dw.url, saveDir = file.path(".", "XevaSet"),
                       XevaSetFileName = NULL, verbose = TRUE)
{
  if (!file.exists(saveDir))
  {
    dir.create(saveDir, recursive = TRUE)
  }
  if (is.null(XevaSetFileName))
  { XevaSetFileName <- rev(strsplit(dw.url, "/")[[1]])[1] }

  destfile=file.path(saveDir,XevaSetFileName)
  downloader::download(url = dw.url, destfile=destfile,
                       quiet = !verbose)

  xs <- readRDS(destfile)
  return(xs)

}


#' Download a XevaSet object or table of available XevaSet objects
#'
#' @description This function allows you to see the available XevaSet object and
#' download them for use with this package. The XevaSet have been extensively
#' curated and organised within a XevaSet class, enabling use with all the
#' analysis tools provided in Xeva.
#'
#' @examples
#' # downloadXevaSet()
#' ### to download a dataset
#' # library(Xeva)
#' # PDXE_BRCA = downloadXevaSet(name="PDXE_BRCA", saveDir="XevaSet")
#'
#' @param name Character string, the name of the XevaSet to download.
#' @param saveDir	\code{Character} string with the folder path where the XevaSet should be saved. Defaults to './XevaSet/'. Will create directory if it does not exist.
#' @param XevaSetFileName \code{character} string, the file name to save the dataset under
#' @param verbose	\code{bool} Should status messages be printed during download. Defaults to TRUE.
#' @return A data.frame if name is NULL, showing all the available XevaSet objects. If name is specified, it will download the dataset from our server
#' @export
downloadXevaSet <- function(name=NULL, saveDir = file.path(".", "XevaSet"),
                            XevaSetFileName = NULL, verbose = TRUE)
{
  axs <- availableXevaSet()
  if(!is.null(name))
  {
    if(name %in% axs$XevaSet.Name)
    {
      dw.url <- axs[axs$XevaSet.Name==name, "URL"]
      rtx <- getXevaSet(dw.url, saveDir = saveDir,
                        XevaSetFileName = XevaSetFileName, verbose = verbose)
      return(rtx)
    } else
    {
      stop("Unknown dataset name. Please use the downloadXevaSet() function for the table of available XevaSet.")
    }
  }
  return(axs)
}
