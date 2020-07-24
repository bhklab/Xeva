.onAttach <- function(libname, pkgname)
{
  version = packageDescription(pkgname, fields = "Version")
  msg = paste0("#-----------*****-----------
", pkgname, " version ", version, "
Bioconductor: http://bioconductor.org/packages/Xeva
      Github: https://github.com/bhklab/Xeva

Citation:
Mer, A.S., et al. Integrative pharmacogenomics analysis of patient-derived xenografts. Cancer Research 79.17 (2019): 4539-4550

DOI: http://dx.doi.org/10.1158/0008-5472.CAN-19-0349

To suppress this message use:
suppressPackageStartupMessages(library(Xeva))
#-----------*****-----------
")

  packageStartupMessage(msg)
}
