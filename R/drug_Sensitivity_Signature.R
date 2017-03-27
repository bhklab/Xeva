drugSensitivitySig <- function(object, mDataType, drugs, features,
                               sensitivity.measure="slop",
                               molecular.summary.stat=c("mean", "median", "first", "last", "or", "and"),
                               sensitivity.summary.stat=c("mean", "median", "first", "last"),
                               returnValues=c("estimate", "pvalue", "fdr"),
                               sensitivity.cutoff, standardize=c("SD", "rescale", "none"), nthread=1, verbose=TRUE, ...)

{
  mDType = "RNASeq"
  z = getMolecularProfiles(object, mDType)

  rr <- PharmacoGx:::rankGeneDrugSensitivity(data=expr, drugpheno=dd, type=type,
                                             batch=batch, single.type=FALSE,
                                             standardize=standardize, nthread=nthread,
                                             verbose=verbose)


}




