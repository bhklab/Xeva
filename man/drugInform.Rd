% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/access_slot_drug.R
\name{drugInform}
\alias{drugInform}
\alias{drugInform,XevaSet-method}
\title{Get drug information
Get the drug information slot from a {XevaSet} object.}
\usage{
drugInform(object)

\S4method{drugInform}{XevaSet}(object)
}
\arguments{
\item{object}{The \code{XevaSet} to retrieve drug information from.}
}
\value{
A \code{data.frame} with the drug annotations.
}
\description{
Get drug information
Get the drug information slot from a {XevaSet} object.
}
\examples{
data(brca)
head(drugInform(brca))
}
