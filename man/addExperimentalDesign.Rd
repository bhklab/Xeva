% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expDesign_manuplation.R
\name{addExperimentalDesign}
\alias{addExperimentalDesign}
\alias{addExperimentalDesign,XevaSet-method}
\title{Add a new experimental design}
\usage{
addExperimentalDesign(
  object,
  treatment = NULL,
  control = NULL,
  batch.id = NULL,
  replace = FALSE
)

\S4method{addExperimentalDesign}{XevaSet}(
  object,
  treatment = NULL,
  control = NULL,
  batch.id = NULL,
  replace = FALSE
)
}
\arguments{
\item{object}{The \code{Xeva} dataset.}

\item{treatment}{The \code{model.id} of treatment.}

\item{control}{The \code{model.id} of control.}

\item{batch.id}{The \code{batch.id} for a new batch.}

\item{replace}{If \code{TRUE}, replace an old batch with new values.}
}
\value{
Returns \code{Xeva} dataset with new experimental design added.
}
\description{
Add a new experimental design in the \code{expDesign} slot.
}
\examples{
data(brca)
brca <- addExperimentalDesign(object=brca, treatment=c("X.6047.LL71"),
        control=c("X.6047.uned"), batch.id="new.batch", replace=FALSE)

}
