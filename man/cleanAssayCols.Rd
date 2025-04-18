% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{cleanAssayCols}
\alias{cleanAssayCols}
\title{Remove columns/cells/samples with NAs exceeding a threshold}
\usage{
cleanAssayCols(se, colmax = 0.8, assay = c("array", "bisulfite"))
}
\arguments{
\item{se}{Input SummarizedExperiment object}

\item{colmax}{The maximum number of NAs allowed as a fraction}

\item{assay}{The type of assay we are working with}
}
\value{
A filtered matrix
}
\description{
Remove columns/cells/samples with NAs exceeding a threshold
}
\examples{
if (requireNamespace("minfi", quietly = TRUE)) {
  data("array_data_chr14", package = "compartmap")
  cleanAssayCols(array.data.chr14, assay = "array")
}
}
