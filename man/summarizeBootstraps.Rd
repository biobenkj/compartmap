% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarizeBootstraps.R
\name{summarizeBootstraps}
\alias{summarizeBootstraps}
\title{Summarize the bootstrap compartment estimates and compute Agresti-Coull confidence intervals}
\usage{
summarizeBootstraps(boot.list, est.ab, q = 0.95, assay = c("rna", "atac"))
}
\arguments{
\item{boot.list}{List of bootstraps as GRanges objects}

\item{est.ab}{The original compartment calls}

\item{q}{Confidence interval to compute (0.95 for 95 percent CI)}

\item{assay}{Type of assay we are working with}
}
\value{
A GRanges object with bootstraps summarized
}
\description{
Summarize the bootstrap compartment estimates and compute Agresti-Coull confidence intervals
}
