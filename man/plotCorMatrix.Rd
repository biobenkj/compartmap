% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMT.R
\name{plotCorMatrix}
\alias{plotCorMatrix}
\title{Plot a denoised correlation matrix}
\usage{
plotCorMatrix(
  denoised.cor.mat,
  midpoint = 0.3,
  return.plot.obj = FALSE,
  uppertri = FALSE,
  lowertri = FALSE
)
}
\arguments{
\item{denoised.cor.mat}{The denoised correlation matrix object from getDenoisedMatrix}

\item{midpoint}{The midpoint for the coloring (default is 0.3)}

\item{return.plot.obj}{Whether to return the ggplot object}

\item{uppertri}{Whether to keep the upper triangle of the matrix}

\item{lowertri}{Whether to keep the lower triangle of the matrix}
}
\value{
Either a ggplot object or plot
}
\description{
Plot a denoised correlation matrix
}
\examples{
dummy <- matrix(rnorm(10000), ncol=25)
set.seed(1000)
my_plot <- plotCorMatrix(dummy, return.plot.obj = TRUE)
}
