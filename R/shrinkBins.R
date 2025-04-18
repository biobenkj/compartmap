#' @title Employ an eBayes shrinkage approach for bin-level estimates for A/B inference
#'
#' @description
#' \code{shrinkBins} returns shrunken bin-level estimates
#'
#' @details This function computes shrunken bin-level estimates using a James-Stein estimator, reformulated as an eBayes procedure
#'
#' @param x Input SummarizedExperiment object
#' @param original.x Full sample set SummarizedExperiment object
#' @param prior.means The means of the bin-level prior distribution
#' @param chr The chromosome to operate on
#' @param res Resolution to perform the binning
#' @param targets The column/sample/cell names to shrink towards
#' @param jse Whether to use a James-Stein estimator (default is TRUE)
#' @param assay What assay type this is ("rna", "atac", "array")
#' @param genome What genome are we working with ("hg19", "hg38", "mm9", "mm10")
#'
#' @return A list object to pass to getCorMatrix
#'
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom stats median sd
#'
#' @export
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' shrunken.bin.scrna <- shrinkBins(
#'   x = k562_scrna_chr14,
#'   original.x = k562_scrna_chr14,
#'   chr = "chr14", assay = "rna"
#' )
#'
shrinkBins <- function(
  x,
  original.x,
  prior.means = NULL,
  chr = NULL,
  res = 1e6,
  targets = NULL,
  jse = TRUE,
  assay = c("rna", "atac", "array"),
  genome = c("hg19", "hg38", "mm9", "mm10")
) {
  # match the assay args
  assay <- match.arg(assay)

  # match the genome if given
  genome <- match.arg(genome)

  # get the prior means
  if (is.null(prior.means)) {
    prior.means <- getGlobalMeans(obj = original.x, targets = targets, assay = assay)
  }

  # helper function for summary
  # not used if JSE is set to TRUE
  atac_fun <- function(x) {
    return(sqrt(mean(x)) * length(x))
  }

  is.atac_or_rna <- assay %in% c("atac", "rna")
  input.fun <- if (jse) {
    mean
  } else if (is.atac_or_rna) {
    atac_fun
  } else {
    median
  }

  input.assay <- if (is.atac_or_rna) {
    assays(original.x)$counts
  } else {
    flogit(assays(original.x)$Beta) # make sure we are with betas or we will double flogit
  }

  # bin the input
  bin.mat <- getBinMatrix(
    mat = as.matrix(cbind(input.assay, prior.means)),
    genloc = rowRanges(x),
    chr = chr,
    res = res,
    FUN = input.fun,
    genome = genome
  )

  # shrink the bins using a James-Stein Estimator
  x.shrink <- t(apply(bin.mat$x, 1, function(r) {
    r.samps <- r[!names(r) %in% "globalMean"]
    r.prior.m <- r["globalMean"]

    if (!is.null(targets) & length(r.samps[targets]) == 1) {
      stop("Cannot perform targeted bin-level shrinkage with one target sample.")
    }

    if (jse) {
      .jse(x = r.samps, grand.mean = r.prior.m, targets = targets)
    } else {
      .ebayes(x = r.samps, prior = r.prior.m, targets = targets)
    }
  }))

  # drop things that are zeroes as global means
  # this can and does crop up in resampling when you have something sparse
  # for instance single-cell data...
  # the correlation will break otherwise
  if (any(bin.mat$x[, "globalMean"] == 0)) {
    bin.mat$gr <- bin.mat$gr[bin.mat$x[, "globalMean"] != 0, ]
    x.shrink <- x.shrink[bin.mat$x[, "globalMean"] != 0, ]
    bin.mat$x <- bin.mat$x[bin.mat$x[, "globalMean"] != 0, ]
  }

  return(list(gr = bin.mat$gr, x = x.shrink[, colnames(x)], gmeans = bin.mat$x[, "globalMean"]))
}

# helper functions for computing shrunken means
.ebayes <- function(x, prior = NULL, targets = NULL) {
  if (!is.null(targets)) {
    C <- sd(x[targets])
  } else {
    C <- sd(x)
  }

  # convert back to beta values
  return(prior + C * (x - prior))
}

.jse <- function(x, grand.mean = NULL, targets = NULL) {
  ## see if we have enough means...
  ## this also assumes we are just computing a straight mean
  if (is.null(targets)) {
    ## typical shrinkage
    c <- 1 - ((length(x) - 3) * (sd(x)^2) / sum((x - grand.mean)^2))
  } else if (length(targets) < 4) {
    message("Number of means fewer than 4. Using Bayes instead.")
    ## this falls back to using Bayes rule which will probably not be great
    ## but it won't explode and may provide some reasonable results anyway
    c <- sd(x[targets])
  } else {
    ## targeted shrinkage
    c <- 1 - ((length(x) - 3) * (sd(x[targets])^2) / sum((x - grand.mean)^2))
  }

  return(grand.mean + c * (x - grand.mean))
}
