#' @title Employ an eBayes shrinkage approach for bin-level estimates for A/B inference
#'
#' @description
#' \code{shrinkBins} returns shrunken bin-level estimates
#'
#' @details This function computes shrunken bin-level estimates using a
#' James-Stein estimator (JSE), reformulated as an eBayes procedure. JSE can be
#' used only if at least 4 targets are provided - any less and `shrinkBins`
#' will fall back to using Bayes rule which will probably not be great but it
#' won't explode and may provide some reasonable results anyway
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
  assay <- match.arg(assay)
  genome <- match.arg(genome)
  verifySE(x)

  target.count <- length(targets)
  if (target.count == 1) {
    stop("Cannot perform targeted bin-level shrinkage with one target sample.")
  } else if (target.count < 4) {
    message("Number of means fewer than 4. Using Bayes instead of JSE.")
    jse <- FALSE
  }

  # get the prior means
  prior.means <- prior.means %||% getGlobalMeans(
    obj = original.x,
    targets = targets,
    assay = assay
  )

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

# helper function for summary when JSE == FALSE
atac_fun <- function(x) {
  sqrt(mean(x)) * length(x)
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

#' James-Stein estimator
#' @param x input vector of binned means across samples
#' @param grand.mean The global mean across samples
#' @param targets Samples to shrink towards
#'
#' \eqn{\hat{\theta}_{JS+} = \left(1 - \frac{(m - 3)\sigma^2}{||\textbf{y} - \nu||^2}\right)}
.jse <- function(x, grand.mean = NULL, targets = NULL) {
  ## see if we have enough means...
  ## this also assumes we are just computing a straight mean
  if (is.null(targets)) {
    ## typical shrinkage
    c <- 1 - ((length(x) - 3) * (sd(x)^2) / sum((x - grand.mean)^2))
  } else {
    ## targeted shrinkage
    c <- 1 - ((length(x) - 3) * (sd(x[targets])^2) / sum((x - grand.mean)^2))
  }

  return(grand.mean + c * (x - grand.mean))
}
