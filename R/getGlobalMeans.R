# Get the rowMeans of a matrix
computeGlobalMean <- function(mat) {
  globalMean <- matrix(rowMeans(mat, na.rm = TRUE), ncol = 1)
  colnames(globalMean) <- "globalMean"
  return(globalMean)
}

#' Get the global means of a matrix
#'
#' @name getGlobalMeans
#'
#' @param obj Input SummarizedExperiment object
#' @param targets Column names or indices to indicate which samples to shrink towards
#' @param assay What type of assay the data are from
#'
#' @return A vector of global or targeted means
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' scrna.global.means <- getGlobalMeans(k562_scrna_chr14, assay = "rna")
getGlobalMeans <- function(obj, targets = NULL, assay = c("atac", "rna", "array")) {
  # match the assay arg
  assay <- match.arg(assay)

  is.array <- assay == "array"

  # get the global means
  # check if shrinkage targets are being used
  if (!is.null(targets)) {
    stargets <- getShrinkageTargets(obj, targets)
    message("Using ", paste(shQuote(targets), collapse = ", "), " as shrinkage targets...")
    globalMean.input <- stargets
  } else {
    globalMean.input <- obj
  }

  assay.data <- .getAssay(globalMean.input, is.array)
  globalMean <- computeGlobalMean(assay.data)
  # coercion to get the rownames to be the GRanges coordinates
  # allows to flip back and forth between a matrix and GRanges for subsetting, etc.
  rownames(globalMean) <- as.character(granges(obj))
  return(globalMean)
}

#' Pre-compute the global means for bootstrapping compartments
#'
#' @name precomputeBootstrapMeans
#'
#' @param obj Input SummarizedExperiment object
#' @param targets Optional targets to shrink towards
#' @param num.bootstraps The number of bootstraps to compute
#' @param assay What type of assay the data are from
#' @param parallel Whether to run in parallel
#' @param num.cores How many cores to use for parallel processing
#'
#' @return A matrix of bootstrapped global means
#'
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data("k562_scrna_chr14", package = "compartmap")
#' scrna.bootstrap.global.means <- precomputeBootstrapMeans(
#'   k562_scrna_chr14,
#'   assay = "rna",
#'   num.bootstraps = 2
#' )
precomputeBootstrapMeans <- function(
  obj,
  targets = NULL,
  num.bootstraps = 100,
  assay = c("atac", "rna", "array"),
  parallel = FALSE,
  num.cores = 1
) {
  # this function precomputes the bootstrapped global means
  # as a default we will make 100 bootstraps
  assay <- match.arg(assay)
  is.array <- assay == "array"

  if (!is.null(targets)) {
    if (length(targets) < 5) stop("Need more than 5 samples for targeted bootstrapping to work.")
    obj <- getShrinkageTargets(obj, targets)
  }
  bootMean <- mclapply(1:num.bootstraps, function(b) {
    message("Working on bootstrap ", b)
    assay.data <- .getAssay(obj, is.array)
    resamp.mat <- .resampleMatrix(assay.data)
    computeGlobalMean(resamp.mat)
  }, mc.cores = ifelse(parallel, num.cores, 1))

  bootResult <- do.call("cbind", bootMean)
  rownames(bootResult) <- as.character(granges(obj))
  return(bootResult)
}

# Get $counts of $Beta depening on whether the input is an array experiment or rna/atac
.getAssay <- function(obj, is.array) {
  assay.name <- ifelse(is.array, "Beta", "counts")

  if (!assay.name %in% names(assays(obj))) {
    msg <- paste(shQuote(assay.name), "not found in the input object's assays")
    stop(msg)
  }

  assay.data <- assays(obj)[[assay.name]]
  if (is.array) {
    assay.data <- flogit(assay.data)
  }
  return(assay.data)
}
