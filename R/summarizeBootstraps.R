#' Summarize the bootstrap compartment estimates and compute Agresti-Coull confidence intervals
#'
#' @name summarizeBootstraps
#'
#' @param boot.list List of bootstraps as GRanges objects
#' @param est.ab The original compartment calls
#' @param q Confidence interval to compute (0.95 for 95 percent CI)
#' @param assay Type of assay we are working with
#'
#' @return A GRanges object with bootstraps summarized
#' @import SummarizedExperiment
#' @export
#'
#' @examples
summarizeBootstraps <- function(boot.list, est.ab, q = 0.95, assay = c("rna", "atac")) {
  # go through the estimated A/B compartments and compute proportions from the boot.list
  est.ab$score <- est.ab$pc

  is.atac_or_rna <- assay %in% c("atac", "rna")

  message("Summarizing bootstraps.")

  # initialize open and closed counts to enumerate
  # est.ab$boot.open <- 0
  # est.ab$boot.closed <- 0
  # filter out failed compartment estimates
  boot.list <- removeEmptyBoots(boot.list)

  # summarize
  # create a dummy matrix and then rowSum them up
  boot.summary.mat.lst <- lapply(boot.list, function(b) {
    # add the pc to the granges object
    b$score <- b$pc

    # generate a dummy GRanges object
    est.ab.dummy <- est.ab
    est.ab.dummy$boot.open <- 0
    est.ab.dummy$boot.closed <- 0

    # determine whether compartment is open and convert the boolean to 1/0 binary result for proportions
    b.isOpen <- .isCompartmentOpen(is.atac_or_rna, b$score)
    b$open <- as.integer(b.isOpen)
    b$closed <- as.integer(!b.isOpen)

    # overlap by common intervals
    ol <- findOverlaps(b, est.ab.dummy)
    mcols(est.ab.dummy)$boot.open[subjectHits(ol)] <- mcols(est.ab.dummy)$boot.open[subjectHits(ol)] + mcols(b)$open[queryHits(ol)]
    mcols(est.ab.dummy)$boot.closed[subjectHits(ol)] <- mcols(est.ab.dummy)$boot.closed[subjectHits(ol)] + mcols(b)$closed[queryHits(ol)]

    # return the dummy mcols for bootstrapped open and closed calls
    return(as.matrix(cbind(
      mcols(est.ab.dummy)$boot.open,
      mcols(est.ab.dummy)$boot.closed
    )))
  })

  # eunumerate bootstraps
  est.ab$boot.open <- rowSums(do.call("cbind", lapply(boot.summary.mat.lst, function(x) x[, 1])))
  est.ab$boot.closed <- rowSums(do.call("cbind", lapply(boot.summary.mat.lst, function(x) x[, 2])))

  message("Computing Agresti-Coull 95% confidence intervals.")
  est.ab$conf.est <- 0
  est.ab$conf.est.upperCI <- 0
  est.ab$conf.est.lowerCI <- 0

  conf.int <- lapply(1:length(est.ab), function(e) {
    compartment.call <- est.ab[e, ]
    ab.score <- compartment.call$score

    # check if the compartment is open
    is.open <- .isCompartmentOpen(is.atac_or_rna, ab.score)

    # get ones and zeroes input for agrestiCoullCI
    ones <- ifelse(is.open, compartment.call$boot.open, compartment.call$boot.closed)
    zeroes <- ifelse(is.open, compartment.call$boot.closed, compartment.call$boot.open)
    agrestiCoullCI(ones, zeroes, q = 0.95)
  }

  # combine the conf.est results into something sensible and bind with est.ab
  conf.int.ests <- do.call("rbind", conf.int)
  # should be of the form:
  # conf.est.lowerCI conf.est conf.est.upperCI
  est.ab$conf.est.lowerCI <- conf.int.ests[, 1]
  est.ab$conf.est <- conf.int.ests[, 2]
  est.ab$conf.est.upperCI <- conf.int.ests[, 3]
  return(est.ab)
}

# Check if a compartment is open based on assay type and eigenvalue
#
# For ATAC/RNA:
# eigen < 0 - closed
# eigen > 0 - open
#
# For methylation the logic is flipped:
# eigen < 0 - open
# eigen > 0 - closed
.isCompartmentOpen <- function(is.atac_or_rna, eigen) {
  (is.atac_or_rna & eigen > 0) | (!is.atac_or_rna & eigen < 0)
}
