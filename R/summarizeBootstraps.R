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
#'
#' @import SummarizedExperiment
#' @importFrom S4Vectors subjectHits queryHits
#'
#' @export
#'
#' @examples
summarizeBootstraps <- function(boot.list, est.ab, q = 0.95, assay = c("rna", "atac")) {
  # go through the estimated A/B compartments and compute proportions from the boot.list
  est.ab$score <- est.ab$pc

  is.atac_or_rna <- assay %in% c("atac", "rna")

  message("Summarizing bootstraps.")

  # filter out failed compartment estimates
  boot.list <- removeEmptyBoots(boot.list)

  # summarize
  # create a dummy matrix and then rowSum them up
  boot.summary.mat.lst <- lapply(boot.list, .getSummary, est.ab = est.ab)

  # eunumerate bootstraps
  .getBootRowSums <- function(index) {
    rowSums(do.call("cbind", lapply(boot.summary.mat.lst, function(x) x[, index])))
  }
  est.ab$boot.open <- .getBootRowSums(1)
  est.ab$boot.closed <- .getBootRowSums(2)

  message("Computing Agresti-Coull 95% confidence intervals.")
  est.ab$conf.est <- 0
  est.ab$conf.est.upperCI <- 0
  est.ab$conf.est.lowerCI <- 0

  conf.int <- lapply(1:length(est.ab), .getCI, q = q)

  # combine the conf.est results into something sensible and bind with est.ab
  conf.int.ests <- do.call("rbind", conf.int)
  # should be of the form:
  # conf.est.lowerCI conf.est conf.est.upperCI
  est.ab$conf.est.lowerCI <- conf.int.ests[, 1]
  est.ab$conf.est <- conf.int.ests[, 2]
  est.ab$conf.est.upperCI <- conf.int.ests[, 3]
  return(est.ab)
}

.getSummary <- function(gr.boot, est.ab) {
  # add the pc to the granges object
  gr.boot$score <- gr.boot$pc

  # generate a dummy GRanges object
  est.ab.dummy <- est.ab
  est.ab.dummy$boot.open <- 0
  est.ab.dummy$boot.closed <- 0

  # determine whether compartment is open and convert the boolean to 1/0 binary result for proportions
  gr.boot.isOpen <- .isCompartmentOpen(is.atac_or_rna, gr.boot$score)
  gr.boot$open <- as.integer(gr.boot.isOpen)
  gr.boot$closed <- as.integer(!gr.boot.isOpen)

  # overlap by common intervals
  ol <- findOverlaps(gr.boot, est.ab.dummy)
  mcols(est.ab.dummy)$boot.open[subjectHits(ol)] <- mcols(est.ab.dummy)$boot.open[subjectHits(ol)] + mcols(gr.boot)$open[queryHits(ol)]
  mcols(est.ab.dummy)$boot.closed[subjectHits(ol)] <- mcols(est.ab.dummy)$boot.closed[subjectHits(ol)] + mcols(gr.boot)$closed[queryHits(ol)]

  # return the dummy mcols for bootstrapped open and closed calls
  return(as.matrix(cbind(
    mcols(est.ab.dummy)$boot.open,
    mcols(est.ab.dummy)$boot.closed
  )))
}

.getCI <- function(gr.row, est.ab, q) {
  compartment.call <- est.ab[gr.row, ]
  ab.score <- compartment.call$score

  # check if the compartment is open
  is.open <- .isCompartmentOpen(is.atac_or_rna, ab.score)

  # get ones and zeroes input for agrestiCoullCI
  ones <- ifelse(is.open, compartment.call$boot.open, compartment.call$boot.closed)
  zeroes <- ifelse(is.open, compartment.call$boot.closed, compartment.call$boot.open)
  agrestiCoullCI(ones, zeroes, q = q)
}
