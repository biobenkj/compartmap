#' Filter compartments using confidence estimates and eigenvalue thresholds
#'
#' @name filterCompartments
#'
#' @param obj Output of condenseSE or fixCompartments
#' @param min.conf Minimum confidence estimate to use when filtering
#' @param min.eigen Minimum absolute eigenvalue to use when filtering
#'
#' @return A filtered/subset of the input object/list
#' @import SummarizedExperiment
#' @export
#'
#' @examples
filterCompartments <- function(obj, min.conf = 0.7, min.eigen = 0.02) {
  # filter compartments
  message("Filtering compartments based on a minimum confidence of ", min.conf * 100, "%")
  message("Filtering compartments based on a minimum absolute eigen value of ", min.eigen)
  if (is(obj, "list")) {
    lapply(obj, function(x) {
      filt <- filterer(x, min.conf, min.eigen)
      return(x[filt, ])
    })
  } else {
    filt <- filterer(obj, min.conf, min.eigen)
    return(obj[filt, ])
  }
}

filterer <- function(obj, min.conf, min.eigen) {
  apply(mcols(obj), 1, function(x) {
    # check if we have "fixed" things
    flip <- "flip.score" %in% names(x)
    conf.name <- ifelse(flip, "filp.conf.est", "conf.est")
    score.name <- ifelse(flip, "filp.score", "score")
    filt.score <- as.numeric(x[conf.name]) >= min.conf & abs(as.numeric(x[score.name])) >= min.eigen
    return(filt.score)
  })
}
