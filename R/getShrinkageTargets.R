#' Get the specified samples to shrink towards instead of the global mean
#'
#' @name getShrinkageTargets
#'
#' @param obj Input matrix
#' @param group Sample names, colnames or column indices to use for targeted
#' shrinkage
#'
#' @return A matrix composed of samples to shrink towards
#' @export
#'
#' @examples
#' dummy <- matrix(rnorm(1000), ncol=25)
#' dummy.sub <- getShrinkageTargets(dummy, group = c(1,5,8,10))
getShrinkageTargets <- function(obj, group) {
  tryCatch(
    obj[, unique(group)],
    error = function(e) {
      no_colnames <- is.null(colnames(obj))
      group.is_colnames <- is.character(group)
      column.type <- ifelse(group.is_colnames, "names", "indices")

      msg <- paste(
        "Error while subsetting targets: provided column", column.type, "not found:"
      )
      if (no_colnames & group.is_colnames) {
        message("The provided object does not have any column names - use column indices instead.")
      }
      message(msg)
      stop(e)
    }
  )
}
