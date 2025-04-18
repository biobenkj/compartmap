#' Compute the SVD of a matrix using irlba
#'
#' @name getSVD
#'
#' @param matrix A p x n input matrix
#' @param k Number of singular vectors to return
#' @param sing.vec Whether to return the right or left singular vector
#'
#' @return A singular vector or matrix with sign corresponding to positive values
#'
#' @importFrom stats cor
#' @importFrom BiocSingular IrlbaParam runSVD
#' @export
#'
#' @examples
#'
#' dummy <- matrix(rnorm(10000), ncol=25)
#' sing_vec <- getSVD(dummy, k = 1, sing.vec = "right")
getSVD <- function (matrix, k = 1, sing.vec = c("left", "right")) {
  sing.vec <- match.arg(sing.vec)

  #center the matrix
  matrix <- t(scale(t(matrix), center = TRUE, scale = FALSE))

  p.mat <- .getUV(matrix, k, sing.vec)
  csums <- colSums(matrix, na.rm=TRUE)

  # check for negative correlation
  # flip sign as necessary to ensure signal is associated with pos. values
  if (cor(csums, p.mat) < 0) {
    p.mat <- -p.mat
  }

  p.mat * sqrt(length(p.mat)) # Chromosome length normalization
}

#' Run SVD and get left or right singular vectors
#'
#' @param matrix A p x n input matrix
#' @param k Number of singular vectors to return
#' @param sing.vec Whether to return the right or left singular vector
#'
#' @return left or right singular vectors
#'
#' @keywords internal
.getUV <- function(matrix, k, sing.vec) {
  SVD <- runSVD(matrix, k=k, BSPARAM=IrlbaParam())
  switch(sing.vec,
    left = SVD$u,
    right = SVD$v
  )
}
