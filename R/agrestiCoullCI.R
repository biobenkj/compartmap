#' Normal alpha/2 quantile
#' 
#' @param q   the quantile at which to extract Z
#' 
#' @return    Z
#'
#' @details
#' \eqn{z_\alpha = \Phi^{-1}(1 - \frac{\alpha}{2})}
#' @param q   the quantile at which to extract Z
#'
#' @return    Z
#' @keywords internal
.z <- function(q) {
  qnorm(1 - ((1 - q) / 2))
}

#' n_tilde in AC
#'
#' @details
#' \eqn{\tilde{n} = n_{\text{successes}} + n_{\text{failures}} + z^2_\alpha}
#'
#' @param n1  number of successes/ones 
#' @param n0  number of failures/zeroes
#' @param q   quantile for eventual CI (e.g. 0.95 for a 95 percent binomial CI)
#'
#' @return    the effective sample size for smoothed CIs
#' @keywords internal
.n_approx <- function(n1, n0, q) {
  (n0 + n1) + (.z(q)^2)
}

#' p_tilde in AC 
#'
#' @details
#' \eqn{\tilde{p} = \frac{1}{\tilde{n}}(n_{\text{success}} + \frac{z^2_\alpha}{2})}
#'
#' @param n1  number of successes/ones 
#' @param n0  number of failures/zeroes
#' @param q   quantile for eventual CI (e.g. 0.95 for a 95 percent binomial CI)
#'
#' @return    the approximate success probability for a smoothed CIs
#' @keywords internal
.p_approx <- function(n1, n0, q) {
  (1/.n_approx(n1, n0, q)) * (n1 + ((.z(q)^2)/2))
}

#' Agresti-Coull confidence interval for a binomial proportion
#'
#' @details
#' \eqn{z_\alpha = \Phi^{-1}(1 - \frac{\alpha}{2})}
#'
#' \eqn{\tilde{n} = n_{\text{successes}} + n_{\text{failures}} + z^2_\alpha}
#'
#' \eqn{\tilde{p} = \frac{1}{\tilde{n}}(n_{\text{success}} + \frac{z^2_\alpha}{2})}
#'
#' \eqn{p \approx \tilde{p} \pm z_\alpha \times \sqrt{\frac{\tilde{p}}{\tilde{n}} \times (1 - \tilde{p})}}
#'
#' @param n1  number of successes/ones 
#' @param n0  number of failures/zeroes
#' @param q   quantile for eventual CI (e.g. 0.95 for a 95 percent binomial CI)
#'
#' @return    the approximate (q x 100) percent confidence interval for (p|n1,n0,q)
#' @export
#'
#' @examples 
#' agrestiCoullCI(10, 3, 0.95)
agrestiCoullCI <- function(n1, n0, q) {
  p_apx <- .p_approx(n1, n0, q)
  n_apx <- .n_approx(n1, n0, q)
  width <- .z(q) * sqrt( (p_apx/n_apx) * (1 - p_apx) )
  data.frame(
    conf.est = p_apx,
    conf.est.lowerCI = pmax(0, (p_apx - width)),
    conf.est.upperCI = pmin((p_apx + width), 1)
  )
}
