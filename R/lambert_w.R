#' Evaluation of the Lambert W Function (Branch -1)
#'
#' @description
#' Evaluates the non-principal (lower) branch \eqn{W_{-1}(z)} of the
#' Lambert W function for \eqn{z} in \eqn{[-1/e, 0)}, i.e., the unique
#' solution \eqn{w \le -1} of \eqn{w e^w = z}. Fully vectorized: every
#' element of \code{z} is refined in parallel via array operations, with
#' already-converged elements skipped on each iteration.
#'
#' @details
#' \code{\link{qautorelevate}} requires \eqn{W_{-1}} to invert the
#' autorelevated cumulative distribution function in closed form (no
#' baseline-specific quantile function is needed), and the
#' \code{"powerlindley"} baseline (see \link{autorelevate-package}) uses
#' it again internally for its own quantile inversion. The evaluation
#' proceeds in two steps:
#' \enumerate{
#'   \item \strong{Starting value.} For \eqn{z} close to \eqn{0^-}, the
#'   classical asymptotic series (see e.g. Corless, Gonnet, Hare, Jeffrey,
#'   & Knuth, 1996) gives \eqn{w_0 = L_1 - L_2 + L_2/L_1}, where
#'   \eqn{L_1 = \log(-z)} and \eqn{L_2 = \log(-L_1)}.
#'   \item \strong{Refinement.} \eqn{w_0} is refined by Halley's method
#'   applied to \eqn{f(w) = w e^w - z}, which converges cubically and is
#'   markedly more robust near the branch point \eqn{z = -1/e} than plain
#'   Newton iteration:
#'   \deqn{w_{n+1} = w_n - \frac{f(w_n)}{f'(w_n) - \dfrac{f(w_n) f''(w_n)}{2 f'(w_n)}}.}
#' }
#' \strong{Convergence.} Near either end of the domain, \eqn{f'(w) = e^w(w+1)}
#' becomes very small, so a small residual \eqn{|f(w)|} alone does not
#' guarantee an accurate \eqn{w}. Convergence is therefore declared when
#' either the residual or the relative step size falls below \code{tol}.
#'
#' Elements outside \eqn{[-1/e, 0)} (other than \code{NA}) return
#' \code{NaN} rather than raising an error, so a single out-of-range
#' element does not abort a whole batch of evaluations.
#'
#' @param z Numeric vector with values in \eqn{[-1/e, 0)}.
#' @param tol Tolerance for convergence.
#' @param max_iter Maximum iterations allowed.
#' @return Numeric vector of evaluated Lambert W_-1 values, the same
#'   length as \code{z}.
#' @references
#' Corless, R. M., Gonnet, G. H., Hare, D. E. G., Jeffrey, D. J., &
#' Knuth, D. E. (1996). On the Lambert W function. \emph{Advances in
#' Computational Mathematics}, 5(1), 329-359.
#' @keywords internal
.lambert_w_minus1 <- function(z, tol = 1e-12, max_iter = 100) {
  z <- as.numeric(z)
  w <- rep(NA_real_, length(z))
  valid <- !is.na(z) & z < 0 & z >= -1 / exp(1)
  w[!is.na(z) & !valid] <- NaN
  if (!any(valid)) return(w)

  zc <- z[valid]
  at_branch_point <- zc == -1 / exp(1)

  L1 <- log(-zc)
  L2 <- log(-L1)
  wc <- L1 - L2 + (L2 / L1)
  wc[at_branch_point] <- -1

  active <- !at_branch_point
  for (iter in seq_len(max_iter)) {
    if (!any(active)) break
    idx <- which(active)
    ew <- exp(wc[idx])
    f <- wc[idx] * ew - zc[idx]
    f1 <- ew * (wc[idx] + 1)
    f2 <- ew * (wc[idx] + 2)
    denom <- f1 - (f * f2) / (2 * f1)
    ok <- denom != 0 & is.finite(denom)
    delta <- numeric(length(idx))
    delta[ok] <- f[ok] / denom[ok]
    wc[idx] <- wc[idx] - delta
    converged <- (abs(f) < tol) | (abs(delta) < tol * pmax(1, abs(wc[idx])))
    active[idx[converged]] <- FALSE
  }
  if (any(active)) {
    warning(sprintf(
      "Lambert W_-1 did not converge to tolerance %.1e for %d value(s) after %d iterations; results may be inaccurate.",
      tol, sum(active), max_iter), call. = FALSE)
  }
  w[valid] <- wc
  w
}

