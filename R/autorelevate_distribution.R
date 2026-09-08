#' Probability Density Function (PDF) for the Autorelevated Family
#'
#' @description
#' Density function of the autorelevated family built from a baseline
#' distribution specified by \code{dist}, \code{p1}, and \code{p2}.
#'
#' @details
#' The autorelevation transform of Krakowski (1973) and
#' Dileepkumar and Sankaran (2022), applied to a baseline density
#' \eqn{f_0(x)} with cumulative hazard \eqn{\Lambda_0(x)}, gives the PDF
#' \deqn{f(x) = f_0(x)\, \Lambda_0(x), \qquad x > 0.}
#' No extra parameter is introduced relative to the baseline model.
#'
#' For \code{dist = "weibull"} (\eqn{\Lambda_0(x) = p1\, x^{p2}}), this
#' reduces to the Autorelevated Weibull density of Dileep Kumar,
#' Shabeer, and Sankaran (2025, Eq. 2.3),
#' \deqn{f(x) = p2\, p1^2\, x^{2\,p2 - 1}\, e^{-p1\, x^{p2}},}
#' with \code{p1} playing the role of the Weibull rate \eqn{\lambda} and
#' \code{p2} the shape \eqn{\beta}. A convenient by-product (their
#' Theorem 2.2) is that if \eqn{X} is Autorelevated Weibull with
#' parameters \code{p1}, \code{p2}, then \eqn{Y = X^{p2}} follows a Gamma
#' distribution with shape 2 and scale \code{1 / p1}; this gives an exact
#' distributional check for \code{dist = "weibull"}.
#'
#' See \link{autorelevate-package} for the cumulative hazard of all ten
#' supported baseline distributions.
#'
#' @param x Vector of quantiles.
#' @param dist Baseline distribution: \code{"weibull"},
#'   \code{"lomax"}, \code{"burr"}, \code{"gompertz"},
#'   \code{"loglogistic"}, \code{"chen"}, \code{"expexp"},
#'   \code{"powerlindley"}, \code{"lognormal"}, or \code{"gamma"}. See
#'   \link{autorelevate-package} for the role of \code{p1} and \code{p2}
#'   within each distribution.
#' @param p1 Baseline parameter 1 (interpretation depends on \code{dist}).
#' @param p2 Baseline parameter 2 (interpretation depends on \code{dist}).
#' @param log Logical; if TRUE, densities are returned as log(f).
#' @return Numeric vector of density values.
#' @references
#' Krakowski, M. (1973). The relevation transform and a generalization of
#' the gamma distribution function. \emph{Revue francaise d'automatique,
#' informatique, recherche operationnelle. Recherche operationnelle},
#' 7(V2), 107-120. \doi{10.1051/ro/197307V201071}
#'
#' Dileepkumar, M., & Sankaran, P. G. (2022). Some results of
#' auto-relevation transform in reliability analysis. \emph{Statistics and
#' Applications}, 20(2), 251-263.
#'
#' Dileep Kumar, M., Shabeer, A. M., & Sankaran, P. G. (2025). Reliability
#' properties and applications of autorelevated Weibull distribution.
#' \emph{American Journal of Mathematical and Management Sciences}, 44(3-4),
#' 215-237. \doi{10.1080/01966324.2026.2665479}
#' @family autorelevate distribution functions
#' @examples
#' x <- seq(0.1, 5, by = 0.1)
#'
#' # Weibull baseline (Autorelevated Weibull)
#' fx <- dautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5)
#' plot(x, fx, type = "l", ylab = "Density", main = "Autorelevated Weibull")
#'
#' # Log-density, useful for likelihood-based estimation
#' dautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5, log = TRUE)
#'
#' # The five newer baseline distributions
#' dautorelevate(x, dist = "chen", p1 = 0.3, p2 = 1.4)
#' dautorelevate(x, dist = "expexp", p1 = 0.8, p2 = 2.3)
#' dautorelevate(x, dist = "powerlindley", p1 = 1.2, p2 = 1.7)
#' dautorelevate(x, dist = "lognormal", p1 = 0.5, p2 = 0.8)
#' dautorelevate(x, dist = "gamma", p1 = 1.1, p2 = 2.4)
#' @export
dautorelevate <- function(x, dist = "weibull", p1 = 0.5, p2 = 1.5, log = FALSE) {
  .validate_ar_inputs(dist, p1, p2)
  pdf_val <- numeric(length(x))
  ok <- !is.na(x) & x > 0
  if (any(ok)) {
    Lambda <- .baseline_cum_hazard(x[ok], dist, p1, p2)
    f0 <- .baseline_pdf(x[ok], dist, p1, p2)
    pdf_val[ok] <- pmax(f0 * Lambda, 0)
  }
  pdf_val[is.na(x)] <- NA_real_
  if (log) return(log(pdf_val))
  pdf_val
}

#' Cumulative Distribution Function (CDF) for the Autorelevated Family
#'
#' @description
#' Distribution function of the autorelevated family built from a
#' baseline distribution specified by \code{dist}, \code{p1}, and
#' \code{p2}.
#'
#' @details
#' Using the baseline cumulative hazard \eqn{\Lambda_0(x)}, the survival
#' function of the autorelevation transform is
#' \deqn{\bar F(x) = e^{-\Lambda_0(x)}\left(1 + \Lambda_0(x)\right)}
#' (see \code{\link{sautorelevate}} for the derivation), so the CDF is
#' \deqn{F(x) = 1 - e^{-\Lambda_0(x)}\left(1 + \Lambda_0(x)\right).}
#' Equivalently, \eqn{\bar F(x) = q(\bar F_0(x))} for the concave
#' distortion function \eqn{q(t) = t(1 - \log t)} on \eqn{[0, 1]}: the
#' autorelevated family is, for every choice of \code{dist}, a distorted
#' version of its own baseline distribution (Dileepkumar & Sankaran, 2022).
#'
#' @param q Vector of quantiles.
#' @param dist Baseline distribution.
#' @param p1 Baseline parameter 1.
#' @param p2 Baseline parameter 2.
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X <= x], otherwise P[X > x].
#' @param log.p Logical; if TRUE, probabilities are returned as log(p).
#' @return Numeric vector of cumulative probability values.
#' @references
#' Dileepkumar, M., & Sankaran, P. G. (2022). Some results of
#' auto-relevation transform in reliability analysis. \emph{Statistics and
#' Applications}, 20(2), 251-263.
#' @family autorelevate distribution functions
#' @examples
#' q <- seq(0.1, 5, by = 0.5)
#' pautorelevate(q, dist = "weibull", p1 = 0.5, p2 = 1.5)
#'
#' # Upper tail P[X > x]
#' pautorelevate(q, dist = "weibull", p1 = 0.5, p2 = 1.5, lower.tail = FALSE)
#'
#' # CDF + survival function sum to 1, for every baseline distribution
#' cdf <- pautorelevate(q, dist = "gamma", p1 = 1.1, p2 = 2.4)
#' surv <- sautorelevate(q, dist = "gamma", p1 = 1.1, p2 = 2.4)
#' all.equal(cdf + surv, rep(1, length(q)))
#' @export
pautorelevate <- function(q, dist = "weibull", p1 = 0.5, p2 = 1.5,
                           lower.tail = TRUE, log.p = FALSE) {
  .validate_ar_inputs(dist, p1, p2)
  cdf_val <- numeric(length(q))
  ok <- !is.na(q) & q > 0
  if (any(ok)) {
    Lambda <- .baseline_cum_hazard(q[ok], dist, p1, p2)
    cdf_val[ok] <- 1 - exp(-Lambda) * (1 + Lambda)
  }
  cdf_val[is.na(q)] <- NA_real_
  if (!lower.tail) cdf_val <- 1 - cdf_val
  if (log.p) return(log(cdf_val))
  cdf_val
}

#' Survival Function for the Autorelevated Family
#'
#' @description
#' Survival (reliability) function of the autorelevated family built
#' from a baseline distribution specified by \code{dist}, \code{p1}, and
#' \code{p2}.
#'
#' @details
#' Let \eqn{\bar F_0(x)} and \eqn{\Lambda_0(x) = -\log \bar F_0(x)} be the
#' baseline survival function and cumulative hazard. For identically
#' distributed populations, the relevation transform of Krakowski (1973)
#' reduces to the autorelevation survival function
#' \deqn{\bar F(x) = \bar F_0(x)\left(1 - \log \bar F_0(x)\right)
#'   = e^{-\Lambda_0(x)}\left(1 + \Lambda_0(x)\right).}
#' For \code{dist = "weibull"} this is the Autorelevated Weibull survival function of
#' Dileep Kumar, Shabeer, and Sankaran (2025, Eq. 2.2),
#' \eqn{\bar F(x) = e^{-p1\,x^{p2}}\left(1 + p1\,x^{p2}\right)}.
#'
#' @param x Vector of quantiles.
#' @param dist Baseline distribution.
#' @param p1 Baseline parameter 1.
#' @param p2 Baseline parameter 2.
#' @param log Logical; if TRUE, survival values are returned as log(S).
#' @return Numeric vector of survival values.
#' @references
#' Krakowski, M. (1973). The relevation transform and a generalization of
#' the gamma distribution function. \emph{Revue francaise d'automatique,
#' informatique, recherche operationnelle. Recherche operationnelle},
#' 7(V2), 107-120. \doi{10.1051/ro/197307V201071}
#' @family autorelevate distribution functions
#' @examples
#' x <- seq(0.1, 5, by = 0.5)
#' sautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5)
#' sautorelevate(x, dist = "gompertz", p1 = 0.3, p2 = 0.5, log = TRUE)
#' @export
sautorelevate <- function(x, dist = "weibull", p1 = 0.5, p2 = 1.5, log = FALSE) {
  .validate_ar_inputs(dist, p1, p2)
  surv_val <- rep(1, length(x))
  ok <- !is.na(x) & x > 0
  if (any(ok)) {
    Lambda <- .baseline_cum_hazard(x[ok], dist, p1, p2)
    surv_val[ok] <- exp(-Lambda) * (1 + Lambda)
  }
  surv_val[is.na(x)] <- NA_real_
  if (log) return(log(surv_val))
  surv_val
}

#' Hazard Rate Function for the Autorelevated Family
#'
#' @description
#' Hazard (failure) rate function of the autorelevated family built from
#' a baseline distribution specified by \code{dist}, \code{p1}, and
#' \code{p2}.
#'
#' @details
#' The hazard rate is the ratio of the PDF to the survival function,
#' \eqn{h(x) = f(x) / \bar F(x)}. Substituting
#' \eqn{f(x) = f_0(x)\Lambda_0(x)} and
#' \eqn{\bar F(x) = e^{-\Lambda_0(x)}(1 + \Lambda_0(x))}, and using
#' \eqn{f_0(x) = h_0(x)\, e^{-\Lambda_0(x)}} for the baseline hazard
#' \eqn{h_0(x)}, gives the compact identity
#' \deqn{h(x) = h_0(x)\, \frac{\Lambda_0(x)}{1 + \Lambda_0(x)}.}
#' Since \eqn{\Lambda_0/(1+\Lambda_0) \in (0, 1)} is increasing in
#' \eqn{\Lambda_0}, the autorelevated hazard is always a *damped* version
#' of the baseline hazard, but the damping factor grows with \eqn{x},
#' which is what allows a monotone baseline hazard to become non-monotone
#' after autorelevation. Dileep Kumar, Shabeer, and Sankaran (2025,
#' Theorem 4.1) show that for \code{dist = "weibull"} with shape
#' \code{p2 = }\eqn{\beta}: the hazard is increasing (IHR) for
#' \eqn{\beta > 1}, upside-down bathtub (UBT) for \eqn{1/2 < \beta < 1}
#' (a strict, open interval -- at \eqn{\beta = 1} exactly, direct
#' calculus on \eqn{h(x) = p1^2 x/(1 + p1 x)} shows
#' \eqn{h'(x) = p1^2/(1+p1x)^2 > 0} for all \eqn{x}, i.e. strictly
#' increasing with no decreasing phase, so \eqn{\beta = 1} belongs with
#' the IHR case, not UBT), and decreasing (DHR) for \eqn{0 < \beta \le 1/2}.
#'
#' @param x Vector of quantiles.
#' @param dist Baseline distribution.
#' @param p1 Baseline parameter 1.
#' @param p2 Baseline parameter 2.
#' @return Numeric vector of hazard rate values.
#' @references
#' Dileep Kumar, M., Shabeer, A. M., & Sankaran, P. G. (2025). Reliability
#' properties and applications of autorelevated Weibull distribution.
#' \emph{American Journal of Mathematical and Management Sciences}, 44(3-4),
#' 215-237. \doi{10.1080/01966324.2026.2665479}
#' @family autorelevate distribution functions
#' @examples
#' x <- seq(0.1, 5, by = 0.1)
#'
#' # Upside-down bathtub hazard (1/2 < beta <= 1)
#' h_ubt <- haautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 0.8)
#' plot(x, h_ubt, type = "l", ylab = "Hazard rate", main = "UBT-shaped hazard")
#'
#' # Increasing hazard (beta > 1)
#' h_ihr <- haautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5)
#' lines(x, h_ihr, col = "blue")
#' @export
haautorelevate <- function(x, dist = "weibull", p1 = 0.5, p2 = 1.5) {
  .validate_ar_inputs(dist, p1, p2)
  f_ar <- dautorelevate(x, dist, p1, p2)
  s_ar <- sautorelevate(x, dist, p1, p2)
  ifelse(is.na(x), NA_real_,
         ifelse(x <= 0, 0, ifelse(s_ar > 0, f_ar / s_ar, 0)))
}

#' Quantile Function via Exact Lambert W_-1 Inversion
#'
#' @description
#' Quantile function of the autorelevated family, obtained analytically
#' via the negative branch of the Lambert W function rather than
#' numerical root-finding.
#'
#' @details
#' We want \eqn{x_p} solving \eqn{F(x_p) = p}, i.e.,
#' \eqn{\bar F(x_p) = 1 - p}. Writing \eqn{S = 1 - p} and
#' \eqn{y = 1 + \Lambda_0(x_p)}, the survival identity
#' \eqn{\bar F(x) = e^{-\Lambda_0(x)}(1 + \Lambda_0(x))} becomes
#' \eqn{y\, e^{-y} = S/e}. Since \eqn{W_{-1}(z)} is defined by
#' \eqn{W_{-1}(z)\, e^{W_{-1}(z)} = z} and \eqn{y \ge 1}, setting
#' \eqn{w = -y} gives \eqn{w e^w = -S/e}, so
#' \deqn{w = W_{-1}(-S/e), \qquad \Lambda_0(x_p) = y - 1 = -1 - w.}
#' Applying the baseline inverse cumulative hazard then gives the exact
#' quantile \eqn{x_p = \Lambda_0^{-1}(-1 - W_{-1}(-S/e))}. For most
#' baselines \eqn{\Lambda_0^{-1}} is elementary algebra; for
#' \code{"powerlindley"} it is itself obtained via a second application
#' of \eqn{W_{-1}} (see \link{baseline_internal}), so no numerical
#' root-finding is used anywhere in this function for any baseline.
#' This generalizes Theorem 2.3 of Dileep Kumar, Shabeer, and Sankaran
#' (2025), which gives the corresponding identity for the Weibull
#' baseline.
#'
#' @param p Vector of probabilities.
#' @param dist Baseline distribution.
#' @param p1 Baseline parameter 1.
#' @param p2 Baseline parameter 2.
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X <= x].
#' @return Numeric vector of evaluated quantiles.
#' @references
#' Dileep Kumar, M., Shabeer, A. M., & Sankaran, P. G. (2025). Reliability
#' properties and applications of autorelevated Weibull distribution.
#' \emph{American Journal of Mathematical and Management Sciences}, 44(3-4),
#' 215-237. \doi{10.1080/01966324.2026.2665479}
#'
#' Ghitany, M. E., Al-Mutairi, D. K., Balakrishnan, N., & Al-Enezi, L. J.
#' (2013). Power Lindley distribution and associated inference.
#' \emph{Computational Statistics & Data Analysis}, 64, 20-33.
#' @family autorelevate distribution functions
#' @examples
#' p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' qautorelevate(p, dist = "weibull", p1 = 0.5, p2 = 1.5)
#'
#' # Round trip: CDF then quantile recovers the original x, for every distribution
#' x <- seq(0.5, 3, by = 0.5)
#' cdf <- pautorelevate(x, dist = "powerlindley", p1 = 1.2, p2 = 1.7)
#' qautorelevate(cdf, dist = "powerlindley", p1 = 1.2, p2 = 1.7)
#' @export
qautorelevate <- function(p, dist = "weibull", p1 = 0.5, p2 = 1.5, lower.tail = TRUE) {
  .validate_ar_inputs(dist, p1, p2)
  if (!lower.tail) p <- 1 - p
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  invalid_p <- ok & (p < 0 | p > 1)
  if (any(invalid_p)) {
    warning("NaNs produced (probabilities outside [0, 1])", call. = FALSE)
    out[invalid_p] <- NaN
    ok[invalid_p] <- FALSE
  }
  out[ok & p <= 0] <- 0
  out[ok & p >= 1] <- Inf
  mid <- ok & p > 0 & p < 1
  if (any(mid)) {
    S <- 1 - p[mid]
    z <- -S / exp(1)
    w <- .lambert_w_minus1(z)
    Lambda_target <- -1 - w
    out[mid] <- .inv_baseline_cum_hazard(Lambda_target, dist, p1, p2)
  }
  out
}

#' Random Generation for the Autorelevated Family
#'
#' @description
#' Draws random variates from the autorelevated family via the inverse
#' transform method, using the exact quantile function
#' \code{\link{qautorelevate}}.
#'
#' @details
#' Because \code{\link{qautorelevate}} is available in closed form (no
#' numerical inversion), generation is exact and fast: a single uniform
#' draw \eqn{U \sim \mathrm{Unif}(0,1)} per variate is transformed via
#' \eqn{X = Q(U)}, where \eqn{Q} is the autorelevated quantile function.
#'
#' @param n Number of observations.
#' @param dist Baseline distribution.
#' @param p1 Baseline parameter 1.
#' @param p2 Baseline parameter 2.
#' @return Numeric vector of simulated random variates.
#' @importFrom stats runif
#' @family autorelevate distribution functions
#' @examples
#' set.seed(1)
#' x <- rautorelevate(500, dist = "weibull", p1 = 0.5, p2 = 1.5)
#' hist(x, breaks = 30, probability = TRUE, main = "Simulated Autorelevated Weibull sample")
#' curve(dautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5),
#'       add = TRUE, col = "blue", lwd = 2)
#' @export
rautorelevate <- function(n, dist = "weibull", p1 = 0.5, p2 = 1.5) {
  .validate_ar_inputs(dist, p1, p2)
  u <- stats::runif(n)
  qautorelevate(u, dist = dist, p1 = p1, p2 = p2)
}

