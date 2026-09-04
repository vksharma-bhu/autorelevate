#' Baseline Functions for Cumulative Hazard and Density
#'
#' @description
#' Internal helpers giving the cumulative hazard \eqn{\Lambda_0(x)}, its
#' inverse, and the density \eqn{f_0(x)} of the ten baseline families
#' supported by this package. See \link{autorelevate-package} for the
#' definition of each family in terms of \code{p1} and \code{p2}.
#'
#' @details
#' Every autorelevated quantity in this package (\code{\link{dautorelevate}},
#' \code{\link{pautorelevate}}, \code{\link{sautorelevate}},
#' \code{\link{haautorelevate}}, \code{\link{qautorelevate}}) is built
#' purely from \eqn{\Lambda_0}, its inverse, and \eqn{f_0}; no other
#' baseline-specific code is required. This keeps all ten families
#' interchangeable and makes it straightforward to add further baselines
#' by extending the three \code{switch} blocks below with a new
#' cumulative hazard, inverse cumulative hazard, and density.
#'
#' The baseline density is recovered from the cumulative hazard via
#' \eqn{f_0(x) = \Lambda_0'(x)\, e^{-\Lambda_0(x)}}, so \eqn{\Lambda_0}
#' alone determines the whole baseline model. For \code{"lognormal"} and
#' \code{"gamma"}, \eqn{\Lambda_0} and its inverse are computed on the log
#' scale via \code{\link[stats]{pnorm}}/\code{\link[stats]{qnorm}} and
#' \code{\link[stats]{pgamma}}/\code{\link[stats]{qgamma}} with
#' \code{log.p = TRUE}, which is both exact and numerically stable in the
#' tails (avoiding cancellation from computing \code{1 - p} for small
#' \code{p}). For \code{"powerlindley"}, the inverse cumulative hazard has
#' no elementary closed form and is instead obtained via the same
#' Lambert \eqn{W_{-1}} solver (\code{\link{.lambert_w_minus1}}) used for
#' \code{\link{qautorelevate}} itself; see Ghitany, Al-Mutairi,
#' Balakrishnan, and Al-Enezi (2013) for the underlying Lindley-family
#' identity that makes this possible.
#'
#' @param x,Lambda Numeric vectors.
#' @param dist Baseline distribution family.
#' @param p1 Parameter 1.
#' @param p2 Parameter 2.
#' @keywords internal
#' @name baseline_internal
NULL

#' @rdname baseline_internal
.baseline_cum_hazard <- function(x, dist, p1, p2) {
  switch(dist,
    weibull      = p1 * (x ^ p2),
    lomax        = p2 * log(1 + p1 * x),
    burr         = p2 * log(1 + x ^ p1),
    gompertz     = p1 * (exp(p2 * x) - 1),
    loglogistic  = log(1 + (x / p1) ^ p2),
    chen         = p1 * (exp(x ^ p2) - 1),
    expexp       = -log(1 - (-expm1(-p1 * x)) ^ p2),
    powerlindley = p1 * x ^ p2 - log(1 + (p1 * x ^ p2) / (p1 + 1)),
    lognormal    = -stats::pnorm((log(x) - p1) / p2, lower.tail = FALSE, log.p = TRUE),
    gamma        = -stats::pgamma(x, shape = p2, rate = p1, lower.tail = FALSE, log.p = TRUE),
    stop("Unsupported distribution family: ", dist)
  )
}

#' @rdname baseline_internal
.inv_baseline_cum_hazard <- function(Lambda, dist, p1, p2) {
  switch(dist,
    weibull      = (Lambda / p1) ^ (1 / p2),
    lomax        = (exp(Lambda / p2) - 1) / p1,
    burr         = (exp(Lambda / p2) - 1) ^ (1 / p1),
    gompertz     = (1 / p2) * log(1 + Lambda / p1),
    loglogistic  = p1 * (exp(Lambda) - 1) ^ (1 / p2),
    chen         = (log(1 + Lambda / p1)) ^ (1 / p2),
    expexp       = -(1 / p1) * log(1 - (-expm1(-Lambda)) ^ (1 / p2)),
    powerlindley = {
      q <- exp(-Lambda)
      z <- -q * (p1 + 1) * exp(-(p1 + 1))
      w <- .lambert_w_minus1(z)
      ((-w - (p1 + 1)) / p1) ^ (1 / p2)
    },
    lognormal    = exp(p1 + p2 * stats::qnorm(-Lambda, lower.tail = FALSE, log.p = TRUE)),
    gamma        = stats::qgamma(-Lambda, shape = p2, rate = p1, lower.tail = FALSE, log.p = TRUE),
    stop("Unsupported distribution family: ", dist)
  )
}

#' @rdname baseline_internal
.baseline_pdf <- function(x, dist, p1, p2) {
  Lambda <- .baseline_cum_hazard(x, dist, p1, p2)
  switch(dist,
    weibull  = p2 * p1 * (x ^ (p2 - 1)) * exp(-Lambda),
    lomax    = p2 * p1 * ((1 + p1 * x) ^ (-p2 - 1)),
    burr     = p1 * p2 * (x ^ (p1 - 1)) * ((1 + x ^ p1) ^ (-p2 - 1)),
    gompertz = p1 * p2 * exp(p2 * x) * exp(-Lambda),
    loglogistic = {
      xa <- (x / p1) ^ p2
      (p2 / p1) * (x / p1) ^ (p2 - 1) / ((1 + xa) ^ 2)
    },
    chen         = p1 * p2 * (x ^ (p2 - 1)) * exp(x ^ p2) * exp(-Lambda),
    expexp       = p1 * p2 * exp(-p1 * x) * ((-expm1(-p1 * x)) ^ (p2 - 1)),
    powerlindley = (p1 ^ 2 / (p1 + 1)) * p2 * (x ^ (p2 - 1)) * (1 + x ^ p2) * exp(-p1 * x ^ p2),
    lognormal    = stats::dnorm((log(x) - p1) / p2) / (x * p2),
    gamma        = stats::dgamma(x, shape = p2, rate = p1),
    stop("Unsupported distribution family: ", dist)
  )
}

