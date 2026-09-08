#' Supported Baseline Distribution Families
#'
#' @description
#' The ten baseline distribution names accepted by the \code{dist} argument
#' throughout this package.
#' @keywords internal
.AR_VALID_DISTS <- c("weibull", "lomax", "burr", "gompertz", "loglogistic",
                      "chen", "expexp", "powerlindley", "lognormal", "gamma")

#' Validate Inputs Shared Across Exported Functions
#'
#' @description
#' Checks \code{dist}, \code{p1}, and \code{p2} (and optionally
#' \code{data}) for validity, and raises an informative error rather than
#' letting invalid input silently propagate as \code{NaN}/\code{Inf}.
#'
#' @details
#' Every exported function in this package calls this validator first.
#' This keeps error messages consistent and means a typo in \code{dist}
#' or a negative parameter is reported immediately, with the offending
#' value and the list of valid options, rather than surfacing later as a
#' cryptic numerical failure inside \code{\link[stats]{optim}} or the
#' Lambert W solver.
#'
#' @param dist Baseline distribution.
#' @param p1 Baseline parameter 1.
#' @param p2 Baseline parameter 2.
#' @param data Optional numeric vector to validate as well.
#' @keywords internal
.validate_ar_inputs <- function(dist, p1, p2, data = NULL) {
  if (!is.character(dist) || length(dist) != 1 || !(dist %in% .AR_VALID_DISTS)) {
    stop("`dist` must be one of: ", paste(shQuote(.AR_VALID_DISTS), collapse = ", "),
         ". Got: ", if (is.character(dist)) shQuote(dist) else class(dist)[1], call. = FALSE)
  }
  if (!is.numeric(p1) || length(p1) != 1 || !is.finite(p1) || p1 <= 0) {
    stop("`p1` must be a single positive, finite number. Got: ", format(p1), call. = FALSE)
  }
  if (!is.numeric(p2) || length(p2) != 1 || !is.finite(p2) || p2 <= 0) {
    stop("`p2` must be a single positive, finite number. Got: ", format(p2), call. = FALSE)
  }
  if (!is.null(data)) {
    if (!is.numeric(data) || length(data) == 0) {
      stop("`data` must be a non-empty numeric vector.", call. = FALSE)
    }
    if (!any(is.finite(data) & data > 0)) {
      stop("`data` must contain at least one strictly positive, finite value ",
           "(the autorelevated family has support on (0, Inf)).", call. = FALSE)
    }
  }
  invisible(TRUE)
}

