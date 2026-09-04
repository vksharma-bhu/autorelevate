#' Data-Adaptive Starting Values for MLE
#'
#' @description
#' Finds a starting point for the MLE optimizer by evaluating the
#' log-likelihood on a coarse grid spanning many orders of magnitude of
#' \code{p1} and \code{p2}, and returning the (log-scale) grid point with
#' the best finite log-likelihood.
#'
#' @details
#' A single fixed starting point (e.g. \code{p1 = 0.5}, \code{p2 = 1})
#' is not safe across arbitrary data scales: for several families
#' (\code{"gompertz"}, \code{"chen"}, \code{"expexp"}) the density
#' involves \code{exp()} of a term that grows with \code{x}, so a
#' starting \code{p2} that is entirely reasonable for data on the order
#' of 1 can silently overflow (or the survival term underflow to exactly
#' 0) for data with much larger values -- for example, \code{p2 = 1} on
#' data containing a value of 79 makes \code{exp(p2*x)} overflow before
#' optimization even begins, leaving the optimizer with no usable
#' gradient anywhere nearby. This grid search guarantees the optimizer
#' starts from a point with a finite, evaluable log-likelihood for
#' whatever scale the data happen to be on, for every supported baseline
#' family, at the cost of a bounded number of cheap, fully-vectorized
#' density evaluations (14x14 = 196 grid points).
#'
#' @param data Vector of positive sample observations.
#' @param dist Baseline distribution family.
#' @return Numeric vector of length 2: \code{c(log(p1), log(p2))}.
#' @keywords internal
.mle_seed_grid <- function(data, dist) {
  grid_vals <- c(1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.5, 1, 2, 5, 10, 50, 100, 1000)
  best_val <- Inf
  best_par <- c(log(0.5), log(1))
  for (p1 in grid_vals) {
    for (p2 in grid_vals) {
      ld <- tryCatch(dautorelevate(data, dist = dist, p1 = p1, p2 = p2, log = TRUE),
                      error = function(e) rep(NA_real_, length(data)))
      if (all(is.finite(ld))) {
        val <- -sum(ld)
        if (is.finite(val) && val < best_val) {
          best_val <- val
          best_par <- c(log(p1), log(p2))
        }
      }
    }
  }
  best_par
}

#' Robust Optimization Helper
#'
#' @description
#' Wraps \code{\link[stats]{optim}} with a fallback chain: BFGS first,
#' then Nelder-Mead if BFGS fails or does not converge, then several
#' randomly perturbed Nelder-Mead restarts if both fail. Never lets a raw
#' optimizer error propagate to the user; instead raises a single clear
#' error only if every attempt fails.
#'
#' @param obj_func Objective function to minimize, taking a length-2
#'   numeric vector.
#' @param init_pars Initial parameter vector (length 2).
#' @param label Short label used in the error message if all attempts fail.
#' @return The best \code{optim()} result found (a list with at least
#'   \code{par}, \code{value}, \code{convergence}).
#' @keywords internal
.robust_optim <- function(obj_func, init_pars, label = "model") {
  attempt <- function(par0, method) {
    tryCatch(stats::optim(par0, obj_func, method = method),
             error = function(e) NULL, warning = function(w) NULL)
  }
  best <- attempt(init_pars, "BFGS")
  if (is.null(best) || !is.finite(best$value) || best$convergence != 0) {
    alt <- attempt(init_pars, "Nelder-Mead")
    if (!is.null(alt) && is.finite(alt$value) &&
        (is.null(best) || !is.finite(best$value) || alt$value < best$value)) {
      best <- alt
    }
  }
  if (is.null(best) || !is.finite(best$value)) {
    for (k in 1:5) {
      par_try <- init_pars + stats::rnorm(length(init_pars), sd = 0.5 * k)
      alt <- attempt(par_try, "Nelder-Mead")
      if (!is.null(alt) && is.finite(alt$value) &&
          (is.null(best) || alt$value < best$value)) {
        best <- alt
      }
    }
  }
  if (is.null(best) || !is.finite(best$value)) {
    stop("Optimization failed to converge while fitting ", label,
         ". Check that `data` are positive, finite, and plausible for ",
         "the chosen `dist`.", call. = FALSE)
  }
  best
}

#' Fit Autorelevated Family Parameters
#'
#' @description
#' Fits the parameters \code{p1}, \code{p2} of an autorelevated
#' distribution to data by Maximum Likelihood (MLE), Maximum Product of
#' Spacings (MPS), Least Squares (LS), Weighted Least Squares (WLS), or
#' Cramer-von Mises (CvM) minimum-distance estimation.
#'
#' @details
#' Let \eqn{x_{(1)} \le \dots \le x_{(n)}} be the ordered, strictly
#' positive observations. Optimization is over \eqn{(\log p1, \log p2)}
#' (via \code{\link[stats]{optim}}, with an automatic Nelder-Mead/restart
#' fallback if the primary optimizer fails to converge; see
#' \code{\link{.robust_optim}}), which enforces positivity without
#' constrained optimization. The starting point itself is chosen
#' adaptively from the data via a coarse log-scale grid search (see
#' \code{\link{.mle_seed_grid}}), rather than a single fixed default,
#' since a fixed starting point can silently overflow or underflow for
#' several families when the data's scale differs greatly from 1.
#' \itemize{
#'   \item \strong{MLE} maximizes the log-likelihood
#'   \eqn{\sum_i \log f(x_{(i)})}, following the ARW log-likelihood
#'   derivation of Dileep Kumar, Shabeer, and Sankaran (2025, Sec. 6.1),
#'   generalized to any baseline via \code{\link{dautorelevate}}.
#'   Standard errors are obtained from the observed Fisher information
#'   (the numerical Hessian at the optimum, evaluated on the
#'   \eqn{(\log p1, \log p2)} scale used internally by the optimizer and
#'   mapped back to the \eqn{(p1, p2)} scale via the delta method,
#'   \eqn{\widehat{\mathrm{Var}}(p) \approx D\, \widehat{\mathrm{Var}}(\log p)\, D}
#'   with \eqn{D = \mathrm{diag}(p1, p2)}).
#'   \item \strong{LS and WLS} minimize, respectively,
#'   \deqn{\sum_{i=1}^n \left(F(x_{(i)}) - \frac{i}{n+1}\right)^2 \quad\text{and}\quad
#'         \sum_{i=1}^n \frac{(n+1)^2(n+2)}{i(n+1-i)}\left(F(x_{(i)}) - \frac{i}{n+1}\right)^2,}
#'   following Swain, Venkatraman, and Wilson (1988).
#'   \item \strong{CvM} minimizes the Cramer-von Mises statistic of Choi
#'   and Bulgren (1968),
#'   \deqn{\frac{1}{12n} + \sum_{i=1}^n \left(F(x_{(i)}) - \frac{2i-1}{2n}\right)^2.}
#'   \item \strong{MPS} maximizes the geometric mean of the spacings
#'   \eqn{D_i = F(x_{(i)}) - F(x_{(i-1)})} (with \eqn{F(x_{(0)}) = 0},
#'   \eqn{F(x_{(n+1)}) = 1}), i.e., minimizes \eqn{-\sum_i \log D_i}.
#' }
#' For MPS/LS/WLS/CvM, the optimizer is seeded at the MLE solution.
#'
#' \strong{Model-selection criteria} (\code{method = "mle"} only, since
#' all four require the maximized log-likelihood \eqn{\hat\ell} and are
#' only asymptotically justified through it), with \eqn{k = 2}
#' parameters for every family in this package:
#' \deqn{\mathrm{AIC} = 2k - 2\hat\ell, \quad
#'       \mathrm{BIC} = k \log n - 2\hat\ell, \quad
#'       \mathrm{CAIC} = \frac{2kn}{n-k-1} - 2\hat\ell, \quad
#'       \mathrm{HQIC} = 2k \log(\log n) - 2\hat\ell,}
#' matching the criteria reported by Dileep Kumar, Shabeer, and Sankaran
#' (2025, Sec. 8.1). \code{CAIC} is \code{NA} when \eqn{n \le k+1} (too
#' few observations for the correction term to be defined); note that
#' this particular \code{CAIC} formula is algebraically identical to the
#' corrected AIC (AICc) of Hurvich and Tsai (1989) — some sources use
#' "CAIC" for a different formula (Bozdogan, 1987), so compare against
#' the exact formula above, not just the acronym, when cross-referencing
#' other software.
#'
#' A Kolmogorov-Smirnov goodness-of-fit test (\code{\link[stats]{ks.test}})
#' comparing the data to the fitted CDF is returned for every method.
#'
#' @param data Vector of positive sample observations. Missing, non-finite,
#'   or non-positive values are dropped with a warning before fitting.
#' @param dist Baseline distribution family: \code{"weibull"},
#'   \code{"lomax"}, \code{"burr"}, \code{"gompertz"},
#'   \code{"loglogistic"}, \code{"chen"}, \code{"expexp"},
#'   \code{"powerlindley"}, \code{"lognormal"}, or \code{"gamma"}.
#' @param method Estimation routine: \code{"mle"}, \code{"mps"},
#'   \code{"ls"}, \code{"wls"}, or \code{"cvm"}.
#' @return Object of class \code{autorelevate_fit}: a list with
#'   components \code{par} (named numeric vector \code{p1}, \code{p2}),
#'   \code{vcov}, \code{value}, \code{convergence}, \code{data},
#'   \code{dist}, \code{method}, \code{n}, \code{log_lik}, \code{aic},
#'   \code{bic}, \code{caic}, \code{hqic}, \code{ks_stat}, and
#'   \code{ks_pval}.
#' @importFrom stats optim ks.test
#' @references
#' Choi, K., & Bulgren, W. (1968). An estimation procedure for mixtures
#' of distributions. \emph{Journal of the Royal Statistical Society
#' Series B}, 30(3), 444-460.
#'
#' Swain, J. J., Venkatraman, S., & Wilson, J. R. (1988). Least-squares
#' estimation of distribution functions in Johnson's translation system.
#' \emph{Journal of Statistical Computation and Simulation}, 29(4),
#' 271-297.
#'
#' Hurvich, C. M., & Tsai, C.-L. (1989). Regression and time series model
#' selection in small samples. \emph{Biometrika}, 76(2), 297-307.
#'
#' Dileep Kumar, M., Shabeer, A. M., & Sankaran, P. G. (2025). Reliability
#' properties and applications of autorelevated Weibull distribution.
#' \emph{American Journal of Mathematical and Management Sciences}, 44(3-4),
#' 215-237. \doi{10.1080/01966324.2026.2665479}
#' @family autorelevate estimation functions
#' @seealso \code{\link{fit_all_methods}}, \code{\link{compare_families}}
#' @examples
#' set.seed(1)
#' x <- rautorelevate(100, dist = "weibull", p1 = 0.5, p2 = 1.5)
#'
#' fit_mle <- fit_autorelevate(x, dist = "weibull", method = "mle")
#' print(fit_mle)
#' summary(fit_mle)
#'
#' fit_ls <- fit_autorelevate(x, dist = "weibull", method = "ls")
#' fit_ls$par
#'
#' # A real dataset
#' data(bladder_cancer)
#' fit_bc <- fit_autorelevate(bladder_cancer, dist = "weibull", method = "mle")
#' summary(fit_bc)
#' @export
fit_autorelevate <- function(data, dist = "weibull", method = "mle") {
  .validate_ar_inputs(dist, 1, 1, data = data)
  valid_methods <- c("mle", "mps", "ls", "wls", "cvm")
  if (!(method %in% valid_methods)) {
    stop("`method` must be one of: ", paste(shQuote(valid_methods), collapse = ", "),
         ". Got: ", shQuote(method), call. = FALSE)
  }

  n_input <- length(data)
  keep <- !is.na(data) & is.finite(data) & data > 0
  n_dropped <- n_input - sum(keep)
  if (n_dropped > 0) {
    warning(sprintf(
      "Dropped %d of %d observation(s) that were missing, non-finite, or not strictly positive; %d observation(s) used.",
      n_dropped, n_input, sum(keep)), call. = FALSE)
  }
  data <- sort(data[keep])
  n <- length(data)
  k <- 2
  i_grid <- seq_len(n)

  obj_func <- function(par_log) {
    p1 <- exp(par_log[1])
    p2 <- exp(par_log[2])
    tryCatch({
      if (method == "mle") {
        log_dens <- dautorelevate(data, dist = dist, p1 = p1, p2 = p2, log = TRUE)
        if (any(!is.finite(log_dens))) return(1e10)
        -sum(log_dens)
      } else if (method == "mps") {
        cdf_vals <- pautorelevate(data, dist = dist, p1 = p1, p2 = p2)
        D <- diff(c(0, cdf_vals, 1))
        if (any(D <= 0 | !is.finite(D))) return(1e10)
        -sum(log(D))
      } else if (method == "ls") {
        cdf_vals <- pautorelevate(data, dist = dist, p1 = p1, p2 = p2)
        sum((cdf_vals - i_grid / (n + 1)) ^ 2)
      } else if (method == "wls") {
        cdf_vals <- pautorelevate(data, dist = dist, p1 = p1, p2 = p2)
        w <- ((n + 1) ^ 2 * (n + 2)) / (i_grid * (n + 1 - i_grid))
        sum(w * (cdf_vals - i_grid / (n + 1)) ^ 2)
      } else {
        cdf_vals <- pautorelevate(data, dist = dist, p1 = p1, p2 = p2)
        1 / (12 * n) + sum((cdf_vals - (2 * i_grid - 1) / (2 * n)) ^ 2)
      }
    }, error = function(e) 1e10)
  }

  obj_mle <- function(par_log) {
    p1 <- exp(par_log[1]); p2 <- exp(par_log[2])
    tryCatch({
      log_dens <- dautorelevate(data, dist = dist, p1 = p1, p2 = p2, log = TRUE)
      if (any(!is.finite(log_dens))) return(1e10)
      -sum(log_dens)
    }, error = function(e) 1e10)
  }

  init_pars <- .mle_seed_grid(data, dist)
  mle_seed <- .robust_optim(obj_mle, init_pars, label = paste0(dist, " (MLE seed)"))
  init_pars <- mle_seed$par

  if (method == "mle") {
    opt <- mle_seed
    par_hat <- exp(opt$par)  # p1, p2 on the original (untransformed) scale
    vcov_mat <- tryCatch({
      # optimHess() is evaluated at (log p1, log p2), since obj_func is
      # parametrized that way (to enforce positivity without constrained
      # optimization). The resulting Hessian -- and its inverse -- are
      # therefore for log(p1), log(p2), NOT for p1, p2 directly. We apply
      # the delta method (Jacobian D = diag(p1, p2), since d(exp(u))/du =
      # exp(u) = the parameter itself) to transform back to the original
      # scale: Var(phi) ~= D %*% Var(theta) %*% D. Verified numerically
      # against a Hessian computed directly on the original scale.
      h <- stats::optimHess(opt$par, obj_func)
      vcov_log <- solve(h)
      D <- diag(par_hat, nrow = 2)
      D %*% vcov_log %*% D
    }, error = function(e) matrix(NA_real_, 2, 2))
    log_lik <- -opt$value
    aic <- 2 * k - 2 * log_lik
    bic <- k * log(n) - 2 * log_lik
    caic <- if (n > k + 1) (2 * k * n) / (n - k - 1) - 2 * log_lik else NA_real_
    hqic <- 2 * k * log(log(n)) - 2 * log_lik
  } else {
    opt <- .robust_optim(obj_func, init_pars, label = paste0(dist, " (", method, ")"))
    vcov_mat <- matrix(NA_real_, 2, 2)
    log_lik <- NA_real_
    aic <- NA_real_
    bic <- NA_real_
    caic <- NA_real_
    hqic <- NA_real_
  }

  estimates <- exp(opt$par)
  names(estimates) <- c("p1", "p2")

  ks_res <- suppressWarnings(
    stats::ks.test(data, pautorelevate, dist = dist, p1 = estimates["p1"], p2 = estimates["p2"])
  )

  res <- list(
    par = estimates, vcov = vcov_mat, value = opt$value,
    convergence = opt$convergence, data = data, dist = dist,
    method = method, n = n, log_lik = log_lik, aic = aic, bic = bic,
    caic = caic, hqic = hqic,
    ks_stat = as.numeric(ks_res$statistic), ks_pval = as.numeric(ks_res$p.value)
  )
  class(res) <- "autorelevate_fit"
  res
}

#' Compare All Estimation Methods for an Autorelevated Model
#'
#' @description
#' Fits a given baseline distribution using all five methods (MLE, MPS,
#' LS, WLS, CvM) and returns a comparison matrix to evaluate parameter
#' stability and Kolmogorov-Smirnov goodness-of-fit across methods.
#'
#' @details
#' Each method targets a different discrepancy between the fitted and
#' empirical distribution (likelihood, spacing products, or a weighted
#' distance between the fitted CDF and the empirical CDF; see
#' \code{\link{fit_autorelevate}} for the exact objective of each
#' method). A method that fails to converge for the data at hand
#' contributes a row of \code{NA} (via \code{\link{.robust_optim}}'s
#' error, caught here) rather than stopping the comparison for the
#' other methods.
#'
#' @param data Vector of positive sample observations.
#' @param dist Baseline distribution family.
#' @return A matrix comparing parameter estimates and KS statistics
#'   across methods (rows: \code{MLE}, \code{MPS}, \code{LS}, \code{WLS},
#'   \code{CVM}; columns: \code{Estimate_p1}, \code{Estimate_p2},
#'   \code{KS_Statistic}).
#' @family autorelevate estimation functions
#' @seealso \code{\link{fit_autorelevate}}, \code{\link{compare_families}}
#' @examples
#' set.seed(1)
#' x <- rautorelevate(100, dist = "weibull", p1 = 0.5, p2 = 1.5)
#' fit_all_methods(x, dist = "weibull")
#' @export
fit_all_methods <- function(data, dist = "weibull") {
  methods <- c("mle", "mps", "ls", "wls", "cvm")
  results <- matrix(NA_real_, nrow = 5, ncol = 3)
  rownames(results) <- toupper(methods)
  colnames(results) <- c("Estimate_p1", "Estimate_p2", "KS_Statistic")

  for (i in seq_along(methods)) {
    fit <- try(fit_autorelevate(data, dist = dist, method = methods[i]), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      results[i, 1] <- fit$par["p1"]
      results[i, 2] <- fit$par["p2"]
      results[i, 3] <- fit$ks_stat
    }
  }
  return(results)
}

#' Compare All Baseline Distribution Families for the Autorelevated Family
#'
#' @description
#' Fits a dataset across all ten supported baseline distribution
#' families using Maximum Likelihood Estimation (MLE) and returns a
#' ranked comparison table based on AIC, BIC, CAIC, HQIC, and
#' goodness-of-fit statistics.
#'
#' @details
#' Since every autorelevated family has exactly two free parameters
#' regardless of \code{dist} (the transform adds none), differences in
#' AIC/BIC/CAIC/HQIC across rows of the returned table reflect purely
#' the baseline family's shape flexibility for the data at hand, with no
#' penalty-term confound from differing parameter counts. Ranking by AIC
#' follows the model-selection approach used for the Autorelevated
#' Weibull member by Dileep Kumar, Shabeer, and Sankaran (2025).
#'
#' @param data Vector of positive sample observations.
#' @param method Estimation method, defaults to \code{"mle"} (required
#'   for AIC/BIC/CAIC/HQIC; see \code{\link{fit_autorelevate}}).
#' @return An object of class \code{autorelevate_compare}: a data frame
#'   with one row per baseline family, ranked by ascending AIC.
#' @references
#' Dileep Kumar, M., Shabeer, A. M., & Sankaran, P. G. (2025). Reliability
#' properties and applications of autorelevated Weibull distribution.
#' \emph{American Journal of Mathematical and Management Sciences}, 44(3-4),
#' 215-237. \doi{10.1080/01966324.2026.2665479}
#' @family autorelevate estimation functions
#' @seealso \code{\link{fit_autorelevate}}, \code{\link{fit_all_methods}}
#' @examples
#' set.seed(1)
#' x <- rautorelevate(150, dist = "weibull", p1 = 0.5, p2 = 1.5)
#' comp <- compare_families(x)
#' print(comp)
#'
#' # On the bundled real dataset
#' data(bladder_cancer)
#' compare_families(bladder_cancer)
#' @export
compare_families <- function(data, method = "mle") {
  results_list <- list()

  for (dist in .AR_VALID_DISTS) {
    fit <- try(fit_autorelevate(data, dist = dist, method = method), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      results_list[[dist]] <- data.frame(
        Family = dist,
        p1 = fit$par["p1"],
        p2 = fit$par["p2"],
        LogLik = fit$log_lik,
        AIC = fit$aic,
        BIC = fit$bic,
        CAIC = fit$caic,
        HQIC = fit$hqic,
        KS_Stat = fit$ks_stat,
        KS_Pval = fit$ks_pval,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results_list) == 0) {
    stop("All model fits failed for the provided data.", call. = FALSE)
  }

  res_df <- do.call(rbind, results_list)
  rownames(res_df) <- NULL
  res_df <- res_df[order(res_df$AIC), ]

  class(res_df) <- c("autorelevate_compare", "data.frame")
  res_df
}

#' Print Method for autorelevate_compare Objects
#' @param x An autorelevate_compare object.
#' @param ... Additional arguments (unused).
#' @return Invisibly, x.
#' @examples
#' set.seed(1)
#' x <- rautorelevate(100, dist = "weibull", p1 = 0.5, p2 = 1.5)
#' print(compare_families(x))
#' @export
print.autorelevate_compare <- function(x, ...) {
  cat("\nAutorelevated Family Comparison Table (Ranked by AIC)\n")
  cat("----------------------------------------------------\n")
  print.data.frame(x, digits = 4, row.names = FALSE)
  invisible(x)
}

#' Print Summary of an Autorelevated Fit
#' @param x An autorelevate_fit object.
#' @param ... Unused additional arguments.
#' @return Invisibly, x.
#' @examples
#' set.seed(1)
#' x <- rautorelevate(100, dist = "weibull", p1 = 0.5, p2 = 1.5)
#' print(fit_autorelevate(x, dist = "weibull", method = "mle"))
#' @export
print.autorelevate_fit <- function(x, ...) {
  cat("\nAutorelevated Family Model Fit\n")
  cat("------------------------------\n")
  cat("Family:   ", x$dist, "\n")
  cat("Method:   ", toupper(x$method), "\n")
  cat("Sample n: ", x$n, "\n\n")
  cat("Parameter Estimates:\n")
  print(x$par)
  invisible(x)
}

#' Summary Table for an Autorelevated Fit
#'
#' @description
#' Reports parameter estimates and, for \code{method = "mle"}, standard
#' errors, log-likelihood, AIC, BIC, CAIC, and HQIC (see Details in
#' \code{\link{fit_autorelevate}} for why these are MLE-only), together
#' with the Kolmogorov-Smirnov goodness-of-fit statistic for every
#' method.
#'
#' @param object An autorelevate_fit object.
#' @param ... Unused additional arguments.
#' @return Matrix of estimates and standard errors (invisibly).
#' @importFrom stats printCoefmat
#' @examples
#' set.seed(1)
#' x <- rautorelevate(100, dist = "weibull", p1 = 0.5, p2 = 1.5)
#' fit <- fit_autorelevate(x, dist = "weibull", method = "mle")
#' summary(fit)
#' @export
summary.autorelevate_fit <- function(object, ...) {
  cat("\nModel Fit Summary (Distribution: ", object$dist, ")\n", sep = "")

  if (object$method == "mle") {
    se <- suppressWarnings(sqrt(diag(object$vcov)))
    mat <- cbind(Estimate = object$par, `Std. Error` = se)
    cat("Estimation Method: MLE\n\n")
    stats::printCoefmat(mat)

    cat("\n--- Goodness-of-Fit & Model Selection ---\n")
    cat(sprintf("Log-Likelihood: %8.3f\n", object$log_lik))
    cat(sprintf("AIC:            %8.3f\n", object$aic))
    cat(sprintf("BIC:            %8.3f\n", object$bic))
    cat(sprintf("CAIC:           %8s\n", if (is.na(object$caic)) "NA" else sprintf("%.3f", object$caic)))
    cat(sprintf("HQIC:           %8.3f\n", object$hqic))
    cat(sprintf("KS Statistic:   %8.4f (p-value: %.4f)\n", object$ks_stat, object$ks_pval))
  } else {
    mat <- cbind(Estimate = object$par)
    cat("Estimation Method: ", toupper(object$method), "\n", sep="")
    cat("(Note: Standard Errors, AIC, BIC, CAIC, and HQIC are derived from\n")
    cat(" Fisher Information and are therefore strictly provided for MLE.)\n\n")
    stats::printCoefmat(mat)

    cat("\n--- Goodness-of-Fit ---\n")
    cat(sprintf("KS Statistic:   %8.4f (p-value: %.4f)\n", object$ks_stat, object$ks_pval))
  }
  invisible(mat)
}

#' Plot Diagnostics for an Autorelevated Fit
#'
#' @description
#' Two-panel diagnostic plot: a histogram of the data with the fitted
#' density overlaid, and the empirical CDF with the fitted CDF overlaid.
#'
#' @param x An autorelevate_fit object.
#' @param ... Unused additional arguments.
#' @return Invisibly, x.
#' @importFrom stats ecdf
#' @importFrom graphics par hist plot lines legend
#' @examples
#' set.seed(1)
#' x <- rautorelevate(200, dist = "weibull", p1 = 0.5, p2 = 1.5)
#' fit <- fit_autorelevate(x, dist = "weibull", method = "mle")
#' plot(fit)
#' @export
plot.autorelevate_fit <- function(x, ...) {
  oldpar <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(oldpar))

  grid <- seq(min(x$data), max(x$data), length.out = 100)

  theo_pdf <- dautorelevate(grid, dist = x$dist, p1 = x$par["p1"], p2 = x$par["p2"])
  graphics::hist(x$data, probability = TRUE, main = paste("Histogram vs PDF (", x$dist, ")", sep=""),
       xlab = "Data", ylab = "Density", col = "lightgray", border = "white")
  graphics::lines(grid, theo_pdf, col = "blue", lwd = 2)
  graphics::legend("topright", legend = c("Histogram", "Fitted PDF"),
         fill = c("lightgray", NA), border = c("white", NA),
         col = c(NA, "blue"), lwd = c(NA, 2), merge = TRUE)

  emp_cdf <- stats::ecdf(x$data)
  theo_cdf <- pautorelevate(grid, dist = x$dist, p1 = x$par["p1"], p2 = x$par["p2"])
  graphics::plot(emp_cdf, main = paste("Empirical vs CDF (", x$dist, ")", sep=""),
       xlab = "Data", ylab = "Cumulative Probability", do.points = FALSE)
  graphics::lines(grid, theo_cdf, col = "red", lwd = 2)
  graphics::legend("bottomright", legend = c("Empirical CDF", "Fitted CDF"),
         col = c("black", "red"), lwd = c(1, 2))

  invisible(x)
}

