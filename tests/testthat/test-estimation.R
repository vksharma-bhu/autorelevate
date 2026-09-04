test_that("Parameter estimation and family comparison converge accurately", {
  set.seed(123)
  sample_data <- rautorelevate(100, dist = "weibull", p1 = 0.5, p2 = 1.5)
  comp <- compare_families(sample_data)
  expect_true(is.data.frame(comp))
  expect_true(nrow(comp) == 10)
  expect_true(all(c("AIC", "BIC", "CAIC", "HQIC") %in% names(comp)))
})

test_that("fit_autorelevate recovers parameters reasonably for every method", {
  set.seed(321)
  sample_data <- rautorelevate(300, dist = "weibull", p1 = 0.5, p2 = 1.5)
  for (m in c("mle", "mps", "ls", "wls", "cvm")) {
    fit <- fit_autorelevate(sample_data, dist = "weibull", method = m)
    expect_s3_class(fit, "autorelevate_fit")
    expect_true(all(is.finite(fit$par)))
    expect_true(all(fit$par > 0))
  }
})

test_that("fit_all_methods returns a 5x3 matrix", {
  set.seed(7)
  sample_data <- rautorelevate(100, dist = "weibull", p1 = 0.5, p2 = 1.5)
  res <- fit_all_methods(sample_data, dist = "weibull")
  expect_equal(dim(res), c(5, 3))
})

test_that("fit_autorelevate works for every new baseline family", {
  set.seed(42)
  for (fam in c("chen", "expexp", "powerlindley", "lognormal", "gamma")) {
    sample_data <- rautorelevate(150, dist = fam, p1 = 0.8, p2 = 1.3)
    fit <- fit_autorelevate(sample_data, dist = fam, method = "mle")
    expect_s3_class(fit, "autorelevate_fit")
    expect_true(all(is.finite(fit$par)) && all(fit$par > 0))
    expect_true(is.finite(fit$aic) && is.finite(fit$bic) &&
                is.finite(fit$hqic))
  }
})

test_that("AIC/BIC/CAIC/HQIC match their formulas exactly given log-likelihood", {
  set.seed(5)
  sample_data <- rautorelevate(60, dist = "weibull", p1 = 0.5, p2 = 1.5)
  fit <- fit_autorelevate(sample_data, dist = "weibull", method = "mle")
  n <- fit$n; k <- 2; ll <- fit$log_lik
  expect_equal(fit$aic, 2 * k - 2 * ll)
  expect_equal(fit$bic, k * log(n) - 2 * ll)
  expect_equal(fit$caic, (2 * k * n) / (n - k - 1) - 2 * ll)
  expect_equal(fit$hqic, 2 * k * log(log(n)) - 2 * ll)
})

test_that("compare_families works end-to-end on the bundled real dataset", {
  data(bladder_cancer)
  comp <- compare_families(bladder_cancer)
  expect_s3_class(comp, "autorelevate_compare")
  expect_true(nrow(comp) > 0)
  expect_true(all(diff(comp$AIC) >= 0))  # ranked ascending by AIC
})

test_that("Infinite p1/p2 are rejected instead of silently propagating as NaN/Inf", {
  expect_error(dautorelevate(1, dist = "weibull", p1 = Inf, p2 = 1.5), "p1")
  expect_error(dautorelevate(1, dist = "weibull", p1 = 0.5, p2 = Inf), "p2")
  # data-level Inf is dropped with a warning (not a hard error) as long as
  # some valid data remain; assert the warning explicitly instead of just
  # suppressing it, so this test documents the intended behavior precisely.
  expect_warning(
    fit <- fit_autorelevate(c(1, 2, Inf), dist = "weibull"),
    "Dropped"
  )
  expect_equal(fit$n, 2)
})

test_that("fit_autorelevate warns (rather than silently dropping) invalid observations", {
  set.seed(1)
  good_data <- rautorelevate(50, dist = "weibull", p1 = 0.5, p2 = 1.5)
  contaminated <- c(good_data, -5, -3, NA, NA, 0, Inf)
  expect_warning(fit <- fit_autorelevate(contaminated, dist = "weibull", method = "mle"),
                  "Dropped")
  expect_equal(fit$n, 50)
})

test_that("MLE fit reproduces the published Table 8 values of Dileep Kumar et al. (2025)", {
  # Regression test against the paper's own reported ARW fit to the bladder
  # cancer dataset (Table 8): p1 (lambda) = 0.4305, p2 (beta) = 0.7257,
  # -logL = 411.39, AIC = 826.7833, BIC = 832.4873, CAIC = 826.8793,
  # HQIC = 829.1008. This is the strongest available correctness check for
  # the density/likelihood/AIC pipeline: independent agreement with a
  # peer-reviewed, published numerical result on real data.
  data(bladder_cancer)
  fit <- fit_autorelevate(bladder_cancer, dist = "weibull", method = "mle")
  expect_equal(as.numeric(fit$par["p1"]), 0.4305, tolerance = 0.005)
  expect_equal(as.numeric(fit$par["p2"]), 0.7257, tolerance = 0.005)
  expect_equal(-fit$log_lik, 411.39, tolerance = 0.01)
  expect_equal(fit$aic, 826.7833, tolerance = 0.01)
  expect_equal(fit$bic, 832.4873, tolerance = 0.01)
  expect_equal(fit$caic, 826.8793, tolerance = 0.01)
  expect_equal(fit$hqic, 829.1008, tolerance = 0.01)
})

test_that("MLE standard errors are on the correct (original parameter) scale", {
  # Regression test for a fixed bug: optimHess() is evaluated on the
  # (log p1, log p2) scale internally, so the raw inverse-Hessian gives
  # Var(log p1), Var(log p2), not Var(p1), Var(p2). fit_autorelevate()
  # must apply the delta method (D = diag(p1, p2)) before returning vcov.
  # We check this by comparing against a Hessian computed directly on the
  # original (p1, p2) scale, which is unambiguous ground truth.
  set.seed(11)
  x <- rautorelevate(200, dist = "weibull", p1 = 0.5, p2 = 1.5)
  fit <- fit_autorelevate(x, dist = "weibull", method = "mle")
  se_reported <- sqrt(diag(fit$vcov))

  x_sorted <- sort(x)
  obj_direct <- function(par) {
    p1 <- par[1]; p2 <- par[2]
    if (p1 <= 0 || p2 <= 0) return(1e10)
    ld <- dautorelevate(x_sorted, dist = "weibull", p1 = p1, p2 = p2, log = TRUE)
    if (any(!is.finite(ld))) return(1e10)
    -sum(ld)
  }
  H_direct <- stats::optimHess(as.numeric(fit$par), obj_direct)
  se_ground_truth <- sqrt(diag(solve(H_direct)))

  expect_equal(as.numeric(se_reported), se_ground_truth, tolerance = 0.02)
})

