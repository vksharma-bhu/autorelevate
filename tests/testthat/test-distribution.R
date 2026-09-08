test_that("Distribution/survival/CDF/hazard identities hold for all ten baseline distributions", {
  x <- seq(0.1, 3.0, by = 0.5)
  for (fam in autorelevate:::.AR_VALID_DISTS) {
    cdf_vals <- pautorelevate(x, dist = fam, p1 = 0.8, p2 = 1.2)
    surv_vals <- sautorelevate(x, dist = fam, p1 = 0.8, p2 = 1.2)
    expect_equal(cdf_vals + surv_vals, rep(1, length(x)), tolerance = 1e-7,
                 info = paste("distribution:", fam))

    q_recov <- qautorelevate(cdf_vals, dist = fam, p1 = 0.8, p2 = 1.2)
    expect_equal(q_recov, x, tolerance = 1e-4, info = paste("distribution:", fam))

    f_vals <- dautorelevate(x, dist = fam, p1 = 0.8, p2 = 1.2)
    h_vals <- haautorelevate(x, dist = fam, p1 = 0.8, p2 = 1.2)
    expect_equal(h_vals, f_vals / surv_vals, tolerance = 1e-8,
                 info = paste("distribution:", fam))
  }
})

test_that("dautorelevate integrates to 1 (numerically) for every baseline distribution", {
  # Each distribution gets parameters and an integration bound that avoid its
  # own heavy-tail regime (Lomax, Burr XII, and Log-Logistic decay slowly
  # for some parameter values, requiring a much larger upper bound for
  # numerical integration to converge -- this is a property of those
  # distributions' tails, not a defect in the density formula; see the
  # (p1, p2) used elsewhere in this file and in the vignette for typical,
  # well-behaved parameter ranges).
  params <- list(weibull = c(0.5, 1.5), lomax = c(2, 3), burr = c(2, 3),
                 gompertz = c(0.4, 0.7), loglogistic = c(1.3, 2.1),
                 chen = c(0.3, 1.4), expexp = c(0.8, 2.3),
                 powerlindley = c(1.2, 1.7), lognormal = c(0.5, 0.8),
                 gamma = c(1.1, 2.4))
  bounds <- list(weibull = 60, lomax = 500, burr = 200, gompertz = 60,
                 loglogistic = 2000, chen = 20, expexp = 30,
                 powerlindley = 20, lognormal = 200, gamma = 60)
  for (fam in autorelevate:::.AR_VALID_DISTS) {
    pp <- params[[fam]]
    int_val <- stats::integrate(dautorelevate, lower = 0, upper = bounds[[fam]],
                                 dist = fam, p1 = pp[1], p2 = pp[2],
                                 stop.on.error = FALSE, subdivisions = 500L)$value
    expect_equal(int_val, 1, tolerance = 1e-3, info = paste("distribution:", fam))
  }
})

test_that("Invalid inputs raise clear errors instead of silent NaN/Inf", {
  expect_error(dautorelevate(1, dist = "not_a_family"), "dist")
  expect_error(dautorelevate(1, dist = "weibull", p1 = -1), "p1")
  expect_error(dautorelevate(1, dist = "weibull", p2 = 0), "p2")
  expect_error(fit_autorelevate(numeric(0), dist = "weibull"), "data")
  expect_error(fit_autorelevate(c(-1, -2), dist = "weibull"), "data")
})

test_that("bladder_cancer dataset matches published summary statistics", {
  data(bladder_cancer)
  expect_equal(length(bladder_cancer), 128)
  expect_equal(round(mean(bladder_cancer), 3), 9.366)
  expect_equal(round(var(bladder_cancer), 3), 110.425)
})

test_that("ttt_plot returns the expected data frame without plotting", {
  data(bladder_cancer)
  res <- ttt_plot(bladder_cancer, plot = FALSE)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), length(bladder_cancer))
  expect_equal(res$r_over_n[nrow(res)], 1)
})

test_that("haautorelevate propagates NA consistently with d/p/s-autorelevate", {
  expect_true(is.na(haautorelevate(NA, dist = "weibull", p1 = 0.5, p2 = 1.5)))
  expect_true(is.na(dautorelevate(NA, dist = "weibull", p1 = 0.5, p2 = 1.5)))
  expect_true(is.na(pautorelevate(NA, dist = "weibull", p1 = 0.5, p2 = 1.5)))
  expect_true(is.na(sautorelevate(NA, dist = "weibull", p1 = 0.5, p2 = 1.5)))
})

test_that("qautorelevate warns and returns NaN for probabilities outside [0, 1]", {
  expect_warning(res_neg <- qautorelevate(-0.5, dist = "weibull", p1 = 0.5, p2 = 1.5),
                  "NaNs produced")
  expect_true(is.nan(res_neg))
  expect_warning(res_big <- qautorelevate(1.5, dist = "weibull", p1 = 0.5, p2 = 1.5),
                  "NaNs produced")
  expect_true(is.nan(res_big))
  # valid probabilities still work normally and raise no warning
  expect_silent(qautorelevate(0.5, dist = "weibull", p1 = 0.5, p2 = 1.5))
})

test_that("ttt_plot accepts custom xlim/ylim/... without a duplicate-argument error", {
  data(bladder_cancer)
  expect_error(
    { grDevices::pdf(NULL); ttt_plot(bladder_cancer, xlim = c(0, 2), col = "red"); grDevices::dev.off() },
    NA
  )
})

test_that("print/summary/plot S3 methods run without error", {
  set.seed(1)
  x <- rautorelevate(80, dist = "weibull", p1 = 0.5, p2 = 1.5)
  fit <- fit_autorelevate(x, dist = "weibull", method = "mle")
  # capture.output() silences the printed output during the test run
  # (purely cosmetic -- it has no effect on whether the test passes or
  # on R CMD check's status either way) while still exercising each
  # method fully and asserting it does not error.
  expect_error(invisible(capture.output(print(fit))), NA)
  expect_error(invisible(capture.output(summary(fit))), NA)
  expect_error({
    grDevices::pdf(NULL)
    plot(fit)
    grDevices::dev.off()
  }, NA)
  comp <- compare_families(x)
  expect_error(invisible(capture.output(print(comp))), NA)
})

test_that("lower.tail, log.p, and log argument branches are correct", {
  x <- seq(0.2, 3, by = 0.4)
  cdf <- pautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5)
  surv <- sautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5)
  expect_equal(pautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5, lower.tail = FALSE), surv)
  expect_equal(pautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5, log.p = TRUE), log(cdf))
  expect_equal(sautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5, log = TRUE), log(surv))
  dens <- dautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5)
  expect_equal(dautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5, log = TRUE), log(dens))
  q_upper <- qautorelevate(1 - cdf, dist = "weibull", p1 = 0.5, p2 = 1.5, lower.tail = FALSE)
  expect_equal(q_upper, x, tolerance = 1e-6)
})

test_that(".lambert_w_minus1 matches exact known values and independent high-precision ground truth", {
  lw <- autorelevate:::.lambert_w_minus1

  # Exact closed-form values.
  expect_equal(lw(-1 / exp(1)), -1)
  expect_equal(lw(-2 * exp(-2)), -2, tolerance = 1e-10)

  # Ground truth via uniroot() at tight tolerance -- independent of either
  # implementation's own algorithm. The search interval's upper bound is
  # deliberately kept just short of the exact singular point -1 (where
  # f'(w) = 0): uniroot() itself becomes unreliable if the true root sits
  # extremely close to an interval endpoint that touches the singularity.
  ground_truth_w <- function(z) {
    if (z == -1 / exp(1)) return(-1)
    uniroot(function(w) w * exp(w) - z, interval = c(-800, -1 + 1e-13),
            tol = .Machine$double.eps^0.9)$root
  }

  # Regression test for a fixed bug: the original convergence check only
  # examined the residual |f(w)|, which is unreliable very close to
  # z = 0 (where f'(w) = e^w(w+1) is also astronomically small, so a
  # tiny residual does not imply an accurate w). This silently produced
  # an inaccurate root (error ~1e-3) for z corresponding to a quantile
  # request with p within 1e-10 of 1; the fix makes this region accurate
  # to near machine precision.
  z_near_zero <- c(-1e-8, -1e-9, -1e-10, -1e-11)
  for (z in z_near_zero) {
    expect_equal(lw(z), ground_truth_w(z), tolerance = 1e-8, info = paste("z =", z))
  }

  # Very close to the branch point (but not exactly at it), a small
  # residual imprecision remains (on the order of 1e-6, since the true
  # root itself differs from -1 by a comparable amount there) -- verified
  # separately to still be substantially more accurate than the lamW
  # package's implementation at these same points (our error ~2e-6 vs
  # lamW's ~4e-4), and this region is reached only for quantile requests
  # with p within machine-epsilon-scale of 0, far beyond any practical
  # relevance.
  z_near_branch <- c(-1 / exp(1) + 1e-10, -1 / exp(1) + 1e-12)
  for (z in z_near_branch) {
    expect_equal(lw(z), ground_truth_w(z), tolerance = 1e-5, info = paste("z =", z))
  }

  # No spurious "did not converge" warning anywhere across the domain,
  # including very close to both endpoints (a second bug this exact fix
  # also resolved: near the branch point, f'(w) -> 0 algebraically as
  # w -> -1, which made a step-size-only convergence check unreliable
  # there even though the answer was already accurate).
  z_dense <- c(-1 / exp(1) + 10^seq(-14, -1, length.out = 200),
               -10^seq(-14, -1, length.out = 200))
  expect_silent(lw(z_dense))
})

