test_that("Autorelevated distribution identities hold mathematically", {
  x <- seq(0.1, 3.0, by = 0.5)
  cdf_vals <- pautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5)
  surv_vals <- sautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5)
  expect_equal(cdf_vals + surv_vals, rep(1, length(x)), tolerance = 1e-10)

  q_recov <- qautorelevate(cdf_vals, dist = "weibull", p1 = 0.5, p2 = 1.5)
  expect_equal(q_recov, x, tolerance = 1e-6)

  f_vals <- dautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5)
  h_vals <- haautorelevate(x, dist = "weibull", p1 = 0.5, p2 = 1.5)
  expect_equal(h_vals, f_vals / surv_vals, tolerance = 1e-8)
})

test_that("Distribution/survival/CDF identities hold for all ten baseline families", {
  x <- seq(0.1, 3.0, by = 0.5)
  for (fam in autorelevate:::.AR_VALID_DISTS) {
    cdf_vals <- pautorelevate(x, dist = fam, p1 = 0.8, p2 = 1.2)
    surv_vals <- sautorelevate(x, dist = fam, p1 = 0.8, p2 = 1.2)
    expect_equal(cdf_vals + surv_vals, rep(1, length(x)), tolerance = 1e-7,
                 info = paste("family:", fam))

    q_recov <- qautorelevate(cdf_vals, dist = fam, p1 = 0.8, p2 = 1.2)
    expect_equal(q_recov, x, tolerance = 1e-4, info = paste("family:", fam))
  }
})

test_that("dautorelevate integrates to 1 (numerically) for every baseline family", {
  # Each family gets parameters and an integration bound that avoid its
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
    expect_equal(int_val, 1, tolerance = 1e-3, info = paste("family:", fam))
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
  expect_equal(round(median(bladder_cancer), 3), 6.395)
  expect_equal(round(var(bladder_cancer), 3), 110.425)
  expect_equal(min(bladder_cancer), 0.08)
  expect_equal(max(bladder_cancer), 79.05)
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

