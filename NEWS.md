# autorelevate 0.1.0

* Initial CRAN submission.
* Implemented the autorelevated family of distributions (density, CDF,
  survival, hazard, exact quantile via Lambert W_-1, and random
  generation) for ten two-parameter baseline lifetime families: Weibull,
  Lomax, Burr XII, Gompertz, Log-Logistic, Chen, Exponentiated
  Exponential, Power Lindley, Log-normal, and Gamma (`dautorelevate()`,
  `pautorelevate()`, `sautorelevate()`, `haautorelevate()`,
  `qautorelevate()`, `rautorelevate()`).
* Added Maximum Likelihood, Maximum Product of Spacings, Least Squares,
  Weighted Least Squares, and Cramer-von Mises estimation via
  `fit_autorelevate()`, reporting AIC, BIC, CAIC, and HQIC for MLE fits.
* Added `fit_all_methods()` and `compare_families()` for comparing
  estimation methods and baseline families.
* Added `ttt_plot()` for Total Time on Test plots (Aarset, 1987), a
  standard pre-fitting diagnostic for hazard shape.
* Added `print()`, `summary()`, and `plot()` methods for
  `autorelevate_fit` objects, and a `print()` method for
  `autorelevate_compare` objects.
* Added a bundled real dataset, `bladder_cancer` (Lee & Wang, 2003;
  n = 128 bladder cancer remission times), with no dependency on
  external data packages. The MLE fit of the Autorelevated Weibull
  family to this dataset reproduces the published values of Dileep
  Kumar, Shabeer, and Sankaran (2025, Table 8) to within optimizer
  tolerance.
* Added an internal input-validation layer used by every exported
  function, giving clear error messages for invalid `dist`, `p1`,
  `p2`, or `data` instead of silent `NaN`/`Inf` propagation, and a
  warning (rather than silent truncation) when non-finite or
  non-positive observations are dropped from `data`.
* Made optimization in `fit_autorelevate()` robust via a data-adaptive
  starting-value grid search and an automatic Nelder-Mead/restart
  fallback if the primary optimizer fails to converge.
* Fully vectorized the internal Lambert W_-1 solver (`.lambert_w_minus1()`)
  for performance at scale, with a convergence check and warning.

