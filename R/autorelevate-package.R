#' autorelevate: The Autorelevated Family of Probability Distributions and Estimation Methods
#'
#' @description
#' Implements the \strong{autorelevated family} of probability
#' distributions -- obtained by applying the autorelevation
#' transform of Krakowski (1973) to a baseline
#' lifetime distribution -- across ten baseline distributions (Weibull,
#' Lomax, Burr XII, Gompertz, Log-Logistic, Chen, Exponentiated
#' Exponential, Power Lindley, Log-normal, and Gamma). Provides
#' vectorized density, distribution, survival, hazard, exact quantile,
#' and random-generation functions for each distribution, together with five
#' estimation methods (Maximum Likelihood, Maximum Product of Spacings,
#' Least Squares, Weighted Least Squares, and Cramer-von Mises),
#' goodness-of-fit testing, a Total Time on Test diagnostic plot, and
#' cross-family model selection by AIC, BIC, CAIC, and HQIC. See
#' Details below for the mathematical derivation of the transform and
#' the exact role of \code{p1}/\code{p2} within each baseline distribution.
#'
#' @details
#' \strong{The relevation and autorelevation transformations.} Let \eqn{X} and
#' \eqn{Y} be non-negative random variables with survival functions
#' \eqn{\bar F_1} and \eqn{\bar F_2}, respectively. The \emph{relevation transform}
#' \eqn{\Psi(X)} of Krakowski (1973) is the survival function of the total
#' lifetime of a system in which an item from population 1 is replaced, at
#' the moment of its failure at age \eqn{x}, by an item of the same age
#' \eqn{x} from population 2. The survival function is
#' \deqn{\bar F(x) = \bar F_1(x) - \int_0^x \frac{\bar F_2(t)}{\bar F_1(t)} \, dF_1(t).}
#' When \eqn{X} and \eqn{Y} are identically distributed with common
#' baseline survival function \eqn{\bar F_0}, the resulting variable is
#' called the \emph{autorelevation} of \eqn{X}, studied systematically by
#' Dileepkumar and Sankaran (2022). Its survival function simplifies to
#' \deqn{\bar F(x) = \bar F_0(x)\left(1 - \log \bar F_0(x)\right) = q(\bar F_0(x)),}
#' where \eqn{q(t) = t(1 - \log t)} is a concave distortion function on
#' \eqn{[0, 1]} with \eqn{q(0) = 0} and \eqn{q(1) = 1}. Writing the baseline
#' cumulative hazard as \eqn{\Lambda_0(x) = -\log \bar F_0(x)}, this becomes
#' \deqn{\bar F(x) = e^{-\Lambda_0(x)}\left(1 + \Lambda_0(x)\right), \qquad
#'       f(x) = f_0(x)\,\Lambda_0(x),}
#' which is the parametrization used throughout this package
#' (\code{\link{sautorelevate}}, \code{\link{dautorelevate}}). Because the
#' transformation adds no extra parameter to the baseline model, every
#' autorelevated member inherits the baseline's parameters \code{p1},
#' \code{p2} with no increase in dimensionality.
#'
#' \strong{Domain.} Throughout this package, and in every formula in this
#' documentation, \eqn{x > 0} (support of the autorelevated family), and
#' both baseline parameters \eqn{p1 > 0} and \eqn{p2 > 0} -- reflecting
#' that every one of the ten baseline distributions is itself a lifetime
#' distribution with positive support and positive parameters. Values
#' of \code{x}, \code{p1}, or \code{p2} outside this range are not part
#' of the model; see \code{\link{.validate_ar_inputs}} for how invalid
#' \code{p1}/\code{p2} are rejected with an explicit error, and the
#' \code{d}/\code{p}/\code{s}/\code{ha} functions' own documentation for
#' how out-of-range \code{x} is handled (mapped to the natural boundary
#' value of \eqn{0}, \eqn{1}, or the density/hazard being \eqn{0}, as
#' appropriate, rather than raising an error).
#'
#' \strong{Baseline distributions.} Ten baseline probability distributions are
#' supported via the \code{dist} argument, each with its own
#' interpretation of \code{p1} and \code{p2} (both must be single
#' positive, finite scalars -- vectorizing over parameter values, the
#' way \code{\link[stats]{dnorm}} vectorizes over \code{mean}, is not
#' supported; only the first argument, e.g. \code{x} or \code{q}, is
#' vectorized):
#' \itemize{
#'   \item \code{"weibull"}: \eqn{\Lambda_0(x) = p1\, x^{p2}}. This is the
#'     Autorelevated Weibull distribution of Dileep Kumar, Shabeer,
#'     and Sankaran (2025).
#'   \item \code{"lomax"}: \eqn{\Lambda_0(x) = p2 \log(1 + p1\, x)}. The
#'     autorelevated Lomax distribution is studied by Sharma, Pal,
#'     Bhardwaj, and Tyagi (2026, submitted), who establish its
#'     upside-down bathtub hazard shape with applications to cancer
#'     survival data.
#'   \item \code{"burr"} (Burr XII): \eqn{\Lambda_0(x) = p2 \log(1 + x^{p1})}.
#'   \item \code{"gompertz"}: \eqn{\Lambda_0(x) = p1\left(e^{p2 x} - 1\right)}.
#'   \item \code{"loglogistic"}: \eqn{\Lambda_0(x) = \log\!\left(1 + (x/p1)^{p2}\right)}.
#'   \item \code{"chen"}: \eqn{\Lambda_0(x) = p1\left(e^{x^{p2}} - 1\right)}
#'     (Chen, 2000).
#'   \item \code{"expexp"} (Exponentiated Exponential):
#'     \eqn{\Lambda_0(x) = -\log\!\left(1 - (1 - e^{-p1 x})^{p2}\right)}
#'     (Gupta & Kundu, 1999).
#'   \item \code{"powerlindley"}: \eqn{\Lambda_0(x) = p1\,x^{p2} -
#'     \log\!\left(1 + \dfrac{p1\,x^{p2}}{p1+1}\right)} (Ghitany,
#'     Al-Mutairi, Balakrishnan, & Al-Enezi, 2013); inverted in closed
#'     form using the same Lambert \eqn{W_{-1}} solver used elsewhere in
#'     this package (see \code{\link{qautorelevate}}).
#'   \item \code{"lognormal"}: \eqn{\Lambda_0(x) = -\log \Phi\!\left(-
#'     \dfrac{\log x - p1}{p2}\right)}, \code{p1} = log-mean, \code{p2} =
#'     log-standard-deviation.
#'   \item \code{"gamma"}: \code{p1} = rate, \code{p2} = shape, using
#'     \code{\link[stats]{pgamma}}/\code{\link[stats]{qgamma}} internally.
#' }
#'
#' \strong{Estimation.} \code{\link{fit_autorelevate}} fits any of the ten
#' members by Maximum Likelihood (MLE), Maximum Product of Spacings
#' (MPS), Least Squares (LS), Weighted Least Squares (WLS), or
#' Cramer-von Mises (CvM) minimum-distance estimation, and reports AIC,
#' BIC, CAIC, and HQIC for \code{method = "mle"}.
#' \code{\link{fit_all_methods}} and \code{\link{compare_families}}
#' compare methods and baseline distributions, respectively.
#' \code{\link{ttt_plot}} produces a Total Time on Test plot, a
#' pre-fitting diagnostic for hazard shape (Aarset, 1987).
#'
#' \strong{Data.} \code{\link{bladder_cancer}} bundles a widely used
#' real lifetime dataset (Lee & Wang, 2003) for illustration.
#'
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
#'
#' Sharma, V. K., Pal, S., Bhardwaj, H., & Tyagi, V. (2026). The
#' autorelevated Lomax distribution: An upside-down bathtub hazard model
#' with properties and applications to cancer survival data. Submitted.
#'
#' Chen, Z. (2000). A new two-parameter lifetime distribution with
#' bathtub shape or increasing failure rate function. \emph{Statistics &
#' Probability Letters}, 49(2), 155-161.
#'
#' Gupta, R. D., & Kundu, D. (1999). Generalized exponential
#' distributions. \emph{Australian & New Zealand Journal of Statistics},
#' 41(2), 173-188.
#'
#' Ghitany, M. E., Al-Mutairi, D. K., Balakrishnan, N., & Al-Enezi, L. J.
#' (2013). Power Lindley distribution and associated inference.
#' \emph{Computational Statistics & Data Analysis}, 64, 20-33.
#'
#' Aarset, M. V. (1987). How to identify a bathtub hazard rate. \emph{IEEE
#' Transactions on Reliability}, R-36(1), 106-108.
#'
#' Lee, E. T., & Wang, J. W. (2003). \emph{Statistical Methods for
#' Survival Data Analysis} (3rd ed.). Wiley.
#'
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

