#' Total Time on Test (TTT) Plot
#'
#' @description
#' Produces the empirical Total Time on Test (TTT) plot of Aarset
#' (1987), a standard graphical diagnostic for identifying the shape of
#' a lifetime data set's hazard rate before fitting any model.
#'
#' @details
#' For ordered data \eqn{x_{(1)} \le \dots \le x_{(n)}}, the empirical TTT
#' statistic is
#' \deqn{G(r/n) = \left(\sum_{i=1}^r x_{(i)} + (n-r)x_{(r)}\right)
#'   \Big/ \sum_{i=1}^n x_{(i)}, \qquad r = 1, \dots, n.}
#' Plotting \eqn{G(r/n)} against \eqn{r/n} and comparing to the 45-degree
#' line indicates the hazard shape: concave indicates an increasing
#' hazard (IHR), convex indicates a decreasing hazard (DHR), convex then
#' concave indicates a bathtub hazard, and concave then convex indicates
#' an upside-down bathtub (UBT) hazard (Aarset, 1987). This is exactly
#' the diagnostic used by Dileep Kumar, Shabeer, and Sankaran (2025,
#' Sec. 8.1) to establish that the bundled \code{\link{bladder_cancer}}
#' dataset has UBT-shaped hazard, motivating the Autorelevated Weibull
#' fit in that section.
#'
#' @param data Vector of positive sample observations. Missing, non-finite,
#'   or non-positive values are dropped with a warning before computing
#'   the plot.
#' @param plot Logical; if TRUE (default) draws the plot. If FALSE, only
#'   the underlying data frame is returned (invisibly either way).
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#' @return Invisibly, a data frame with columns \code{r_over_n} and
#'   \code{TTT}.
#' @importFrom graphics plot abline
#' @importFrom utils modifyList
#' @references
#' Aarset, M. V. (1987). How to identify a bathtub hazard rate. \emph{IEEE
#' Transactions on Reliability}, R-36(1), 106-108.
#' @examples
#' data(bladder_cancer)
#' ttt_plot(bladder_cancer)   # concave-then-convex: UBT hazard
#' @export
ttt_plot <- function(data, plot = TRUE, ...) {
  n_input <- length(data)
  keep <- !is.na(data) & is.finite(data) & data > 0
  n_dropped <- n_input - sum(keep)
  if (n_dropped > 0) {
    warning(sprintf(
      "Dropped %d of %d observation(s) that were missing, non-finite, or not strictly positive; %d observation(s) used.",
      n_dropped, n_input, sum(keep)), call. = FALSE)
  }
  data <- data[keep]
  if (length(data) == 0) stop("`data` must contain at least one strictly positive, finite value.", call. = FALSE)
  x <- sort(data)
  n <- length(x)
  total <- sum(x)
  r <- seq_len(n)
  csum <- cumsum(x)
  G <- (csum + (n - r) * x) / total
  scaled_r <- r / n

  if (plot) {
    plot_args <- utils::modifyList(
      list(x = scaled_r, y = G, type = "l", xlab = "r/n", ylab = "TTT(r/n)",
           main = "Total Time on Test (TTT) Plot", xlim = c(0, 1), ylim = c(0, 1)),
      list(...)
    )
    do.call(graphics::plot, plot_args)
    graphics::abline(0, 1, lty = 2, col = "gray50")
  }
  invisible(data.frame(r_over_n = scaled_r, TTT = G))
}

