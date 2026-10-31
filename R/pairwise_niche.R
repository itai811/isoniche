#' Niche breadth indices from pairwise isodar parameters
#'
#' Summarises the intercepts and slopes of the pairwise isodars directly, as an
#' alternative to the isodar-adjusted inequality ([iai()]). A perfect
#' generalist has intercepts of 0 and slopes of 1 between all habitat pairs.
#'
#' @details
#' Each isodar is first expressed so that its intercept is non-negative: an
#' isodar \eqn{y = a + b x} with \eqn{a < 0} and \eqn{b > 0} is rewritten as
#' \eqn{x = -a/b + (1/b) y} (habitats swapped). [fit_isodars()] does this by
#' default (`flip_intercept = TRUE`); `pairwise_niche()` repeats it so that the
#' result does not depend on that choice.
#'
#' **Mean indices** are the (weighted) mean intercept and mean slope.
#'
#' **Inverse root-mean-square deviation indices** measure deviations from the
#' values that indicate no habitat preference (0 and 1), not from the mean.
#' The slope niche breadth is
#' \deqn{1 / \sqrt{\sum w_{ij} (b_{ij} - 1)^2 / \sum w_{ij}}}{1 / sqrt(sum(w * (b - 1)^2) / sum(w))}
#' and the intercept niche breadth is
#' \deqn{1 / \sqrt{\sum w_{ij} a_{ij}^2 / \sum w_{ij}}}{1 / sqrt(sum(w * a^2) / sum(w))}
#' where \eqn{a_{ij}} and \eqn{b_{ij}} are the intercept and slope of the
#' isodar between habitats \eqn{i} and \eqn{j}, and \eqn{w_{ij}} its weight.
#' Larger values indicate broader niches; extreme specialists approach 0, and
#' a perfect generalist gives `Inf`.
#'
#' Slopes are not symmetric around 1: slopes of 2 and 0.5 describe preferences
#' of equal strength in opposite directions. `log_slope = TRUE` measures slope
#' deviations as \eqn{(\log b)^2}{(log b)^2}, which removes this asymmetry. It
#' is undefined for non-positive slopes (for example, a habitat abandoned as
#' total abundance rises); the slope index is then `NA` and `note` says why.
#'
#' **Units.** The intercept index (and the mean intercept) is in units of
#' abundance per sampling unit, so it depends on survey effort and is
#' comparable among species only when surveys are standardised. The slope
#' indices are unitless.
#'
#' @param isodars Pairwise isodars from `fit_isodars(dim = "pairwise")`, or a
#'   data frame with columns `habitat_x`, `habitat_y`, `intercept`, `slope`
#'   (and `sd` for `weights = "1/var"`).
#' @param weights `NULL` (equal weights), `"1/var"` (inverse variance of each
#'   isodar, \eqn{1/(sd^2(1+b^2))}, computed after orienting the isodars), or a
#'   numeric vector with one non-negative weight per isodar.
#' @param log_slope If `TRUE`, slope deviations are \eqn{(\log b)^2}{(log b)^2}
#'   instead of \eqn{(b - 1)^2}.
#'
#' @return A one-row tibble with `mean_intercept`, `mean_slope`,
#'   `slope_niche_breadth`, `intercept_niche_breadth`, `n_isodars` and `note`.
#'   The oriented isodars are in attribute `"isodars"`.
#'
#' @seealso [fit_isodars()], [iai()]
#' @examples
#' set.seed(1)
#' sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
#' pw <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"))
#' pairwise_niche(pw)
#' pairwise_niche(pw, weights = "1/var", log_slope = TRUE)
#' @export
pairwise_niche <- function(isodars, weights = NULL, log_slope = FALSE) {
  if (!is.data.frame(isodars) ||
      !all(c("habitat_x", "habitat_y", "intercept", "slope") %in% names(isodars))) {
    stop("'isodars' must be pairwise isodars from fit_isodars().", call. = FALSE)
  }
  if (inherits(isodars, "isodar_line")) {
    stop("pairwise_niche() needs pairwise isodars (dim = \"pairwise\").",
         call. = FALSE)
  }
  if (nrow(isodars) == 0) {
    stop("'isodars' has no rows (fewer than two habitats used?).", call. = FALSE)
  }
  log_slope <- isTRUE(as.logical(log_slope))
  iso <- orient_isodars(isodars)
  a <- iso$intercept
  b <- iso$slope
  K <- nrow(iso)
  w <- if (is.null(weights)) {
    rep(1, K)
  } else if (identical(weights, "1/var")) {
    if (!"sd" %in% names(iso)) {
      stop("weights = \"1/var\" needs an 'sd' column in 'isodars'.", call. = FALSE)
    }
    1 / pmax(iso$sd^2 * (1 + b^2), 1e-12)
  } else if (is.numeric(weights) && length(weights) == K &&
             all(is.finite(weights)) && all(weights >= 0) && sum(weights) > 0) {
    weights
  } else {
    stop("'weights' must be NULL, \"1/var\" or one non-negative number per isodar.",
         call. = FALSE)
  }

  note <- ""
  if (log_slope && any(b <= 0)) {
    slope_dev <- rep(NA_real_, K)
    note <- "non-positive slope(s): log-slope index undefined"
  } else {
    slope_dev <- if (log_slope) log(b)^2 else (b - 1)^2
  }
  wmean <- function(x) sum(w * x) / sum(w)
  out <- tibble::tibble(
    mean_intercept = wmean(a),
    mean_slope = wmean(b),
    slope_niche_breadth = 1 / sqrt(wmean(slope_dev)),
    intercept_niche_breadth = 1 / sqrt(wmean(a^2)),
    n_isodars = K,
    note = note
  )
  attr(out, "isodars") <- iso
  out
}

# Express each isodar with a non-negative intercept (swap habitats if needed)
orient_isodars <- function(isodars) {
  iso <- as.data.frame(isodars, stringsAsFactors = FALSE)
  if (is.null(iso$flipped)) iso$flipped <- FALSE
  neg <- is.finite(iso$intercept) & iso$intercept < 0
  if (any(neg & !(iso$slope > 0))) {
    stop("An isodar with a negative intercept and a non-positive slope cannot ",
         "be expressed with a non-negative intercept.", call. = FALSE)
  }
  if (any(neg)) {
    hx <- iso$habitat_x[neg]
    iso$habitat_x[neg] <- iso$habitat_y[neg]
    iso$habitat_y[neg] <- hx
    iso$intercept[neg] <- -iso$intercept[neg] / iso$slope[neg]
    iso$slope[neg] <- 1 / iso$slope[neg]
    iso$flipped[neg] <- !iso$flipped[neg]
  }
  tibble::as_tibble(iso)
}
