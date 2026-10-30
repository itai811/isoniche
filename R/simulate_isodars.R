#' Simulate abundances in three habitats linked by isodars
#'
#' Simulates abundances at `n` sites: habitat 1 is drawn uniformly between 0
#' and `hab1_max`; habitat 2 = `int1` + `slope1` x habitat 1 + noise; habitat 3
#' = `int2` + `slope2` x habitat 2 + noise. Values are rounded to integers and
#' negative values set to 0.
#'
#' @param slope1,int1 Isodar of habitat 2 on habitat 1.
#' @param slope2,int2 Isodar of habitat 3 on habitat 2.
#' @param noise Standard deviation of the Gaussian noise.
#' @param n Number of sites.
#' @param hab1_max Upper bound of abundance in habitat 1. The default,
#'   `100 / slope1^1.2 - int1`, keeps total abundances in a comparable range
#'   across parameter values.
#' @return A tibble with columns `hab1`, `hab2`, `hab3` and `total_abundance`.
#' @examples
#' set.seed(1)
#' simulate_isodars(1, 1.1, 1, 0, noise = 2, n = 10)
#' @export
simulate_isodars <- function(slope1, slope2, int1, int2, noise = 1, n = 30,
                             hab1_max = 100 / (slope1^1.2) - int1) {
  stopifnot(is.numeric(noise), noise >= 0, n > 0)
  if (!is.finite(hab1_max) || hab1_max <= 0) {
    stop("'hab1_max' must be positive; try different parameter values.", call. = FALSE)
  }
  hab1 <- stats::runif(n, 0, hab1_max)
  hab2 <- slope1 * hab1 + int1 + stats::rnorm(n, 0, noise)
  hab3 <- slope2 * hab2 + int2 + stats::rnorm(n, 0, noise)
  X <- pmax(round(cbind(hab1, hab2, hab3)), 0)
  tibble::tibble(hab1 = X[, 1], hab2 = X[, 2], hab3 = X[, 3],
                 total_abundance = rowSums(X))
}
