#' Isodar-adjusted inequality (IAI)
#'
#' Computes niche breadth from fitted isodars by allocating a given total
#' abundance among habitats according to the isodars and measuring the
#' inequality of that allocation (1 - Gini). Evaluating it at several totals
#' gives a rarefaction-like curve of niche breadth against abundance.
#'
#' @details
#' **Allocation.** For pairwise isodars, the predicted abundances are those that
#' minimise the (optionally weighted) sum of squared deviations from all
#' pairwise isodar equations, \eqn{y - b x = a}, subject to summing to the
#' total. This constrained least-squares problem is solved directly; habitats
#' predicted to be negative are set to zero and the rest are solved again. For
#' an n-dimensional isodar (a line), the allocation is the single point at
#' which the line reaches the total; negative habitats are set to zero.
#'
#' **Density-independent niche breadth** (`abundances = NULL`) is evaluated at
#' the smallest total (searched over integers) at which every used habitat
#' holds at least `min_occupancy` individuals. If no total up to `max_search`
#' achieves this, a single row with `NA` is returned.
#'
#' **Absent habitats** (see [fit_isodars()]) are fixed at 0 and count as 0 in
#' the inequality, so they lower niche breadth; with a single used habitat,
#' niche breadth is 0.
#'
#' **Density-dependent niche breadth** is the value at a high total abundance,
#' supplied through `abundances`.
#'
#' @param isodars Output of [fit_isodars()].
#' @param abundances Numeric vector of total abundances, or `NULL` for the
#'   density-independent value.
#' @param weights Pairwise only: `NULL` (equal weights), `"1/var"` (inverse
#'   variance of each isodar's vertical residuals, \eqn{1/(sd^2(1+b^2))}), or a
#'   numeric vector with one weight per isodar.
#' @param scaled If `TRUE`, 1 - Gini is rescaled to 0-1, so that a species
#'   found in a single habitat scores 0 and an even spread scores 1.
#' @param min_occupancy Individuals required in every habitat for the
#'   density-independent value.
#' @param max_search Largest total searched for the density-independent value.
#' @param dim Optional; inferred from `isodars`. Kept for compatibility.
#' @param plot If `TRUE` and several abundances are given, plots the curve.
#' @param method Inequality measure; only `"gini"` is currently implemented.
#'
#' @return A tibble with `total_abundance`, `niche_breadth`, one column per
#'   habitat (predicted abundance), and `note`.
#' @seealso [fit_isodars()], [isoderr()], [pairwise_niche()]
#' @examples
#' set.seed(1)
#' sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
#' pw <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"))
#' iai(pw)                          # density-independent
#' iai(pw, abundances = c(50, 500)) # density-dependent
#' ln <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"), dim = "ndim")
#' iai(ln, abundances = 500)
#' @export
iai <- function(isodars, abundances = NULL, weights = NULL, scaled = TRUE,
                min_occupancy = 1, max_search = 10000, dim = NULL,
                plot = FALSE, method = "gini") {
  method <- match.arg(method, "gini")
  scaled <- isTRUE(as.logical(scaled))
  if (!is.null(abundances) &&
      (!is.numeric(abundances) || !all(is.finite(abundances)) || any(abundances < 0))) {
    stop("'abundances' must be NULL or finite, non-negative numbers.", call. = FALSE)
  }
  kind <- if (inherits(isodars, "isodar_line")) "ndim" else
    if (inherits(isodars, "isodars_pairwise")) "pairwise" else NA
  if (is.na(kind)) {
    # plain data frames from older code
    kind <- if (all(c("habitat_x", "habitat_y") %in% names(isodars))) "pairwise" else
      stop("'isodars' must be the output of fit_isodars().", call. = FALSE)
  }
  if (!is.null(dim) && !identical(dim, kind)) {
    stop("'dim = \"", dim, "\"' does not match the fitted isodars (", kind, ").",
         call. = FALSE)
  }
  alloc <- if (kind == "pairwise") pairwise_allocator(isodars, weights) else
    line_allocator(isodars)
  habitats <- attr(alloc, "habitats")
  H <- length(habitats)
  absent <- attr(isodars, "absent")
  used <- !habitats %in% absent
  absent_note <- if (any(!used)) {
    paste0("absent habitat(s) fixed at 0: ", paste(habitats[!used], collapse = ", "))
  } else ""

  if (is.null(abundances)) {
    first_ok <- NA_real_
    for (Total in seq_len(max_search)) {
      x <- alloc(Total)
      if (all(is.finite(x)) && all(x[used] >= min_occupancy - 1e-9)) {
        first_ok <- Total
        break
      }
    }
    if (is.na(first_ok)) {
      note <- paste0("no total up to ", max_search, " occupies all ",
                     if (any(!used)) "used " else "", "habitats")
      out <- tibble::as_tibble(c(list(total_abundance = NA_real_, niche_breadth = NA_real_),
                                 stats::setNames(as.list(rep(NA_real_, H)), habitats),
                                 list(note = note)))
      return(out)
    }
    abundances <- first_ok
  }
  rows <- lapply(as.numeric(abundances), function(Total) {
    x <- alloc(Total)
    nb <- 1 - gini(x)
    if (scaled) nb <- min(max((nb - 1 / H) * H / (H - 1), 0), 1) # rounding
    note <- if (any(x[used] == 0)) "habitat(s) at zero" else ""
    note <- paste(c(note, absent_note)[c(note, absent_note) != ""], collapse = "; ")
    c(list(total_abundance = Total, niche_breadth = nb),
      stats::setNames(as.list(x), habitats),
      list(note = note))
  })
  out <- tibble::as_tibble(do.call(rbind.data.frame, rows))
  if (plot && nrow(out) > 1) {
    graphics::plot(out$total_abundance, out$niche_breadth, type = "l",
                   xlab = "Total abundance", ylab = "Niche breadth (IAI)")
  }
  out
}

gini <- function(x) {
  x <- sort(as.numeric(x))
  n <- length(x)
  if (sum(x) <= 0) return(NA_real_)
  (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
}

# Allocation along an n-dimensional isodar line
line_allocator <- function(fit) {
  mu <- fit$mean; v <- fit$direction; habs <- fit$habitat
  f <- function(Total) {
    x <- mu + (Total - sum(mu)) / sum(v) * v
    if (any(x < 0)) {
      x <- pmax(x, 0)
      if (sum(x) > 0) x <- x * Total / sum(x)
    }
    stats::setNames(x, habs)
  }
  attr(f, "habitats") <- habs
  f
}

# Constrained least-squares allocation from pairwise isodars
pairwise_allocator <- function(isodars, weights = NULL) {
  habs <- attr(isodars, "habitats")
  if (is.null(habs)) habs <- unique(c(isodars$habitat_x, isodars$habitat_y))
  H <- length(habs); K <- nrow(isodars)
  A <- matrix(0, K, H, dimnames = list(NULL, habs))
  for (k in seq_len(K)) {
    A[k, isodars$habitat_y[k]] <- 1
    A[k, isodars$habitat_x[k]] <- -isodars$slope[k]
  }
  b <- isodars$intercept
  w <- if (is.null(weights)) rep(1, K) else if (identical(weights, "1/var")) {
    1 / pmax(isodars$sd^2 * (1 + isodars$slope^2), 1e-12)
  } else if (is.numeric(weights) && length(weights) == K) weights else
    stop("'weights' must be NULL, \"1/var\" or one number per isodar.", call. = FALSE)
  solve_free <- function(Total, free) {
    Af <- A[, free, drop = FALSE]
    m <- length(free)
    Q <- crossprod(Af, Af * w)
    r <- crossprod(Af, b * w)
    M <- rbind(cbind(2 * Q, 1), c(rep(1, m), 0))
    s <- tryCatch(solve(M, c(2 * r, Total)), error = function(e) {
      qr.solve(M, c(2 * r, Total))
    })
    s[seq_len(m)]
  }
  # absent habitats are fixed at 0
  start <- which(!habs %in% attr(isodars, "absent"))
  f <- function(Total) {
    free <- start
    x <- numeric(H)
    for (iter in seq_len(H)) {
      xf <- solve_free(Total, free)
      if (all(xf >= 0) || length(free) == 1) { x[free] <- pmax(xf, 0); break }
      free <- free[xf >= 0]
      if (length(free) == 0) break
    }
    if (sum(x) > 0) x <- x * Total / sum(x)
    stats::setNames(x, habs)
  }
  attr(f, "habitats") <- habs
  f
}
