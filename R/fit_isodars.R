#' Fit isodars
#'
#' Fits isodars to abundances (or densities) recorded in two or more habitats
#' at a set of sites (or other replicates, such as years), for a single
#' species. Data are in wide format: one row per site, one column per habitat.
#'
#' @details
#' **Pairwise isodars** (`dim = "pairwise"`) are fitted for every pair of
#' habitats by total least squares (TLS), which is identical to major-axis
#' regression in two dimensions. Each pair uses every site where both habitats
#' were recorded (non-missing), so sites lacking some habitats still contribute.
#'
#' **The n-dimensional isodar** (`dim = "ndim"`) is a line. Under the ideal
#' free distribution, equal fitness across \eqn{n} habitats imposes \eqn{n-1}
#' equalities, so abundances across habitats change along a single line as
#' total abundance changes, whatever the number of habitats. The line is
#' described by a point (`mean`) and a unit direction vector (`direction`),
#' with one component per habitat. Each component is the rate at which that
#' habitat gains individuals as total abundance rises (the analogue of the
#' isodar slope); where the line enters each habitat reflects its
#' density-independent preference (the analogue of the intercept). Every
#' pairwise isodar can be derived from it, and with two habitats it is the
#' classic isodar.
#'
#' With `ndim_method = "pca"` the line is the first principal axis of the
#' covariance matrix of complete sites (the TLS line). With
#' `ndim_method = "latent"` it is estimated from all sites with at least two
#' recorded habitats, as a one-factor model,
#' \eqn{N_{hs} = \mu_h + v_h t_s + \varepsilon_{hs}}, where \eqn{t_s} is the
#' position of site \eqn{s} along the isodar (requires \pkg{glmmTMB}).
#' `ndim_method = "auto"` uses `"pca"` when every site is complete and
#' `"latent"` otherwise.
#'
#' @param data A data frame with one row per site and one numeric column per
#'   habitat. `NA` marks a habitat that is absent from, or was not sampled at, a
#'   site.
#' @param hab_cols Character vector of habitat column names.
#' @param n_habitats Alternatively, the number of habitat columns, taken from
#'   the start of `data`. One of `hab_cols` or `n_habitats` is required.
#' @param dim `"pairwise"` or `"ndim"`.
#' @param flip_intercept Pairwise only. If `TRUE`, an isodar with a negative
#'   intercept and positive slope is expressed in the reverse direction
#'   (habitats swapped), so that intercepts are non-negative.
#' @param ndim_method `"auto"`, `"pca"` or `"latent"`; see Details.
#'
#' @return For `dim = "pairwise"`, a tibble of class `isodars_pairwise` with
#'   one row per habitat pair: `habitat_x`, `habitat_y`, `intercept`, `slope`,
#'   `sd` (SD of perpendicular residuals), `n_sites`, and `flipped` (whether
#'   the habitats were swapped). For `dim = "ndim"`, a tibble of class
#'   `isodar_line` with one row per habitat: `habitat`, `mean`, `direction`;
#'   attributes `var_explained` (share of variance along the line, `"pca"`
#'   only), `eigenvalues`, `n_sites` and `ndim_method`.
#'
#' @references Morris, D.W. (1988). Habitat-dependent population regulation
#'   and community structure. *Evolutionary Ecology*, 2, 253–269.
#' @seealso [iai()], [isoderr()]
#' @examples
#' set.seed(1)
#' sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
#' fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"))
#' fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"), dim = "ndim")
#' @export
fit_isodars <- function(data, hab_cols = NULL, n_habitats = NULL,
                        dim = c("pairwise", "ndim"), flip_intercept = TRUE,
                        ndim_method = c("auto", "pca", "latent")) {
  dim <- match.arg(dim)
  ndim_method <- match.arg(ndim_method)
  X <- habitat_matrix(data, hab_cols, n_habitats)
  if (dim == "pairwise") {
    fit_pairwise(X, flip_intercept)
  } else {
    fit_line(X, ndim_method)
  }
}

# --- helpers ----------------------------------------------------------------

habitat_matrix <- function(data, hab_cols, n_habitats) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  if (is.null(hab_cols)) {
    if (is.null(n_habitats)) {
      stop("Specify the habitat columns with 'hab_cols' (or 'n_habitats').",
           call. = FALSE)
    }
    if (n_habitats > ncol(data)) stop("'n_habitats' exceeds the number of columns.", call. = FALSE)
    hab_cols <- names(data)[seq_len(n_habitats)]
  }
  missing_cols <- setdiff(hab_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Habitat column(s) not found: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (length(hab_cols) < 2) stop("At least two habitats are required.", call. = FALSE)
  X <- data[, hab_cols, drop = FALSE]
  if (!all(vapply(X, is.numeric, logical(1)))) {
    stop("Habitat columns must be numeric: ", paste(hab_cols, collapse = ", "), call. = FALSE)
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (any(X < 0, na.rm = TRUE)) stop("Abundances must be non-negative.", call. = FALSE)
  X
}

# TLS (major-axis) line for one pair
tls_pair <- function(x, y) {
  mx <- mean(x); my <- mean(y)
  sxx <- stats::var(x); syy <- stats::var(y); sxy <- stats::cov(x, y)
  theta <- 0.5 * atan2(2 * sxy, sxx - syy)
  slope <- tan(theta)
  intercept <- my - slope * mx
  resid <- ((y - my) * cos(theta) - (x - mx) * sin(theta))
  list(intercept = intercept, slope = slope,
       sd = sqrt(sum(resid^2) / max(length(x) - 2, 1)))
}

fit_pairwise <- function(X, flip_intercept = TRUE) {
  habs <- colnames(X)
  out <- list()
  for (i in seq_len(ncol(X) - 1)) {
    for (j in (i + 1):ncol(X)) {
      ok <- stats::complete.cases(X[, c(i, j)])
      x <- X[ok, i]; y <- X[ok, j]
      if (sum(ok) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
        warning("Isodar ", habs[i], "-", habs[j],
                " not fitted (fewer than 3 shared sites or no variation).", call. = FALSE)
        next
      }
      f <- tls_pair(x, y)
      hx <- habs[i]; hy <- habs[j]; flipped <- FALSE
      if (flip_intercept && is.finite(f$intercept) && f$intercept < 0 &&
          is.finite(f$slope) && f$slope > 0) {
        # TLS is symmetric in x and y: no refit needed
        f$intercept <- -f$intercept / f$slope
        f$slope <- 1 / f$slope
        hx <- habs[j]; hy <- habs[i]; flipped <- TRUE
      }
      out[[length(out) + 1]] <- data.frame(
        habitat_x = hx, habitat_y = hy, intercept = f$intercept, slope = f$slope,
        sd = f$sd, n_sites = sum(ok), flipped = flipped, stringsAsFactors = FALSE)
    }
  }
  if (length(out) == 0) stop("No pairwise isodar could be fitted.", call. = FALSE)
  res <- tibble::as_tibble(do.call(rbind, out))
  attr(res, "habitats") <- habs
  class(res) <- c("isodars_pairwise", class(res))
  res
}

fit_line <- function(X, ndim_method = "auto") {
  habs <- colnames(X)
  complete <- stats::complete.cases(X)
  if (ndim_method == "auto") ndim_method <- if (all(complete)) "pca" else "latent"
  if (ndim_method == "pca") {
    Xc <- X[complete, , drop = FALSE]
    if (nrow(Xc) < ncol(X) + 1) {
      stop("Not enough complete sites for ndim_method = 'pca'; try 'latent'.", call. = FALSE)
    }
    mu <- colMeans(Xc)
    e <- eigen(stats::cov(Xc), symmetric = TRUE)
    v <- e$vectors[, 1]
    ev <- e$values
    var_exp <- ev[1] / sum(ev)
    n_used <- nrow(Xc)
  } else {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("ndim_method = 'latent' requires the 'glmmTMB' package.", call. = FALSE)
    }
    long <- data.frame(site = factor(rep(seq_len(nrow(X)), ncol(X))),
                       habitat = factor(rep(habs, each = nrow(X)), levels = habs),
                       y = as.vector(X))
    long <- long[!is.na(long$y), ]
    keep <- names(which(table(long$site) >= 2))
    long <- droplevels(long[long$site %in% keep, ])
    check_connected(X[as.integer(keep), , drop = FALSE])
    m <- glmmTMB::glmmTMB(y ~ 0 + habitat + rr(0 + habitat | site, d = 1), data = long)
    L <- m$obj$env$report(m$fit$parfull)$fact_load[[1]][, 1]
    v <- L / sqrt(sum(L^2))
    mu <- as.numeric(glmmTMB::fixef(m)$cond)
    ev <- NA_real_; var_exp <- NA_real_
    n_used <- nlevels(long$site)
  }
  if (sum(v) < 0) v <- -v
  res <- tibble::tibble(habitat = habs, mean = as.numeric(mu), direction = as.numeric(v))
  attr(res, "var_explained") <- var_exp
  attr(res, "eigenvalues") <- ev
  attr(res, "n_sites") <- n_used
  attr(res, "ndim_method") <- ndim_method
  class(res) <- c("isodar_line", class(res))
  res
}

# habitats must be linked through sites where they co-occur
check_connected <- function(X) {
  H <- ncol(X)
  co <- crossprod(!is.na(X)) > 0
  reached <- 1
  repeat {
    new <- unique(c(reached, which(colSums(co[reached, , drop = FALSE]) > 0)))
    if (length(new) == length(reached)) break
    reached <- new
  }
  if (length(reached) < H) {
    stop("Habitats are not connected through shared sites: ",
         paste(colnames(X)[-reached], collapse = ", "),
         " never co-occur with the others.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @export
print.isodar_line <- function(x, ...) {
  cat("n-dimensional isodar (line), method:", attr(x, "ndim_method"),
      "| sites:", attr(x, "n_sites"), "\n")
  ve <- attr(x, "var_explained")
  if (!is.null(ve) && is.finite(ve)) cat("Variance along the line:", round(ve, 3), "\n")
  NextMethod()
}
