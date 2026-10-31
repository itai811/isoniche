# isoniche 1.1.0

New maintainer: Itai Granot. Author order is now Granot, Dubiner, Belmaker;
Itai Granot and Shahar Dubiner contributed equally.

## Bug fixes

These change results relative to 1.0.0.

* The n-dimensional isodar (`fit_isodars(dim = "ndim")`) is now a line (the
  first principal axis), not a hyperplane (the smallest eigenvector) combined
  with "closest to equal allocation". The old approach biased niche breadth
  toward generalism.
* Inverse-variance weights (`weights = "1/var"`) were applied twice
  (effectively 1 / sd^4). They are now 1 / (sd^2 (1 + slope^2)), applied once.
* `isoderr(output = "SE")` divided the bootstrap standard deviation by
  `sqrt(n_boot)`. The bootstrap SD is itself the standard error, so "SD" and
  "SE" now both return it.
* `fit_isodars()` no longer recurses forever when an isodar has a negative
  slope and a negative intercept.
* Density-independent niche breadth (`iai(abundances = NULL)`) is now evaluated
  at the first integer total with at least `min_occupancy` (default 1)
  individual in every habitat. Previously the threshold was > 1e-9
  individuals. If no total qualifies, one row of `NA` is returned.
* `hab_cols = NULL` silently used every column of `data`. Now `hab_cols` (or
  `n_habitats`) is required.

## New features

* New `pairwise_niche()` summarises pairwise isodars into niche breadth
  indices (paper, Appendix S1): the mean intercept and mean slope, and the
  inverse root-mean-square deviations of slopes from 1 and of intercepts from
  0, with optional weights (`"1/var"` or numeric) and an option to measure
  slope deviations on the log scale (`log_slope = TRUE`).
* Pairwise isodars use every site where both habitats of a pair were recorded,
  so sites that lack some habitats still contribute.
* `fit_isodars(dim = "ndim", ndim_method = "latent")` estimates the line from
  incomplete sites (every site with at least two habitats) as a one-factor
  model fitted with glmmTMB, and checks that the habitats are connected through
  shared sites. `ndim_method = "auto"` picks `"pca"` when all sites are
  complete and `"latent"` otherwise.
* `fit_isodars()` returns classed tibbles (`isodars_pairwise`, `isodar_line`).
* `iai()` allocates abundance among habitats by weighted least squares on the
  pairwise isodar equations, constrained to the total; habitats predicted
  negative are set to 0 and the rest are re-solved.
* `isoderr(output = "CI")` returns percentile bootstrap intervals.
* `iai()` gives an informative error when `abundances` contains `NA`,
  non-finite or negative values.
* Unit tests cover the paper's reference values, a simulation with a known
  answer, each bug fix above, input validation and the latent method.

# isoniche 1.0.0

* First CRAN release (maintainer Shahar Dubiner).
