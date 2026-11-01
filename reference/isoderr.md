# Bootstrap uncertainty of the isodar-adjusted inequality

Resamples sites with replacement (keeping all habitats of a site
together), refits the isodars and recomputes the IAI for each resample.
Absent habitats are identified once, on the full data, and kept fixed in
every resample.

## Usage

``` r
isoderr(
  data,
  hab_cols = NULL,
  n_habitats = NULL,
  dim = "pairwise",
  abundances = NULL,
  n_boot = 1000,
  output = c("SD", "SE", "CI"),
  conf_level = 0.95,
  weights = NULL,
  flip_intercept = TRUE,
  min_occupancy = 1,
  max_search = 10000,
  seed = NULL,
  ndim_method = "auto",
  absent = "auto",
  absent_max_count = 0,
  absent_min_sites = 5
)
```

## Arguments

- data:

  A data frame. Wide format: one row per site and one numeric column per
  habitat; `NA` marks a habitat that is absent from, or was not sampled
  at, a site. Long format: one row per species x site x habitat (see
  Details).

- hab_cols:

  Wide format: character vector of habitat column names.

- n_habitats:

  Alternatively, the number of habitat columns, taken from the start of
  `data`. One of `hab_cols` or `n_habitats` is required.

- dim:

  `"pairwise"` or `"ndim"`.

- abundances:

  Numeric vector of total abundances, or `NULL` for the
  density-independent value.

- n_boot:

  Number of bootstrap replicates.

- output:

  `"SD"` or `"SE"` (both the bootstrap standard deviation, which is the
  standard error of the estimate) or `"CI"` (percentile interval).

- conf_level:

  Confidence level for `output = "CI"`.

- weights:

  Pairwise only: `NULL` (equal weights), `"1/var"` (inverse variance of
  each isodar's vertical residuals, \\1/(sd^2(1+b^2))\\), or a numeric
  vector with one weight per isodar.

- flip_intercept:

  Pairwise only. If `TRUE`, an isodar with a negative intercept and
  positive slope is expressed in the reverse direction (habitats
  swapped), so that intercepts are non-negative.

- min_occupancy:

  Individuals required in every habitat for the density-independent
  value.

- max_search:

  Largest total searched for the density-independent value.

- seed:

  Optional random-number seed.

- ndim_method:

  `"auto"`, `"pca"` or `"latent"`; see Details.

- absent:

  Habitats the species does not use: `"auto"` (detected from the data,
  see Details), `"none"`, or a character vector of habitat names
  (`character(0)` is the same as `"none"`).

- absent_max_count, absent_min_sites:

  Thresholds for `absent = "auto"`; see Details.

## Value

A tibble with `total_abundance`, `niche_breadth` (estimate from the full
data), `uncertainty` (SD/SE) or `lower` and `upper` (CI), and `n_valid`
(bootstrap replicates that produced a value).

## Examples

``` r
set.seed(1)
sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
isoderr(sim, hab_cols = c("hab1", "hab2", "hab3"), abundances = 100,
        n_boot = 50, seed = 1)
#> # A tibble: 1 × 4
#>   total_abundance niche_breadth uncertainty n_valid
#>             <dbl>         <dbl>       <dbl>   <dbl>
#> 1             100         0.656     0.00364      50
```
