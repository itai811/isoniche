# Density-independent and density-dependent niche breadth

Fits isodars (unless given fitted ones) and returns, for each species,
the density-independent (DI) and density-dependent (DD) niche breadth
with bootstrap confidence intervals, the spillover abundance, preference
reversals and diagnostic flags. With a `species` column, everything runs
per species and the result has one row per species; otherwise the data
are one species.

## Usage

``` r
niche_breadth(
  fit_or_data,
  ...,
  dependent_at = "asymptote",
  min_occupancy = 1,
  reference_effort = NULL,
  absent = "auto",
  n_boot = 999,
  conf_level = 0.95,
  scaled = TRUE,
  weights = "1/var",
  diagnostics = list(),
  seed = NULL,
  max_search = 10000
)
```

## Arguments

- fit_or_data:

  Data in long or wide format (see
  [`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md)),
  or the long-format output of
  [`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md)
  (an `isodar_fits` table).

- ...:

  Passed to
  [`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md)
  when `fit_or_data` is data: the column names (`species`, `site`,
  `habitat`, `abundance`, `effort`) or `hab_cols`, and `dim`,
  `ndim_method`, `flip_intercept`, `absent_max_count`,
  `absent_min_sites`.

- dependent_at:

  `"asymptote"`, `"max_observed"` or a total abundance; see Details.

- min_occupancy:

  Individuals required in every used habitat at the spillover abundance.

- reference_effort, absent:

  Passed to
  [`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md)
  (data only).

- n_boot:

  Number of bootstrap replicates (0 for none).

- conf_level:

  Confidence level of the bootstrap intervals.

- scaled:

  If `TRUE`, niche breadth is scaled to 0-1 (see
  [`iai()`](https://itai811.github.io/isoniche/reference/iai.md)).

- weights:

  Pairwise isodars only: `"1/var"` (default, as in the paper) or `NULL`;
  see [`iai()`](https://itai811.github.io/isoniche/reference/iai.md).

- diagnostics:

  Named list of thresholds for the flags; see Details.

- seed:

  Optional random-number seed (bootstrap and null model).

- max_search:

  Largest total searched for the spillover abundance.

## Value

A tibble of class `niche_breadth`, one row per species, with `species`
(if given), `method`, `n_sites`, `nb_di`, `nb_di_lower`, `nb_di_upper`,
`spillover_abundance`, `spillover_lower`, `spillover_upper`, `nb_dd`,
`nb_dd_lower`, `nb_dd_upper`, `dependent_at` (`Inf` for the asymptote),
`var_explained` (line, `"pca"` only), `intercept_p`, `absent` and
`abandoned` (habitats, comma-separated), `reversals` (list of data
frames: `total_abundance`, `habitat_1` (more abundant before),
`habitat_2`, `niche_breadth`), `n_boot_valid`, `flags` (names of the
flags raised, comma-separated), one logical `flag_*` column per
diagnostic, and `note`. Printing it adds one line giving the number of
flagged species.

## Details

**Niche breadth** is the isodar-adjusted inequality
([`iai()`](https://itai811.github.io/isoniche/reference/iai.md)): total
abundance is allocated among habitats as the isodars predict and niche
breadth is 1 - Gini of that allocation (scaled to 0-1 by default).
Absent habitats count as 0 (see
[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md)).

**Spillover abundance** \\T_0\\ is the smallest total abundance (an
integer, in individuals per reference effort; see
[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md))
at which every used habitat holds at least `min_occupancy` individuals:
the point at which the species starts to spill over into its
least-favoured used habitat. DI niche breadth is the value at \\T_0\\.
It depends on the sampling unit (`reference_effort`); DD niche breadth
does not.

**DD niche breadth.** With `dependent_at = "asymptote"` (default) it is
the limit as total abundance \\T \to \infty\\, where the allocation is
proportional to the direction of the isodar line (or follows the
pairwise slopes) and intercepts no longer matter; it is computed
directly. A number gives the value at that total, and `"max_observed"`
the value at the largest site total in the data.

**The curve between them.** While no habitat changes rank or is
abandoned, the allocation is linear in \\T\\, so niche breadth follows
\$\$NB(T) = NB\_{DD} + (NB\_{DI} - NB\_{DD}) T_0 / T,\$\$ and is halfway
between the two values at \\T = 2 T_0\\. Where the predicted abundances
of two habitats cross, the ranking of habitats changes and the curve can
be non-monotonic; these preference-reversal abundances are returned in
`reversals` (totals from \\T_0\\ upward) and flagged. A habitat used at
low density whose predicted abundance falls to 0 at high totals is
listed in `abandoned` (it still counts for DI niche breadth). Use
[`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md)
to compute and plot whole curves.

**Bootstrap.** Whole sites are resampled with replacement and the
isodars refitted, keeping the absent habitats (and the line method) of
the full data. Intervals are percentile intervals.

**Diagnostics.** Each `flag_*` column is `TRUE` when:

- `few_sites`: fewer than `min_sites` sites;

- `few_pair_sites`: some pair of used habitats shares fewer than
  `min_pair_sites` sites;

- `low_abundance`: the median site total is below `min_median_count`
  individuals (counts, not densities), or more than `max_zero_share` of
  the sampled cells of used habitats are zeros;

- `di_extrapolation`: \\T_0\\ lies outside the observed range of site
  totals by more than a factor `extrapolation_ratio`;

- `dd_extrapolation`: for the asymptote, the largest site total is below
  `extrapolation_ratio` x \\T_0\\ (with the default 2, the data do not
  reach the halfway point of the curve); for a numeric `dependent_at`,
  it exceeds the largest site total by more than that factor;

- `abandoned`: some used habitat is abandoned at high density
  (informational);

- `reversal`: habitats change rank along the curve (informational);

- `poor_fit`: the line explains less than `min_var_explained` of the
  variance (`ndim_method = "pca"`), or the pairwise isodars are
  inconsistent: for some three habitats \\i, j, k\\ the slopes \\b\_{ij}
  b\_{jk}\\ and \\b\_{ik}\\ differ in sign or by more than a factor
  `max_slope_ratio`;

- `disconnected`: habitats are not connected through shared sites;

- `intercept_noise`: no isodar intercept differs from those expected
  from sampling alone. Following the paper (section 1.2.3), `n_null`
  data sets are simulated without density-dependent habitat selection:
  each site's individuals are redistributed among its sampled habitats
  in proportion to the species' overall density in each habitat
  (multinomial). Two-sided p-values of the pairwise TLS intercepts are
  Holm-corrected within species; `intercept_p` is the smallest, and the
  flag is set when it exceeds `alpha`. The test needs integer counts
  (`NA` otherwise). Set `n_null = 0` to skip it;

- `unstable_bootstrap`: more than `max_failed_share` of bootstrap
  replicates failed, or a confidence interval is wider than
  `max_ci_width`;

- `undersampled`: a habitat had no records at too few sites to be
  declared absent (see
  [`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md));

- `no_di`: no total up to `max_search` occupies every used habitat;

- `fit_failed`: the isodars could not be fitted (see `note`).

Defaults: `min_sites = 10`, `min_pair_sites = 5`,
`min_median_count = 5`, `max_zero_share = 0.5`,
`extrapolation_ratio = 2`, `min_var_explained = 0.5`,
`max_slope_ratio = 2`, `n_null = 199`, `alpha = 0.05`,
`max_failed_share = 0.1`, `max_ci_width = 0.5`. Change any of them
through `diagnostics`, e.g. `diagnostics = list(min_sites = 8)`.

**Comparing species.** Absent habitats lower niche breadth, so species
compared with each other should be assessed over the same set of
habitats: the habitats available to them.

## See also

[`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md),
[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md),
[`check_coverage()`](https://itai811.github.io/isoniche/reference/check_coverage.md)

## Examples

``` r
three <- moorea_fish[moorea_fish$species %in% c("Acanthurus nigrofuscus",
  "Scarus schlegeli", "Siganus argenteus"), ]
nb <- niche_breadth(three, species = "species", site = "site",
                    habitat = "habitat", abundance = "abundance",
                    n_boot = 99, seed = 1)
nb
#> # A tibble: 3 × 35
#>   species       method n_sites nb_di nb_di_lower nb_di_upper spillover_abundance
#>   <chr>         <chr>    <int> <dbl>       <dbl>       <dbl>               <dbl>
#> 1 Acanthurus n… pairw…      13 0.468      0.0671       0.499                 133
#> 2 Scarus schle… pairw…      13 0.143      0.0462       0.718                  30
#> 3 Siganus arge… pairw…      13 0.159      0.0450       0.561                  42
#> # ℹ 28 more variables: spillover_lower <dbl>, spillover_upper <dbl>,
#> #   nb_dd <dbl>, nb_dd_lower <dbl>, nb_dd_upper <dbl>, dependent_at <dbl>,
#> #   var_explained <dbl>, intercept_p <dbl>, absent <chr>, abandoned <chr>,
#> #   reversals <list>, n_boot_valid <int>, flags <chr>, flag_few_sites <lgl>,
#> #   flag_few_pair_sites <lgl>, flag_low_abundance <lgl>,
#> #   flag_di_extrapolation <lgl>, flag_dd_extrapolation <lgl>,
#> #   flag_abandoned <lgl>, flag_reversal <lgl>, flag_poor_fit <lgl>, …
#> 3 species flagged (not counting abandoned habitats or reversals); see summary()
nb[c("species", "nb_di", "spillover_abundance", "nb_dd", "flags")]
#> # A tibble: 3 × 5
#>   species                nb_di spillover_abundance  nb_dd flags                 
#>   <chr>                  <dbl>               <dbl>  <dbl> <chr>                 
#> 1 Acanthurus nigrofuscus 0.468                 133 0.365  reversal, poor_fit, u…
#> 2 Scarus schlegeli       0.143                  30 0.399  abandoned, reversal, …
#> 3 Siganus argenteus      0.159                  42 0.0102 abandoned, reversal, …
```
