# Niche breadth indices from pairwise isodar parameters

Summarises the intercepts and slopes of the pairwise isodars directly,
as an alternative to the isodar-adjusted inequality
([`iai()`](https://itai811.github.io/isoniche/reference/iai.md)). A
perfect generalist has intercepts of 0 and slopes of 1 between all
habitat pairs.

## Usage

``` r
pairwise_niche(isodars, weights = NULL, log_slope = FALSE)
```

## Arguments

- isodars:

  Pairwise isodars from `fit_isodars(dim = "pairwise")`, or a data frame
  with columns `habitat_x`, `habitat_y`, `intercept`, `slope` (and `sd`
  for `weights = "1/var"`).

- weights:

  `NULL` (equal weights), `"1/var"` (inverse variance of each isodar,
  \\1/(sd^2(1+b^2))\\, computed after orienting the isodars), or a
  numeric vector with one non-negative weight per isodar.

- log_slope:

  If `TRUE`, slope deviations are \\(\log b)^2\\ instead of \\(b -
  1)^2\\.

## Value

A one-row tibble with `mean_intercept`, `mean_slope`,
`slope_niche_breadth`, `intercept_niche_breadth`, `n_isodars` and
`note`. The oriented isodars are in attribute `"isodars"`.

## Details

Each isodar is first expressed so that its intercept is non-negative: an
isodar \\y = a + b x\\ with \\a \< 0\\ and \\b \> 0\\ is rewritten as
\\x = -a/b + (1/b) y\\ (habitats swapped).
[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md)
does this by default (`flip_intercept = TRUE`); `pairwise_niche()`
repeats it so that the result does not depend on that choice.

**Mean indices** are the (weighted) mean intercept and mean slope.

**Inverse root-mean-square deviation indices** measure deviations from
the values that indicate no habitat preference (0 and 1), not from the
mean. The slope niche breadth is \$\$1 / \sqrt{\sum w\_{ij} (b\_{ij} -
1)^2 / \sum w\_{ij}}\$\$ and the intercept niche breadth is \$\$1 /
\sqrt{\sum w\_{ij} a\_{ij}^2 / \sum w\_{ij}}\$\$ where \\a\_{ij}\\ and
\\b\_{ij}\\ are the intercept and slope of the isodar between habitats
\\i\\ and \\j\\, and \\w\_{ij}\\ its weight. Larger values indicate
broader niches; extreme specialists approach 0, and a perfect generalist
gives `Inf`.

Slopes are not symmetric around 1: slopes of 2 and 0.5 describe
preferences of equal strength in opposite directions. `log_slope = TRUE`
measures slope deviations as \\(\log b)^2\\, which removes this
asymmetry. It is undefined for non-positive slopes (for example, a
habitat abandoned as total abundance rises); the slope index is then
`NA` and `note` says why.

**Units.** The intercept index (and the mean intercept) is in units of
abundance per sampling unit, so it depends on survey effort and is
comparable among species only when surveys are standardised. The slope
indices are unitless.

## See also

[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md),
[`iai()`](https://itai811.github.io/isoniche/reference/iai.md)

## Examples

``` r
set.seed(1)
sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
pw <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"))
pairwise_niche(pw)
#> # A tibble: 1 × 6
#>   mean_intercept mean_slope slope_niche_breadth intercept_niche_brea…¹ n_isodars
#>            <dbl>      <dbl>               <dbl>                  <dbl>     <int>
#> 1           7.79       1.75                1.21                  0.116         3
#> # ℹ abbreviated name: ¹​intercept_niche_breadth
#> # ℹ 1 more variable: note <chr>
pairwise_niche(pw, weights = "1/var", log_slope = TRUE)
#> # A tibble: 1 × 6
#>   mean_intercept mean_slope slope_niche_breadth intercept_niche_brea…¹ n_isodars
#>            <dbl>      <dbl>               <dbl>                  <dbl>     <int>
#> 1           6.14       1.59                2.12                  0.151         3
#> # ℹ abbreviated name: ¹​intercept_niche_breadth
#> # ℹ 1 more variable: note <chr>
```
