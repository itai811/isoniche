# Niche breadth indices from pairwise isodar parameters

Summarises the intercepts and slopes of the pairwise isodars directly
(Appendix S1 of Granot, Dubiner & Belmaker), as an alternative to the
isodar-adjusted inequality
([`iai()`](https://itai811.github.io/isoniche/reference/iai.md)), and
checks that the isodars are consistent with each other. A perfect
generalist has intercepts of 0 and slopes of 1 between all habitat
pairs.

## Usage

``` r
pairwise_niche(
  isodars,
  weights = NULL,
  scale = c("log", "linear"),
  cycle_tolerance = log(1.5),
  log_slope = NULL
)
```

## Arguments

- isodars:

  Pairwise isodars from `fit_isodars(dim = "pairwise")`, or a data frame
  with columns `habitat_x`, `habitat_y`, `intercept`, `slope` (and `sd`,
  the SD of perpendicular residuals, for `weights = "1/var"`).

- weights:

  `NULL` (equal weights), `"1/var"` (inverse variance of each isodar;
  see Details), or a numeric vector with one non-negative weight per
  isodar.

- scale:

  `"log"` (default; orientation-invariant indices) or `"linear"` (the
  earlier indices); see Details.

- cycle_tolerance:

  Largest acceptable \\\|\log(product)\|\\ of slopes around a cycle of
  three habitats (default `log(1.5)`).

- log_slope:

  Deprecated: use `scale`. `TRUE` is `scale = "log"` and `FALSE` is
  `scale = "linear"`.

## Value

An object of class `pairwise_niche`, a list with

- `indices`: a one-row tibble. For `scale = "log"`:
  `mean_abs_log_slope`, `mean_distance`, `slope_niche_breadth`,
  `intercept_niche_breadth`, `n_isodars`, `n_slope` (isodars used for
  the slope indices) and `note`. For `scale = "linear"`:
  `mean_intercept`, `mean_slope`, `slope_niche_breadth`,
  `intercept_niche_breadth`, `n_isodars` and `note`.

- `cycles`: a tibble with one row per cycle of three habitats with
  isodars between all three pairs: `habitat_1`, `habitat_2`,
  `habitat_3`, `product` (\\b\_{12} b\_{23} b\_{31}\\), `log_product`
  and `inconsistent`.

- `inconsistent`: `TRUE` if any cycle is inconsistent.

- `isodars`: the isodars used, with their `weight`, and for
  `scale = "log"` the `distance` \\d\\ and `log_slope`; for
  `scale = "linear"` they are oriented to non-negative intercepts.

- `scale`, `cycle_tolerance`.

## Details

**Orientation invariance.** An isodar between habitats \\i\\ and \\j\\
can be written with either habitat on the x-axis: \\N_j = a + b N_i\\
or, equivalently, \\N_i = -a/b + (1/b) N_j\\. With `scale = "log"`
(default) every index is the same whichever way each isodar is written:

- slope deviations are \\(\log b)^2\\, since \\\log(1/b) = -\log b\\;

- intercept deviations use the perpendicular distance of the isodar from
  the origin, \\d = \|a\| / \sqrt{1 + b^2}\\, which is the same for both
  ways of writing the line.

The variance-based indices are \$\$1 / \sqrt{\sum w (\log b)^2 / \sum
w}\$\$ for slopes and \$\$1 / \sqrt{\sum w d^2 / \sum w}\$\$ for
intercepts, where \\w\\ are the isodar weights. They measure deviations
from the values that indicate no habitat preference (\\b = 1\\, \\d =
0\\). Larger values indicate broader niches; extreme specialists
approach 0, and a perfect generalist gives `Inf`. The mean indices are
the weighted means of \\\|\log b\|\\ and of \\d\\.

With `weights = "1/var"` and `scale = "log"`, each isodar is weighted by
\\1/sd^2\\, the inverse variance of its perpendicular residuals, which
is also orientation-invariant.

**Non-positive slopes** (for example, a habitat that loses individuals
as another gains them) have no logarithm. Such isodars are left out of
the slope indices, with a warning naming them, but still used in the
intercept indices.

**Linear scale.** `scale = "linear"` gives the earlier indices: each
isodar is first expressed so that its intercept is non-negative
(habitats swapped if needed), then the mean intercept and mean slope,
and \\1 / \sqrt{\sum w (b - 1)^2 / \sum w}\\ and \\1 / \sqrt{\sum w a^2
/ \sum w}\\, with `"1/var"` weights \\1/(sd^2(1+b^2))\\. These depend on
how each isodar is oriented, and slopes of 2 and 0.5 (preferences of
equal strength in opposite directions) count as different deviations.

**Consistency check.** If the isodars describe one consistent pattern of
habitat selection, the slopes around any cycle of three habitats
multiply to 1: \\b\_{AB} b\_{BC} b\_{CA} = 1\\, where \\b\_{AB}\\ is the
slope of habitat \\B\\ on habitat \\A\\. Each cycle's product is
returned in `cycles`, and the result is flagged as inconsistent when
\\\|\log(product)\|\\ exceeds `cycle_tolerance` for any cycle (or the
product is not positive).

**Units.** The intercept indices are in units of abundance per sampling
unit (see `reference_effort` in
[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md)),
so they depend on survey effort and are comparable among species only
when surveys are standardised. The slope indices are unitless.

## See also

[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md),
[`iai()`](https://itai811.github.io/isoniche/reference/iai.md)

## Examples

``` r
set.seed(1)
sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
pw <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"))
pn <- pairwise_niche(pw, weights = "1/var")
pn
#> Niche breadth from pairwise isodars (log scale)
#> # A tibble: 1 × 7
#>   mean_abs_log_slope mean_distance slope_niche_breadth intercept_niche_breadth
#>                <dbl>         <dbl>               <dbl>                   <dbl>
#> 1              0.486          3.37                1.95                   0.285
#> # ℹ 3 more variables: n_isodars <int>, n_slope <int>, note <chr>
#> 1 cycle(s) of three habitats; 0 inconsistent (|log(product of slopes)| > 0.405); see $cycles
pn$cycles
#> # A tibble: 1 × 6
#>   habitat_1 habitat_2 habitat_3 product log_product inconsistent
#>   <chr>     <chr>     <chr>       <dbl>       <dbl> <lgl>       
#> 1 hab1      hab2      hab3        0.999   -0.000649 FALSE       
pairwise_niche(pw, scale = "linear")$indices
#> # A tibble: 1 × 6
#>   mean_intercept mean_slope slope_niche_breadth intercept_niche_brea…¹ n_isodars
#>            <dbl>      <dbl>               <dbl>                  <dbl>     <int>
#> 1           7.79       1.75                1.21                  0.116         3
#> # ℹ abbreviated name: ¹​intercept_niche_breadth
#> # ℹ 1 more variable: note <chr>
```
