# Isodar-adjusted inequality (IAI)

Computes niche breadth from fitted isodars by allocating a given total
abundance among habitats according to the isodars and measuring the
inequality of that allocation (1 - Gini). Evaluating it at several
totals gives a rarefaction-like curve of niche breadth against
abundance.

## Usage

``` r
iai(
  isodars,
  abundances = NULL,
  weights = NULL,
  scaled = TRUE,
  min_occupancy = 1,
  max_search = 10000,
  dim = NULL,
  plot = FALSE,
  method = "gini"
)
```

## Arguments

- isodars:

  Output of
  [`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md).

- abundances:

  Numeric vector of total abundances, or `NULL` for the
  density-independent value.

- weights:

  Pairwise only: `NULL` (equal weights), `"1/var"` (inverse variance of
  each isodar's vertical residuals, \\1/(sd^2(1+b^2))\\), or a numeric
  vector with one weight per isodar.

- scaled:

  If `TRUE`, 1 - Gini is rescaled to 0-1, so that a species found in a
  single habitat scores 0 and an even spread scores 1.

- min_occupancy:

  Individuals required in every habitat for the density-independent
  value.

- max_search:

  Largest total searched for the density-independent value.

- dim:

  Optional; inferred from `isodars`. Kept for compatibility.

- plot:

  If `TRUE` and several abundances are given, plots the curve.

- method:

  Inequality measure; only `"gini"` is currently implemented.

## Value

A tibble with `total_abundance`, `niche_breadth`, one column per habitat
(predicted abundance), and `note`.

## Details

**Allocation.** For pairwise isodars, the predicted abundances are those
that minimise the (optionally weighted) sum of squared deviations from
all pairwise isodar equations, \\y - b x = a\\, subject to summing to
the total. This constrained least-squares problem is solved directly;
habitats predicted to be negative are set to zero and the rest are
solved again. For an n-dimensional isodar (a line), the allocation is
the single point at which the line reaches the total; negative habitats
are set to zero.

**Density-independent niche breadth** (`abundances = NULL`) is evaluated
at the smallest total (searched over integers) at which every used
habitat holds at least `min_occupancy` individuals. If no total up to
`max_search` achieves this, a single row with `NA` is returned.

**Absent habitats** (see
[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md))
are fixed at 0 and count as 0 in the inequality, so they lower niche
breadth; with a single used habitat, niche breadth is 0.

**Density-dependent niche breadth** is the value at a high total
abundance, supplied through `abundances`.

## See also

[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md),
[`isoderr()`](https://itai811.github.io/isoniche/reference/isoderr.md),
[`pairwise_niche()`](https://itai811.github.io/isoniche/reference/pairwise_niche.md)

## Examples

``` r
set.seed(1)
sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
pw <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"))
iai(pw)                          # density-independent
#> # A tibble: 1 × 6
#>   total_abundance niche_breadth  hab1  hab2  hab3 note 
#>             <dbl>         <dbl> <dbl> <dbl> <dbl> <chr>
#> 1              23         0.381  1.03  6.70  15.3 ""   
iai(pw, abundances = c(50, 500)) # density-dependent
#> # A tibble: 2 × 6
#>   total_abundance niche_breadth   hab1  hab2  hab3 note 
#>             <dbl>         <dbl>  <dbl> <dbl> <dbl> <chr>
#> 1              50         0.574   6.72  15.2  28.0 ""   
#> 2             500         0.721 102.   157.  241.  ""   
ln <- fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"), dim = "ndim")
iai(ln, abundances = 500)
#> # A tibble: 1 × 6
#>   total_abundance niche_breadth  hab1  hab2  hab3 note 
#>             <dbl>         <dbl> <dbl> <dbl> <dbl> <chr>
#> 1             500         0.721  102.  157.  241. ""   
```
