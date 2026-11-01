# Convert a community matrix to long format

Turns a community matrix (one row per sample, one column per species, as
used by vegan and most community-ecology software) and a table
describing the samples into the long format used by
[`check_coverage()`](https://itai811.github.io/isoniche/reference/check_coverage.md),
[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md)
and
[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md):
one row per species x site x habitat.

## Usage

``` r
community_to_long(comm, samples, site, habitat, effort = NULL, keep = NULL)
```

## Arguments

- comm:

  A matrix or data frame of counts: one row per sample, one column per
  species (column names are the species names). Counts must be
  non-negative and not missing.

- samples:

  A data frame with one row per row of `comm`, in the same order, giving
  the site, habitat and (optionally) effort of each sample.

- site, habitat:

  Names of the site and habitat columns of `samples`.

- effort:

  Optional name of the effort column of `samples` (area, transect
  length, trap nights...; positive).

- keep:

  Optional names of other site-level columns of `samples` to carry over,
  such as longitude and latitude for
  [`group_sites()`](https://itai811.github.io/isoniche/reference/group_sites.md).
  They must have one value per site.

## Value

A tibble with columns `species`, the site and habitat columns (named as
in `samples`), `abundance`, `effort` (named as in `samples`, or
`"effort"` = number of samples) and any `keep` columns, ready for
`niche_breadth(x, species = "species", site = ..., habitat = ..., abundance = "abundance", effort = ...)`.

## Details

Each row of `comm` is one sample (a transect, quadrat, trap, count...)
and the matching row of `samples` says at which site and in which
habitat it was taken. Samples taken at the same site and habitat are
pooled: counts are summed, and so are efforts. Without an `effort`
column, every sample has effort 1, so the effort of a site x habitat is
its number of samples; sites and habitats sampled unevenly are then
still compared as densities per sample.

Every sampled site x habitat appears for every species, with 0 where the
species was not recorded (which is what makes zeros informative).
Species never recorded are dropped.

## See also

[`check_coverage()`](https://itai811.github.io/isoniche/reference/check_coverage.md),
[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md),
[`vignette("data-preparation")`](https://itai811.github.io/isoniche/articles/data-preparation.md)

## Examples

``` r
comm <- matrix(c(5, 0, 2,
                 3, 1, 0,
                 0, 4, 6,
                 1, 2, 2), nrow = 4, byrow = TRUE,
               dimnames = list(NULL, c("sp_a", "sp_b", "sp_c")))
samples <- data.frame(site = c("s1", "s1", "s1", "s2"),
                      habitat = c("reef", "reef", "lagoon", "reef"))
community_to_long(comm, samples, site = "site", habitat = "habitat")
#> # A tibble: 9 × 5
#>   species site  habitat abundance effort
#>   <chr>   <chr> <chr>       <dbl>  <dbl>
#> 1 sp_a    s1    lagoon          0      1
#> 2 sp_a    s1    reef            8      2
#> 3 sp_a    s2    reef            1      1
#> 4 sp_b    s1    lagoon          4      1
#> 5 sp_b    s1    reef            1      2
#> 6 sp_b    s2    reef            2      1
#> 7 sp_c    s1    lagoon          6      1
#> 8 sp_c    s1    reef            2      2
#> 9 sp_c    s2    reef            2      1
```
