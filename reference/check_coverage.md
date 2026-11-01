# Sampling coverage of each species and habitat

Summarises, for each species and habitat, how many of the species' sites
sampled the habitat, how often the species was found there, and whether
[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md)
would treat the habitat as used, absent or undersampled. Run it before
fitting, to decide which species and habitats can be analysed and
whether sites should be grouped
([`group_sites()`](https://itai811.github.io/isoniche/reference/group_sites.md)).

## Usage

``` r
check_coverage(
  data,
  species = NULL,
  site,
  habitat,
  abundance,
  effort = NULL,
  absent_max_count = 0,
  absent_min_sites = 5
)
```

## Arguments

- data:

  Long-format data: one row per species x site x habitat.

- species, site, habitat, abundance, effort:

  Long format: names of the columns holding species (optional; without
  it the data are one species), site (or other replicate), habitat,
  abundance (counts) and sampling effort (optional; without it every row
  has effort 1).

- absent_max_count, absent_min_sites:

  Thresholds for `absent = "auto"`; see Details.

## Value

A tibble with one row per species x habitat: `species` (if given),
`habitat`, `n_sites` (sites where the species was recorded in any
habitat), `n_complete` (of those, sites where every habitat was
sampled), `n_sampled` (sites where this habitat was sampled),
`n_occupied` (sites where the species was recorded in this habitat),
`total_count` and `status` (`"used"`, `"absent"` or `"undersampled"`).

## See also

[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md),
[`group_sites()`](https://itai811.github.io/isoniche/reference/group_sites.md)

## Examples

``` r
cov <- check_coverage(moorea_fish, species = "species", site = "site",
                      habitat = "habitat", abundance = "abundance")
cov[cov$status != "used", ]
#> # A tibble: 325 × 8
#>    species    habitat n_sites n_complete n_sampled n_occupied total_count status
#>    <chr>      <chr>     <int>      <int>     <dbl>      <dbl>       <dbl> <chr> 
#>  1 Abudefduf… outer_…      10         10        10          0           0 absent
#>  2 Abudefduf… outer_…       9          9         9          0           0 absent
#>  3 Acanthuru… barrie…       7          7         7          0           0 absent
#>  4 Acanthuru… fringi…       7          7         7          0           0 absent
#>  5 Acanthuru… barrie…       1          1         1          0           0 under…
#>  6 Acanthuru… outer_…       1          1         1          0           0 under…
#>  7 Acanthuru… barrie…       3          3         3          0           0 under…
#>  8 Acanthuru… barrie…       1          1         1          0           0 under…
#>  9 Acanthuru… fringi…       1          1         1          0           0 under…
#> 10 Acanthuru… barrie…       1          1         1          0           0 under…
#> # ℹ 315 more rows
```
