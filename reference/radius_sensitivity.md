# Sensitivity of niche breadth to the grouping radius

Groups sites with
[`group_sites()`](https://itai811.github.io/isoniche/reference/group_sites.md)
at each radius and estimates niche breadth with
[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md),
showing how the number of sampling units, their completeness and the
estimates change with the radius.

## Usage

``` r
radius_sensitivity(
  data,
  species = NULL,
  site,
  habitat,
  abundance,
  effort = NULL,
  lon,
  lat,
  radii_km = c(0, 0.5, 1, 2, 5),
  linkage = "complete",
  ...
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

- lon, lat:

  Names of the longitude and latitude columns (decimal degrees; one
  location per site).

- radii_km:

  Radii (km) to compare.

- linkage:

  `"complete"` or `"single"`.

- ...:

  Passed to
  [`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md)
  (for example `dim`, `n_boot`; the default here is `n_boot = 0`).

## Value

A tibble with one row per radius x species: `radius_km`, `species` (if
given), `n_units` (sampling units with the species), `complete_share`
(share of those units in which every habitat was sampled), `nb_di`,
`spillover_abundance`, `nb_dd` and `flags`.

## See also

[`group_sites()`](https://itai811.github.io/isoniche/reference/group_sites.md)

## Examples

``` r
set.seed(1)
d <- simulate_isodars(1.5, 1.5, 5, 5, n = 20, p_missing = 0.3,
                      format = "long")
xy <- data.frame(site = unique(d$site), lon = stats::runif(20, 0, 0.1),
                 lat = stats::runif(20, 0, 0.1))
d <- merge(d, xy)
radius_sensitivity(d, site = "site", lon = "lon", lat = "lat",
                   habitat = "habitat", abundance = "abundance",
                   radii_km = c(0, 2, 5), diagnostics = list(n_null = 0))
#> # A tibble: 3 × 7
#>   radius_km n_units complete_share nb_di spillover_abundance nb_dd flags        
#>       <dbl>   <int>          <dbl> <dbl>               <dbl> <dbl> <chr>        
#> 1         0      20          0.35  0.383                  25 0.743 ""           
#> 2         2      14          0.643 0.447                  16 0.736 "di_extrapol…
#> 3         5       6          1     0.390                 160 0.251 "few_sites, …
```
