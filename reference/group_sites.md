# Group nearby sites into sampling units

Merges sites within a given distance of each other into one sampling
unit, summing counts and efforts. Grouping helps when single sites lack
some habitats or hold too few individuals, at the cost of fewer units.

## Usage

``` r
group_sites(
  data,
  site,
  lon,
  lat,
  radius_km,
  linkage = "complete",
  species = NULL,
  habitat,
  abundance,
  effort = NULL
)
```

## Arguments

- data:

  Long-format data: one row per species x site x habitat.

- lon, lat:

  Names of the longitude and latitude columns (decimal degrees; one
  location per site).

- radius_km:

  Grouping radius in km: one number, or a named vector with one radius
  per species.

- linkage:

  `"complete"` or `"single"`.

- species, site, habitat, abundance, effort:

  Long format: names of the columns holding species (optional; without
  it the data are one species), site (or other replicate), habitat,
  abundance (counts) and sampling effort (optional; without it every row
  has effort 1).

## Value

A tibble in long format with the original column names: the species (if
given), site (the group name), habitat, abundance and effort (named
`"effort"` if there was no effort column), plus `lon` and `lat` (mean of
the sites) and `n_sites` (number of sites merged). Every sampled group x
habitat appears for every species, with 0 where the species was not
recorded.

## Details

Sites are clustered by complete linkage on great-circle distance, so
every pair of sites in a group lies within `radius_km` of each other
(single linkage, `linkage = "single"`, can chain distant sites
together). With `radius_km = 0` only sites with identical coordinates
are merged.

Counts are summed over the sites of a group, and so are efforts: the
effort of a group in a habitat is the summed effort of the sites that
sampled that habitat. Without an `effort` column, every site x habitat
has effort 1, so the group's effort is the number of sites that sampled
the habitat. Groups are named after their sites (joined with `"+"`).

The radius is not chosen automatically: choosing it to maximise
completeness or niche breadth would be circular. Use
[`radius_sensitivity()`](https://itai811.github.io/isoniche/reference/radius_sensitivity.md)
to see how the number of units, their completeness and the niche breadth
estimates change with the radius.

## See also

[`radius_sensitivity()`](https://itai811.github.io/isoniche/reference/radius_sensitivity.md),
[`check_coverage()`](https://itai811.github.io/isoniche/reference/check_coverage.md)

## Examples

``` r
set.seed(1)
d <- simulate_isodars(1.5, 1.5, 5, 5, n = 12, format = "long")
sites <- unique(d$site)
xy <- data.frame(site = sites, lon = stats::runif(12, 0, 0.05),
                 lat = stats::runif(12, 0, 0.05))
d <- merge(d[c("site", "habitat", "abundance")], xy)
g <- group_sites(d, site = "site", lon = "lon", lat = "lat", radius_km = 2,
                 habitat = "habitat", abundance = "abundance")
head(g)
#> # A tibble: 6 × 7
#>   site               habitat abundance effort    lon    lat n_sites
#>   <chr>              <chr>       <dbl>  <dbl>  <dbl>  <dbl>   <int>
#> 1 site1+site12       hab1           25      2 0.0438 0.0168       2
#> 2 site1+site12       hab2           49      2 0.0438 0.0168       2
#> 3 site1+site12       hab3           85      2 0.0438 0.0168       2
#> 4 site10+site5+site8 hab1           51      3 0.0382 0.0423       3
#> 5 site10+site5+site8 hab2           93      3 0.0382 0.0423       3
#> 6 site10+site5+site8 hab3          152      3 0.0382 0.0423       3
```
