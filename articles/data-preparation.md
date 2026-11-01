# Preparing your data

``` r

library(isoniche)
```

This vignette covers the decisions to make before estimating niche
breadth: the data format, sampling effort, what counts as a replicate,
which habitats to include, and what to do with incomplete or small
sites.

## Long format

All the main functions
([`check_coverage()`](https://itai811.github.io/isoniche/reference/check_coverage.md),
[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md),
[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md),
[`group_sites()`](https://itai811.github.io/isoniche/reference/group_sites.md))
take data in long format, with one row per species x site x habitat and
columns named by arguments:

| argument | column | required |
|----|----|----|
| `species` | species name | no (without it the data are one species) |
| `site` | replicate (site, or survey) | yes |
| `habitat` | habitat | yes |
| `abundance` | number of individuals | yes |
| `effort` | sampling effort | no (default: 1 for every row) |

Rules:

- A habitat not sampled at a site has no row. A habitat sampled where
  the species was not found has a row with abundance 0. If any species
  has a row for a site and habitat, that site and habitat count as
  sampled, and missing rows of other species there are treated as 0.
- There is at most one row per species, site and habitat; pool repeated
  samples first (see below).
- Sites where a species has no individuals at all are dropped for that
  species: they carry no information on how it distributes among
  habitats.
- If `habitat` is a factor, its levels set the order of habitats.

## From a community matrix

Community data are often stored as a matrix of samples x species, with a
separate table of sample attributes.
[`community_to_long()`](https://itai811.github.io/isoniche/reference/community_to_long.md)
converts them, pooling samples from the same site and habitat and
summing their effort:

``` r

comm <- matrix(c(12, 0, 3,
                  8, 2, 0,
                  5, 1, 1,
                  0, 5, 9,
                  2, 7, 4),
               nrow = 5, byrow = TRUE,
               dimnames = list(NULL, c("sp_a", "sp_b", "sp_c")))
samples <- data.frame(site = c("north", "north", "north", "south", "south"),
                      habitat = c("reef", "reef", "lagoon", "reef", "lagoon"),
                      area_m2 = c(100, 50, 100, 100, 100),
                      lon = c(-149.80, -149.80, -149.80, -149.85, -149.85),
                      lat = c(-17.48, -17.48, -17.48, -17.55, -17.55))
long <- community_to_long(comm, samples, site = "site", habitat = "habitat",
                          effort = "area_m2", keep = c("lon", "lat"))
long
#> # A tibble: 12 × 7
#>    species site  habitat abundance area_m2   lon   lat
#>    <chr>   <chr> <chr>       <dbl>   <dbl> <dbl> <dbl>
#>  1 sp_a    north lagoon          5     100 -150. -17.5
#>  2 sp_a    north reef           20     150 -150. -17.5
#>  3 sp_a    south lagoon          2     100 -150. -17.6
#>  4 sp_a    south reef            0     100 -150. -17.6
#>  5 sp_b    north lagoon          1     100 -150. -17.5
#>  6 sp_b    north reef            2     150 -150. -17.5
#>  7 sp_b    south lagoon          7     100 -150. -17.6
#>  8 sp_b    south reef            5     100 -150. -17.6
#>  9 sp_c    north lagoon          1     100 -150. -17.5
#> 10 sp_c    north reef            3     150 -150. -17.5
#> 11 sp_c    south lagoon          4     100 -150. -17.6
#> 12 sp_c    south reef            9     100 -150. -17.6
```

The result goes straight into the other functions, e.g.
`niche_breadth(long, species = "species", site = "site", habitat = "habitat", abundance = "abundance", effort = "area_m2")`.
Without an `effort` column, each sample has effort 1, so pooled cells
have effort equal to their number of samples.

## Effort and the unit of abundance

Isodars relate *densities*, so counts are divided by effort: density =
abundance / effort. Effort can be area, transect length (with constant
width), trap nights or any relative measure, as long as it is the same
kind of measure in every row.

Densities are then expressed per `reference_effort`, by default the
median effort of the data. This is the sampling unit in which “one
individual” and all reported abundances (spillover abundance, the x-axis
of niche breadth curves) are measured. DI niche breadth is evaluated
where the least favoured habitat holds one individual per reference
effort, so **DI niche breadth and spillover abundance depend on this
unit, while isodar slopes and DD niche breadth do not**. Use the same
`reference_effort` for all species you compare, and report it.

``` r

fish <- moorea_fish[moorea_fish$species == "Chaetodon ornatissimus", ]
fish$effort <- 1
unit1 <- niche_breadth(fish, site = "site", habitat = "habitat",
                       abundance = "abundance", effort = "effort",
                       n_boot = 0, diagnostics = list(n_null = 0))
unit2 <- niche_breadth(fish, site = "site", habitat = "habitat",
                       abundance = "abundance", effort = "effort",
                       reference_effort = 2, n_boot = 0,
                       diagnostics = list(n_null = 0))
data.frame(reference_effort = c(1, 2),
           rbind(unit1[c("nb_di", "spillover_abundance", "nb_dd")],
                 unit2[c("nb_di", "spillover_abundance", "nb_dd")]))
#>   reference_effort     nb_di spillover_abundance     nb_dd
#> 1                1 0.3347681                  27 0.6369444
#> 2                2 0.2742806                  51 0.6369444
```

With `reference_effort = 2`, abundances are counted per two sampling
units, so “one individual” in the least-favoured habitat is a different
point along the isodars: the spillover abundance and DI niche breadth
change (here the spillover abundance almost doubles, but the direction
and size of the change depend on the isodars), while DD niche breadth is
unchanged.

The default `reference_effort`, the median effort of the pooled site x
habitat cells, depends on how many samples happened to be pooled. Set it
to a natural unit instead, for example the area of one standard
transect, and use the same value for every species you compare. The unit
is printed with the results of
[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md).

## Replicates: sites or surveys

Isodars are estimated across replicates that differ in abundance. A
replicate can be a site, or one survey (a year, a season) of a site; put
its identifier in the `site` column. With several sites surveyed
repeatedly, each site x survey can be a replicate (combine the two into
one identifier), but repeated surveys of a site are not independent, and
the bootstrap, which resamples replicates, will then give intervals that
are too narrow. Pooling surveys within a site, or analysing years at a
single site, avoids this.

## Which habitats

Include the habitats *available* to the species: those it could reach. A
habitat the species never uses lowers its niche breadth, so species
compared with each other should be assessed over the same set of
habitats. Rare or poorly sampled habitats can be merged with similar
ones before the analysis.

## Absent and undersampled habitats

A habitat in which a species was never recorded is treated in one of two
ways:

- **absent**, when the habitat was sampled at enough sites where the
  species occurred (`absent_min_sites`, default 5) and the species was
  found there no more than `absent_max_count` times (default 0). It is
  fixed at 0, excluded from the isodar fit and counted as 0 in niche
  breadth.
- **undersampled**, when it was sampled at too few such sites. The zeros
  may then be chance, so the habitat is kept and the species is flagged.

``` r

cov <- check_coverage(moorea_fish, species = "species", site = "site",
                      habitat = "habitat", abundance = "abundance")
cov[cov$species == "Chromis viridis", ]
#> # A tibble: 3 × 8
#>   species     habitat n_sites n_complete n_sampled n_occupied total_count status
#>   <chr>       <chr>     <int>      <int>     <dbl>      <dbl>       <dbl> <chr> 
#> 1 Chromis vi… barrie…      13         13        13         13        4617 used  
#> 2 Chromis vi… fringi…      13         13        13         11        2376 used  
#> 3 Chromis vi… outer_…      13         13        13          0           0 absent
```

The rule can be changed (`absent_max_count`, `absent_min_sites`),
switched off (`absent = "none"`) or replaced by naming the absent
habitats (`absent = "outer_slope"`), in
[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md)
and
[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md).

With the default `absent_max_count = 0`, a single individual recorded in
a habitat makes it used. That one record can move the spillover
abundance a long way (the species must now reach one individual in that
habitat), so if a few stray individuals should not count, raise
`absent_max_count`.

A habitat used at low abundance but predicted to be empty at high
abundance is *not* absent: it counts for DI niche breadth and is listed
in the `abandoned` column of
[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md).

## Incomplete sites

Not every habitat needs to be sampled at every site. Pairwise isodars
(the default, `dim = "pairwise"`) use every site where both habitats of
a pair were sampled. The n-dimensional isodar (`dim = "ndim"`) can use
sites with at least two habitats through a one-factor model
(`ndim_method = "latent"`, needs the **glmmTMB** package); with
`ndim_method = "auto"` it is chosen whenever some sites are incomplete.
Neither is preferred in general; when many sites are incomplete, fitting
both and comparing them is a useful check. `n_sites` in the results
counts the sites where the species was recorded; the latent fit uses
those with at least two sampled habitats.

``` r

set.seed(1)
inc <- simulate_isodars(slopes = c(1.5, 1, 0.8), intercepts = c(5, 0, 2),
                        n = 60, p_missing = 0.3, format = "long")
head(inc)
#> # A tibble: 6 × 3
#>   site  habitat abundance
#>   <chr> <chr>       <dbl>
#> 1 site1 hab1           15
#> 2 site1 hab2           29
#> 3 site1 hab3           28
#> 4 site1 hab4           25
#> 5 site2 hab1           21
#> 6 site2 hab3           38
niche_breadth(inc, site = "site", habitat = "habitat",
              abundance = "abundance", dim = "ndim", n_boot = 0,
              diagnostics = list(n_null = 0))[c("method", "n_sites", "nb_di",
                                                "nb_dd", "flags")]
#> # A tibble: 1 × 5
#>   method        n_sites nb_di nb_dd flags                     
#>   <chr>           <int> <dbl> <dbl> <chr>                     
#> 1 line (latent)      60 0.721 0.887 di_extrapolation, reversal
#> Abundances are per site (equal effort in every site and habitat).
```

## Grouping nearby sites

When single sites lack habitats or hold few individuals, nearby sites
can be merged into larger sampling units with
[`group_sites()`](https://itai811.github.io/isoniche/reference/group_sites.md).
It uses complete linkage on great-circle distance, so every pair of
sites in a group lies within the radius; counts and efforts are summed.

``` r

set.seed(2)
d <- simulate_isodars(1.5, 1.5, 5, 5, n = 30, p_missing = 0.3,
                      format = "long")
xy <- data.frame(site = unique(d$site), lon = stats::runif(30, 0, 0.1),
                 lat = stats::runif(30, 0, 0.1))
d <- merge(d, xy)
g <- group_sites(d, site = "site", lon = "lon", lat = "lat", radius_km = 2,
                 habitat = "habitat", abundance = "abundance")
head(g)
#> # A tibble: 6 × 7
#>   site          habitat abundance effort     lon    lat n_sites
#>   <chr>         <chr>       <dbl>  <dbl>   <dbl>  <dbl>   <int>
#> 1 site1+site26  hab1           28      1 0.0247  0.0673       2
#> 2 site1+site26  hab2           64      2 0.0247  0.0673       2
#> 3 site1+site26  hab3           35      1 0.0247  0.0673       2
#> 4 site10+site15 hab1           31      1 0.00364 0.0263       2
#> 5 site10+site15 hab2           92      2 0.00364 0.0263       2
#> 6 site10+site15 hab3           66      1 0.00364 0.0263       2
```

The radius is not chosen automatically: choosing it to improve the
results would be circular. Instead,
[`radius_sensitivity()`](https://itai811.github.io/isoniche/reference/radius_sensitivity.md)
shows how the number of units, their completeness and the estimates
change with the radius; choose a radius on ecological grounds (e.g. the
scale of movement of the species) and check that the conclusions do not
hinge on it:

``` r

radius_sensitivity(d, site = "site", lon = "lon", lat = "lat",
                   habitat = "habitat", abundance = "abundance",
                   radii_km = c(0, 2, 5), diagnostics = list(n_null = 0))
#> # A tibble: 3 × 7
#>   radius_km n_units complete_share nb_di spillover_abundance nb_dd flags        
#>       <dbl>   <int>          <dbl> <dbl>               <dbl> <dbl> <chr>        
#> 1         0      30          0.233 0.391                  23 0.734 di_extrapola…
#> 2         2      16          0.688 0.384                  36 0.489 reversal, po…
#> 3         5       7          1     0.297                  46 0.695 few_sites, d…
```

`radius_km` can also be a named vector with one radius per species.

## One species in wide format

For a single species, data can also be given in wide format, one row per
site and one column per habitat (`NA` for habitats not sampled), with
the habitat columns named in `hab_cols`:

``` r

set.seed(1)
sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
head(sim)
#> # A tibble: 6 × 4
#>    hab1  hab2  hab3 total_abundance
#>   <dbl> <dbl> <dbl>           <dbl>
#> 1    15    27    45              87
#> 2    21    37    60             118
#> 3    32    54    87             173
#> 4    51    83   129             263
#> 5    11    23    40              74
#> 6    51    82   128             261
niche_breadth(sim, hab_cols = c("hab1", "hab2", "hab3"), n_boot = 0,
              diagnostics = list(n_null = 0))[c("nb_di", "nb_dd")]
#> # A tibble: 1 × 2
#>   nb_di nb_dd
#>   <dbl> <dbl>
#> 1 0.381 0.738
#> Abundances are per site (equal effort in every site and habitat).
```

The lower-level functions
[`iai()`](https://itai811.github.io/isoniche/reference/iai.md) and
[`isoderr()`](https://itai811.github.io/isoniche/reference/isoderr.md)
work on fitted isodars of one species in this format; see their help
pages.
