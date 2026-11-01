# Fit isodars

Fits isodars to abundances (or densities) recorded in two or more
habitats at a set of sites (or other replicates, such as years). Data
are either in wide format (one species; one row per site, one column per
habitat) or in long format (one or many species; one row per species x
site x habitat).

## Usage

``` r
fit_isodars(
  data,
  hab_cols = NULL,
  n_habitats = NULL,
  dim = c("pairwise", "ndim"),
  flip_intercept = TRUE,
  ndim_method = c("auto", "pca", "latent"),
  absent = "auto",
  absent_max_count = 0,
  absent_min_sites = 5,
  species = NULL,
  site = NULL,
  habitat = NULL,
  abundance = NULL,
  effort = NULL,
  reference_effort = NULL
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

- flip_intercept:

  Pairwise only. If `TRUE`, an isodar with a negative intercept and
  positive slope is expressed in the reverse direction (habitats
  swapped), so that intercepts are non-negative.

- ndim_method:

  `"auto"`, `"pca"` or `"latent"`; see Details.

- absent:

  Habitats the species does not use: `"auto"` (detected from the data,
  see Details), `"none"`, or a character vector of habitat names
  (`character(0)` is the same as `"none"`).

- absent_max_count, absent_min_sites:

  Thresholds for `absent = "auto"`; see Details.

- species, site, habitat, abundance, effort:

  Long format: names of the columns holding species (optional; without
  it the data are one species), site (or other replicate), habitat,
  abundance (counts) and sampling effort (optional; without it every row
  has effort 1).

- reference_effort:

  Long format: the sampling unit in which densities are expressed
  (default: median effort); see Details.

## Value

For `dim = "pairwise"`, a tibble of class `isodars_pairwise` with one
row per pair of used habitats: `habitat_x`, `habitat_y`, `intercept`,
`slope`, `sd` (SD of perpendicular residuals), `n_sites`, and `flipped`
(whether the habitats were swapped). For `dim = "ndim"`, a tibble of
class `isodar_line` with one row per habitat: `habitat`, `mean`,
`direction` (both 0 for absent habitats); attributes `var_explained`
(share of variance along the line, `"pca"` only), `eigenvalues`,
`n_sites` and `ndim_method`. Both have attributes `habitats` (all
habitats), `absent` and `undersampled`.

In long format, a tibble of class `isodar_fits` with one row per
species: `species` (if given), `fit` (the fit as above, or `NULL` if it
failed), `counts` and `effort` (site x habitat matrices), `n_sites` and
`note` (warnings and errors from fitting). One summary warning is given
if any species could not be fitted.

## Details

**Input formats.** In wide format, give `hab_cols` (or `n_habitats`);
the values are used as they are, and the result is a single fit. In long
format, give the column names `site`, `habitat` and `abundance` (and
optionally `species` and `effort`); the result is an `isodar_fits` table
with one row per species, which
[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md)
and
[`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md)
take. Replicates may be sites or time steps (for example, years at one
site); in long format they go in the `site` column. A site x habitat
combination that appears for any species is taken as sampled, and a
species without a row there has 0 individuals there; a habitat missing
at a site simply has no rows. Sites where a species has no individuals
in any habitat carry no information on its habitat selection and are
dropped for that species.

**Effort and densities.** Isodars relate densities. In long format,
density = abundance / effort x `reference_effort`. Effort may be an
area, a transect length (constant width) or a unitless relative effort,
but must be the same kind of measure in every row; it must be the same
for every row of a site x habitat combination. `reference_effort`
(default: the median effort of sampled site x habitat combinations) is
the sampling unit in which abundances are expressed: "one individual" in
the density-independent niche breadth, the spillover abundance and the
x-axis of
[`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md)
all mean one individual per `reference_effort`. Density-independent
niche breadth and spillover abundance therefore depend on this unit;
slopes and density-dependent niche breadth do not. The default depends
on how samples happen to be pooled, so set `reference_effort` to a
natural unit (for example, the area of one standard transect) and use
the same value for every species you compare.

In wide format, rows are used as given: unlike long format, a site with
0 in every habitat is kept (it pulls the isodars toward the origin), so
remove such rows if they mean the species was absent from the site.

**Pairwise isodars** (`dim = "pairwise"`) are fitted for every pair of
habitats by total least squares (TLS), which is identical to major-axis
regression in two dimensions. Each pair uses every site where both
habitats were recorded (non-missing), so sites lacking some habitats
still contribute. When many sites are incomplete, both pairwise isodars
and the latent line (below) use all of them; fitting both is a useful
check.

**The n-dimensional isodar** (`dim = "ndim"`) is a line. Under the ideal
free distribution, equal fitness across \\n\\ habitats imposes \\n-1\\
equalities, so abundances across habitats change along a single line as
total abundance changes, whatever the number of habitats. The line is
described by a point (`mean`) and a unit direction vector (`direction`),
with one component per habitat. Each component is the rate at which that
habitat gains individuals as total abundance rises (the analogue of the
isodar slope); where the line enters each habitat reflects its
density-independent preference (the analogue of the intercept). Every
pairwise isodar can be derived from it, and with two habitats it is the
classic isodar.

With `ndim_method = "pca"` the line is the first principal axis of the
covariance matrix of complete sites (the TLS line). With
`ndim_method = "latent"` it is estimated from all sites with at least
two recorded habitats, as a one-factor model, \\N\_{hs} = \mu_h + v_h
t_s + \varepsilon\_{hs}\\, where \\t_s\\ is the position of site \\s\\
along the isodar (requires glmmTMB). `ndim_method = "auto"` uses `"pca"`
when every site is complete and `"latent"` otherwise.

**Absent habitats.** A habitat that the species never uses carries no
information about habitat selection, and an isodar cannot be fitted to
it. Such a habitat is excluded from the fit and fixed at 0 abundance;
[`iai()`](https://itai811.github.io/isoniche/reference/iai.md) still
counts it, as 0, in the inequality, so it lowers niche breadth. With
`absent = "auto"` a habitat is declared absent if the species was
recorded in it no more than `absent_max_count` times in total (in the
units of wide `data`; as counts, before conversion to densities, in long
format), and the habitat was sampled at no fewer than `absent_min_sites`
sites where the species was present in another habitat. A habitat with
(almost) no records but fewer such sites is *undersampled*: it is kept
in the fit, with a warning, because its zeros may reflect too little
sampling rather than avoidance. A habitat used at low density and
abandoned at high density is not absent; it is fitted, and
[`iai()`](https://itai811.github.io/isoniche/reference/iai.md) clips it
to 0 at totals where the isodars predict negative abundance. With the
default `absent_max_count = 0`, a single individual makes a habitat
used, which can change the spillover abundance a lot; raise
`absent_max_count` if a few stray individuals should not count.

Absent habitats lower niche breadth, so species should be compared over
the same set of habitats, namely the habitats available to them.

## References

Morris, D.W. (1987). Spatial scale and the cost of density-dependent
habitat selection. *Evolutionary Ecology*, 1, 379–388.

## See also

[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md),
[`iai()`](https://itai811.github.io/isoniche/reference/iai.md),
[`isoderr()`](https://itai811.github.io/isoniche/reference/isoderr.md)

## Examples

``` r
set.seed(1)
sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"))
#> # A tibble: 3 × 7
#>   habitat_x habitat_y intercept slope    sd n_sites flipped
#>   <chr>     <chr>         <dbl> <dbl> <dbl>   <int> <lgl>  
#> 1 hab1      hab2           5.18  1.50 0.649      30 FALSE  
#> 2 hab1      hab3          13.0   2.24 0.958      30 FALSE  
#> 3 hab2      hab3           5.24  1.50 0.723      30 FALSE  
fit_isodars(sim, hab_cols = c("hab1", "hab2", "hab3"), dim = "ndim")
#> n-dimensional isodar (line), method: pca | sites: 30 
#> Variance along the line: 0.999 
#> # A tibble: 3 × 3
#>   habitat  mean direction
#>   <chr>   <dbl>     <dbl>
#> 1 hab1     28.8     0.348
#> 2 hab2     48.3     0.521
#> 3 hab3     77.5     0.780

# long format, several species
fish <- moorea_fish[moorea_fish$species %in% c("Acanthurus nigrofuscus",
  "Scarus schlegeli", "Zebrasoma scopas"), ]
fits <- fit_isodars(fish, species = "species", site = "site",
                    habitat = "habitat", abundance = "abundance")
fits
#> Isodars (pairwise) for 3 species; habitats: barrier_reef, fringing_reef, outer_slope; reference effort: 1
#> # A tibble: 3 × 6
#>   species                fit                counts   effort   n_sites note 
#>   <chr>                  <list>             <list>   <list>     <int> <chr>
#> 1 Acanthurus nigrofuscus <isdrs_pr [3 × 7]> <dbl[…]> <dbl[…]>      13 ""   
#> 2 Scarus schlegeli       <isdrs_pr [3 × 7]> <dbl[…]> <dbl[…]>      13 ""   
#> 3 Zebrasoma scopas       <isdrs_pr [3 × 7]> <dbl[…]> <dbl[…]>      13 ""   
fits$fit[[which(fits$species == "Scarus schlegeli")]]
#> # A tibble: 3 × 7
#>   habitat_x     habitat_y    intercept    slope    sd n_sites flipped
#>   <chr>         <chr>            <dbl>    <dbl> <dbl>   <int> <lgl>  
#> 1 fringing_reef barrier_reef     27.0   0.654   52.3       13 TRUE   
#> 2 barrier_reef  outer_slope       1.78 -0.00259  1.37      13 FALSE  
#> 3 fringing_reef outer_slope       1.94 -0.00330  1.34      13 FALSE  
```
