# Summarise niche breadth results and their diagnostics

Summarises a
[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md)
result: how many species have estimates, the range of niche breadth
values, and for each diagnostic flag how many species raised it, what it
means and what to do about it.

## Usage

``` r
# S3 method for class 'niche_breadth'
summary(object, ...)

# S3 method for class 'summary.niche_breadth'
print(x, ...)
```

## Arguments

- object:

  The output of
  [`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md).

- ...:

  Ignored.

- x:

  The output of [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

An object of class `summary.niche_breadth`: a list with `n_species`,
`n_estimated` (species with DI and DD niche breadth), `n_clean`
(estimated species with no warning flags), `clean` (their names, or row
numbers without a `species` column), `estimates` (quartiles of `nb_di`,
`spillover_abundance` and `nb_dd`), `flags` (a data frame with `flag`,
`type`, `n`, `meaning` and `advice`, for the flags raised) and
`settings`. Flags of type `"info"` (abandoned habitats, preference
reversals) describe the species and are not warnings.

## See also

[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md),
[`vignette("diagnostics")`](https://itai811.github.io/isoniche/articles/diagnostics.md)

## Examples

``` r
three <- moorea_fish[moorea_fish$species %in% c("Acanthurus nigrofuscus",
  "Scarus schlegeli", "Siganus argenteus"), ]
nb <- niche_breadth(three, species = "species", site = "site",
                    habitat = "habitat", abundance = "abundance",
                    n_boot = 0, seed = 1)
summary(nb)
#> Niche breadth for 3 species
#>   DD niche breadth at: asymptote; bootstrap replicates: 0
#>   with DI and DD estimates: 3
#>   without warning flags:    1
#> 
#> Estimates (quartiles):
#>                        min    q25 median    q75     max
#> nb_di                0.143  0.151  0.159  0.314   0.468
#> spillover_abundance 30.000 36.000 42.000 87.500 133.000
#> nb_dd                0.010  0.188  0.365  0.382   0.399
#> 
#> Flags raised (species):
#>   poor_fit               2  a line does not describe the data well, or the pairwise isodars disagree
#>                             -> check the isodars (fit_isodars()); try the other method
#>   reversal               3  habitats change rank along the niche breadth curve
#>                             -> informational: compare species at a stated abundance (see iai_curve())
#>   abandoned              2  a used habitat is abandoned at high density
#>                             -> informational: the habitat counts for DI but not for DD niche breadth
#> 
#> See vignette("diagnostics") for details. Thresholds are set with niche_breadth(diagnostics = ).
```
