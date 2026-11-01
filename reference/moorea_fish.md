# Reef fishes of Moorea in three habitats

Counts of 325 reef fish species in three reef habitats at 13 sites
around Moorea, French Polynesia: the data of Granot, Dubiner & Belmaker
(Table S2). Fish were counted on belt transects by CRIOBE between 2004
and 2018; there is one count per species, site and habitat. Only sites
where a species was recorded are included for that species; all three
habitats were surveyed at every site, and zeros are explicit.

## Usage

``` r
moorea_fish
```

## Format

A tibble with one row per species x site x habitat:

- species:

  Species name.

- site:

  Site (13 sites).

- habitat:

  `"barrier_reef"`, `"fringing_reef"` or `"outer_slope"`.

- abundance:

  Number of individuals counted.

## Source

"MCR LTER Coral Reef Long-term Population and Community Dynamics:
Fishes" in BioTIME 2.0: Dornelas, M., Antão, L. H., Bates, A. E. et al.
(2025). BioTIME 2.0: expanding and improving a database of biodiversity
time series. *Global Ecology and Biogeography*, 34(5), e70003.

## Examples

``` r
head(moorea_fish)
#> # A tibble: 6 × 4
#>   species                   site      habitat       abundance
#>   <chr>                     <chr>     <chr>             <int>
#> 1 Abudefduf septemfasciatus Afareaitu barrier_reef          0
#> 2 Abudefduf septemfasciatus Afareaitu fringing_reef         1
#> 3 Abudefduf septemfasciatus Afareaitu outer_slope           0
#> 4 Abudefduf septemfasciatus Aroa      barrier_reef          8
#> 5 Abudefduf septemfasciatus Aroa      fringing_reef         7
#> 6 Abudefduf septemfasciatus Aroa      outer_slope           0
length(unique(moorea_fish$species))
#> [1] 325
```
