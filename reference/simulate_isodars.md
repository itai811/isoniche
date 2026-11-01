# Simulate abundances in habitats linked by isodars

Simulates abundances at `n` sites along a chain of isodars: habitat 1 is
drawn uniformly between 0 and `hab1_max`, and each following habitat
\\k + 1\\ = `intercepts[k]` + `slopes[k]` x habitat \\k\\ + noise.
Values are rounded to integers and negative values set to 0. The default
is three habitats, set by `slope1`, `slope2`, `int1` and `int2`; the
sequence of random numbers is the same as in version 1.0.0.

## Usage

``` r
simulate_isodars(
  slope1 = 1,
  slope2 = 1,
  int1 = 0,
  int2 = 0,
  noise = 1,
  n = 30,
  hab1_max = 100/(slopes[1]^1.2) - intercepts[1],
  slopes = c(slope1, slope2),
  intercepts = c(int1, int2),
  p_missing = 0,
  min_habitats = 2,
  format = c("wide", "long")
)
```

## Arguments

- slope1, int1:

  Isodar of habitat 2 on habitat 1.

- slope2, int2:

  Isodar of habitat 3 on habitat 2.

- noise:

  Standard deviation of the Gaussian noise.

- n:

  Number of sites.

- hab1_max:

  Upper bound of abundance in habitat 1. The default,
  `100 / slopes[1]^1.2 - intercepts[1]`, keeps total abundances in a
  comparable range across parameter values.

- slopes, intercepts:

  Any number of habitats: vectors of length \\H - 1\\ giving the isodar
  of each habitat on the previous one (default: `c(slope1, slope2)` and
  `c(int1, int2)`).

- p_missing:

  Probability that a habitat is missing (not sampled) at a site; at
  least `min_habitats` habitats are kept at every site.

- min_habitats:

  Minimum number of habitats recorded per site.

- format:

  `"wide"` (one column per habitat, `NA` where missing) or `"long"`
  (columns `site`, `habitat`, `abundance`; missing habitats have no
  row).

## Value

Wide: a tibble with columns `hab1`, ..., `habH` and `total_abundance`.
Long: a tibble with `site`, `habitat` and `abundance`.

## Examples

``` r
set.seed(1)
simulate_isodars(1, 1.1, 1, 0, noise = 2, n = 10)
#> # A tibble: 10 × 4
#>     hab1  hab2  hab3 total_abundance
#>    <dbl> <dbl> <dbl>           <dbl>
#>  1    26    26    28              80
#>  2    37    39    43             119
#>  3    57    59    67             183
#>  4    90    92   103             285
#>  5    20    20    24              64
#>  6    89    93   104             286
#>  7    94    95   106             295
#>  8    65    65    72             202
#>  9    62    59    61             182
#> 10     6     9    12              27
simulate_isodars(slopes = c(1.5, 1, 0.8), intercepts = c(5, 0, 2), n = 5,
                 p_missing = 0.3)
#> # A tibble: 5 × 5
#>    hab1  hab2  hab3  hab4 total_abundance
#>   <dbl> <dbl> <dbl> <dbl>           <dbl>
#> 1    NA    NA    44    38              82
#> 2    49    NA    79    65             193
#> 3    NA    42    43    NA              85
#> 4    NA    25    26    24              75
#> 5    NA    11    10    10              31
```
