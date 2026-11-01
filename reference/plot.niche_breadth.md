# Plot niche breadth against abundance

Computes the niche breadth curves of a
[`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md)
result with
[`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md)
and plots them with
[`plot.iai_curve()`](https://itai811.github.io/isoniche/reference/plot.iai_curve.md).

## Usage

``` r
# S3 method for class 'niche_breadth'
plot(x, n_boot = 0, seed = NULL, ...)
```

## Arguments

- x:

  Output of
  [`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md).

- n_boot, seed:

  Passed to
  [`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md)
  (bootstrap bands; default none).

- ...:

  Passed to
  [`plot.iai_curve()`](https://itai811.github.io/isoniche/reference/plot.iai_curve.md),
  e.g. `species` or `x_axis`.

## Value

The
[`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md)
result, invisibly.

## See also

[`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md),
[`plot.iai_curve()`](https://itai811.github.io/isoniche/reference/plot.iai_curve.md)

## Examples

``` r
three <- moorea_fish[moorea_fish$species %in% c("Acanthurus nigrofuscus",
  "Scarus schlegeli", "Siganus argenteus"), ]
nb <- niche_breadth(three, species = "species", site = "site",
                    habitat = "habitat", abundance = "abundance",
                    n_boot = 0, diagnostics = list(n_null = 0))
plot(nb)
```
