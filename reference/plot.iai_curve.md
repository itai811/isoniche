# Plot niche breadth curves

Plots the curves from
[`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md):
one line per species with bootstrap bands (dotted), a point at each
species' density-independent value (at its own spillover abundance), the
asymptote (density-dependent value) as a short horizontal dash at the
right, and preference reversals (crosses).

## Usage

``` r
# S3 method for class 'iai_curve'
plot(
  x,
  species = NULL,
  x_axis = c("log_abundance", "relative"),
  abundances = NULL,
  crossings = FALSE,
  legend = TRUE,
  col = NULL,
  ...
)
```

## Arguments

- x:

  Output of
  [`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md).

- species:

  Optional species to plot (default: all).

- x_axis:

  `"log_abundance"` (total abundance, log scale) or `"relative"` (total
  abundance / spillover abundance, log scale; every species' DI value is
  then at 1).

- abundances:

  Optional total abundances to mark with vertical lines
  (`"log_abundance"` only).

- crossings:

  If `TRUE`, marks the totals where two species' curves cross (open
  circles; `"log_abundance"` only).

- legend:

  `TRUE` (default) adds a legend (up to 12 species) in a corner free of
  curves and bands, or above the plot if there is none; a position such
  as `"topleft"` or `"bottomright"` (see
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html)), or
  `"outside"` (above the plot), places it there; `FALSE` omits it.

- col:

  Colours, one per species (default: the current palette).

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.

## See also

[`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md)
