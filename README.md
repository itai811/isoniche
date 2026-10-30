# isoniche

<!-- badges: start -->
[![R-CMD-check](https://github.com/itai811/isoniche/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/itai811/isoniche/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/isoniche)](https://CRAN.R-project.org/package=isoniche)
<!-- badges: end -->

Abundance-based indices of niche breadth are biased: a rare species looks
like a specialist simply because it has too few individuals to fill every
habitat. **isoniche** corrects habitat niche breadth for abundance using
*isodars*, the relationships between a species' densities in pairs of
habitats under density-dependent habitat selection.

The package fits pairwise and n-dimensional isodars (including data in which
some sites lack some habitats) and computes the isodar-adjusted inequality
(IAI): niche breadth (1 - Gini) of the allocation of a given total abundance
among habitats predicted by the isodars. This gives

* **density-independent niche breadth**, at the smallest total that places at
  least one individual in every habitat, and
* **density-dependent niche breadth**, at high total abundance,

with bootstrap uncertainty.

## Installation

```r
# development version (1.1.0)
# install.packages("remotes")
remotes::install_github("itai811/isoniche")

# CRAN (1.0.0; see NEWS.md for the fixes in 1.1.0)
install.packages("isoniche")
```

## Example

```r
library(isoniche)

set.seed(1)
sim <- simulate_isodars(1.5, 1.5, 5, 5, noise = 1, n = 30)
habs <- c("hab1", "hab2", "hab3")

# pairwise isodars (total least squares)
pw <- fit_isodars(sim, hab_cols = habs)
iai(pw, weights = "1/var")                       # density-independent
iai(pw, abundances = 1000, weights = "1/var")    # density-dependent

# n-dimensional isodar: a line through habitat space
ln <- fit_isodars(sim, hab_cols = habs, dim = "ndim")
iai(ln)
iai(ln, abundances = 1000)

# bootstrap over sites
isoderr(sim, hab_cols = habs, abundances = 1000, n_boot = 200,
        output = "CI", seed = 1)
```

## Version 1.1.0

Version 1.1.0 rewrites the core and fixes several bugs in 1.0.0, including the
n-dimensional isodar, which is now a line. Results can differ from 1.0.0; see
[NEWS.md](NEWS.md).

## Citation

Granot, I., Dubiner, S. & Belmaker, J. Why abundance-based indices of niche
breadth are biased, and what can be done to improve them. *Ecology Letters*
(in revision). Itai Granot and Shahar Dubiner contributed equally.

## License

GPL-3
