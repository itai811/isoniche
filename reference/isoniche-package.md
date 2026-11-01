# isoniche: abundance-independent habitat niche breadth

Estimates the habitat niche breadth of species from counts in several
habitats, corrected for abundance with isodars (Granot, Dubiner &
Belmaker 2026). For each species it reports density-independent (DI)
niche breadth, how widely the species uses habitats when rare;
density-dependent (DD) niche breadth, how widely it uses them when
abundant; and the curve of niche breadth against abundance between the
two.

## A typical analysis

1.  Put the data in long format, one row per species x site x habitat
    ([`community_to_long()`](https://itai811.github.io/isoniche/reference/community_to_long.md)
    converts a community matrix).

2.  Check sampling coverage with
    [`check_coverage()`](https://itai811.github.io/isoniche/reference/check_coverage.md);
    group nearby sites with
    [`group_sites()`](https://itai811.github.io/isoniche/reference/group_sites.md)
    if sites are small or incomplete.

3.  Estimate niche breadth with
    [`niche_breadth()`](https://itai811.github.io/isoniche/reference/niche_breadth.md)
    and check the diagnostics with
    [`summary()`](https://rdrr.io/r/base/summary.html).

4.  Plot niche breadth against abundance with
    [`iai_curve()`](https://itai811.github.io/isoniche/reference/iai_curve.md)
    and
    [`plot.iai_curve()`](https://itai811.github.io/isoniche/reference/plot.iai_curve.md).

Start with
[`vignette("isoniche")`](https://itai811.github.io/isoniche/articles/isoniche.md);
see also
[`vignette("data-preparation")`](https://itai811.github.io/isoniche/articles/data-preparation.md),
[`vignette("diagnostics")`](https://itai811.github.io/isoniche/articles/diagnostics.md)
and
[`vignette("method")`](https://itai811.github.io/isoniche/articles/method.md).

## Lower-level functions

[`fit_isodars()`](https://itai811.github.io/isoniche/reference/fit_isodars.md)
fits pairwise or n-dimensional isodars;
[`iai()`](https://itai811.github.io/isoniche/reference/iai.md) computes
niche breadth from fitted isodars at given total abundances;
[`isoderr()`](https://itai811.github.io/isoniche/reference/isoderr.md)
bootstraps it for one species;
[`pairwise_niche()`](https://itai811.github.io/isoniche/reference/pairwise_niche.md)
gives indices from the isodar parameters;
[`simulate_isodars()`](https://itai811.github.io/isoniche/reference/simulate_isodars.md)
simulates data from known isodars.

## References

Granot, I., Dubiner, S. & Belmaker, J. (2026). Why abundance-based
indices of niche breadth are biased, and what can be done to improve
them. *Ecology Letters* (in revision).

## See also

Useful links:

- <https://itai811.github.io/isoniche/>

- <https://github.com/itai811/isoniche>

- Report bugs at <https://github.com/itai811/isoniche/issues>

## Author

**Maintainer**: Itai Granot <itai811@gmail.com> (Contributed equally
with Shahar Dubiner)

Authors:

- Shahar Dubiner <dubiner@mail.tau.ac.il>
  ([ORCID](https://orcid.org/0000-0003-1771-0087)) (Contributed equally
  with Itai Granot)

- Jonathan Belmaker ([ORCID](https://orcid.org/0000-0002-5618-7359))
