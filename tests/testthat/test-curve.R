test_that("iai_curve() starts at the DI value and approaches the asymptote", {
  nb <- do.call(niche_breadth, c(list(moorea_long(ref_species)), long_cols,
                                 list(n_boot = 0, dim = "ndim")))
  cu <- iai_curve(nb)
  expect_s3_class(cu, "iai_curve")
  expect_named(cu, c("species", "total_abundance", "relative_abundance",
                     "niche_breadth", "lower", "upper"))
  pts <- attr(cu, "points")
  expect_equal(pts$spillover_abundance, nb$spillover_abundance)
  expect_equal(pts$nb_dd, nb$nb_dd)
  for (i in 1:3) {
    d <- cu[cu$species == ref_species[i], ]
    expect_equal(d$relative_abundance[1], 1)
    expect_equal(d$niche_breadth[1], nb$nb_di[i])
    fit <- attr(nb, "fits")$fit[[i]]
    expect_equal(d$niche_breadth[50], iai(fit, abundances = d$total_abundance[50])$niche_breadth)
  }
  expect_true(all(is.na(cu$lower)))
})

test_that("crossings between species curves are found", {
  nb <- do.call(niche_breadth, c(list(moorea_long(ref_species)), long_cols,
                                 list(n_boot = 0, dim = "ndim")))
  cu <- iai_curve(nb, totals = 10^seq(1, 4, length.out = 300))
  cr <- attr(cu, "crossings")
  expect_gt(nrow(cr), 0)
  fits <- attr(nb, "fits")
  for (k in seq_len(nrow(cr))) {
    a <- iai(fits$fit[[match(cr$species_1[k], ref_species)]],
             abundances = cr$total_abundance[k])$niche_breadth
    b <- iai(fits$fit[[match(cr$species_2[k], ref_species)]],
             abundances = cr$total_abundance[k])$niche_breadth
    expect_equal(a, b, tolerance = 0.01)
  }
})

test_that("bootstrap bands and plots work", {
  fits <- do.call(fit_isodars, c(list(moorea_long(ref_species[2:3])), long_cols,
                                 list(dim = "ndim")))
  cu <- iai_curve(fits, n_boot = 30, seed = 1)
  expect_true(all(cu$lower <= cu$upper + 1e-12, na.rm = TRUE))
  expect_gt(mean(is.finite(cu$lower)), 0.9)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(cu))
  expect_invisible(plot(cu, x_axis = "relative", species = ref_species[3]))
  expect_invisible(plot(cu, abundances = c(100, 1000), crossings = TRUE))
  expect_error(plot(cu, species = "Nemo"), "not in")
  # a single wide fit, no bootstrap
  single <- iai_curve(fits$fit[[1]])
  expect_false("species" %in% names(single))
  expect_invisible(plot(single))
  expect_error(iai_curve(fits$fit[[1]], n_boot = 5), "need the data")
  expect_error(iai_curve(fits, totals = -1), "totals")
})
