iso <- tibble::tibble(
  habitat_x = c("a", "a", "b"), habitat_y = c("b", "c", "c"),
  intercept = c(2, 0, 4), slope = c(2, 0.5, 1), sd = c(1, 2, 0.5)
)

test_that("indices match the Appendix S1 formulas", {
  out <- pairwise_niche(iso)
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("mean_intercept", "mean_slope", "slope_niche_breadth",
                      "intercept_niche_breadth", "n_isodars", "note"))
  expect_equal(out$mean_intercept, 2)
  expect_equal(out$mean_slope, 3.5 / 3)
  expect_equal(out$slope_niche_breadth, 1 / sqrt((1 + 0.25 + 0) / 3))
  expect_equal(out$intercept_niche_breadth, 1 / sqrt((4 + 0 + 16) / 3))
  expect_equal(out$n_isodars, 3)
})

test_that("weights are applied to means and deviations", {
  w <- c(1, 2, 3)
  out <- pairwise_niche(iso, weights = w)
  expect_equal(out$mean_intercept, sum(w * iso$intercept) / 6)
  expect_equal(out$mean_slope, sum(w * iso$slope) / 6)
  expect_equal(out$slope_niche_breadth,
               1 / sqrt(sum(w * (iso$slope - 1)^2) / 6))

  v <- 1 / (iso$sd^2 * (1 + iso$slope^2))
  expect_equal(pairwise_niche(iso, weights = "1/var"),
               pairwise_niche(iso, weights = v))
  # zero weight drops an isodar
  expect_equal(pairwise_niche(iso, weights = c(1, 1, 0))$mean_slope, 1.25)
})

test_that("log_slope treats reciprocal slopes symmetrically", {
  sym <- tibble::tibble(habitat_x = c("a", "a"), habitat_y = c("b", "c"),
                        intercept = c(0, 0), slope = c(2, 0.5))
  out <- pairwise_niche(sym, log_slope = TRUE)
  expect_equal(out$slope_niche_breadth, 1 / log(2))
  # on the linear scale 2 counts as a larger deviation than 0.5
  expect_equal(pairwise_niche(sym)$slope_niche_breadth,
               1 / sqrt((1 + 0.25) / 2))
})

test_that("isodars are oriented to non-negative intercepts", {
  neg <- iso
  neg$intercept[1] <- -2
  neg$slope[1] <- 2
  out <- pairwise_niche(neg)
  o <- attr(out, "isodars")
  expect_equal(c(o$habitat_x[1], o$habitat_y[1]), c("b", "a"))
  expect_equal(o$intercept[1], 1)
  expect_equal(o$slope[1], 0.5)
  expect_true(o$flipped[1])
  expect_true(all(o$intercept >= 0))

  # same answer whichever way fit_isodars() expressed the isodars
  d <- data.frame(a = 3:20, b = 2 * (3:20) - 5 + c(0.3, -0.3),
                  c = 3:20 + c(1, -1))
  expect_equal(
    pairwise_niche(fit_isodars(d, hab_cols = c("a", "b", "c"))),
    pairwise_niche(fit_isodars(d, hab_cols = c("a", "b", "c"),
                               flip_intercept = FALSE)),
    ignore_attr = TRUE
  )
})

test_that("a perfect generalist gives infinite indices", {
  gen <- tibble::tibble(habitat_x = "a", habitat_y = "b", intercept = 0,
                        slope = 1)
  out <- pairwise_niche(gen)
  expect_equal(out$slope_niche_breadth, Inf)
  expect_equal(out$intercept_niche_breadth, Inf)
})

test_that("non-positive slopes make the log-slope index NA", {
  pw <- fit_isodars(moorea_species("Scarus schlegeli"), hab_cols = moorea_habs)
  expect_true(any(pw$slope < 0))
  out <- pairwise_niche(pw, weights = "1/var", log_slope = TRUE)
  expect_true(is.na(out$slope_niche_breadth))
  expect_match(out$note, "non-positive slope")
  expect_true(is.finite(out$intercept_niche_breadth))
  expect_true(is.finite(pairwise_niche(pw)$slope_niche_breadth))
})

test_that("pairwise_niche() checks its arguments", {
  ln <- fit_isodars(moorea_species("Scarus schlegeli"), hab_cols = moorea_habs,
                    dim = "ndim")
  expect_error(pairwise_niche(ln), "pairwise isodars")
  expect_error(pairwise_niche(data.frame(x = 1)), "pairwise isodars")
  expect_error(pairwise_niche(iso[0, ]), "no rows")
  expect_error(pairwise_niche(iso, weights = c(1, 2)), "weights")
  expect_error(pairwise_niche(iso, weights = c(1, -1, 1)), "weights")
  expect_error(pairwise_niche(iso[, 1:4], weights = "1/var"), "'sd' column")
  bad <- iso
  bad$intercept[1] <- -1
  bad$slope[1] <- -1
  expect_error(pairwise_niche(bad), "non-negative intercept")
})
