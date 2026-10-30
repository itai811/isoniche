moorea_habs <- c("barrier_reef", "fringing_reef", "outer_slope")

# Three species from the paper's Table S2 (data-raw/moorea_fish.csv)
moorea_species <- function(species) {
  d <- utils::read.csv(test_path("fixtures", "moorea_ref.csv"))
  d[d$Species == species, ]
}

# Line through four habitats with ~45% of habitat-by-site cells missing and
# at least two habitats recorded at every site
sim_incomplete <- function(n = 300, p_miss = 0.57, noise = 0.5) {
  v <- c(1, 2, 3, 4) / sqrt(30)
  mu <- c(2, 5, 8, 11)
  t <- stats::runif(n, 0, 30)
  X <- outer(t, v) + matrix(mu, n, 4, byrow = TRUE) +
    stats::rnorm(n * 4, 0, noise)
  miss <- matrix(stats::runif(n * 4) < p_miss, n, 4)
  for (i in seq_len(n)) {
    while (sum(!miss[i, ]) < 2) miss[i, sample.int(4, 1)] <- FALSE
  }
  X[miss] <- NA
  colnames(X) <- paste0("h", 1:4)
  list(data = as.data.frame(X), direction = v)
}
