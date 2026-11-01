# Builds data/moorea_fish.rda from the paper's Table S2 (moorea_fish.csv):
# one row per species x site where the species was recorded, one column per
# habitat. Reshaped to long format with explicit zeros.
wide <- utils::read.csv("data-raw/moorea_fish.csv", stringsAsFactors = FALSE)
habs <- c("barrier_reef", "fringing_reef", "outer_slope")
moorea_fish <- data.frame(
  species = rep(wide$Species, each = length(habs)),
  site = rep(wide$Site, each = length(habs)),
  habitat = rep(habs, times = nrow(wide)),
  abundance = as.integer(t(as.matrix(wide[habs]))),
  stringsAsFactors = FALSE
)
moorea_fish <- moorea_fish[order(moorea_fish$species, moorea_fish$site), ]
rownames(moorea_fish) <- NULL
moorea_fish <- tibble::as_tibble(moorea_fish)
save(moorea_fish, file = "data/moorea_fish.rda", compress = "xz")
