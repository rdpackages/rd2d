################################################################################
# RD2D Package
# Numerical Illustration
################################################################################
rm(list = ls(all = TRUE))

library(rd2d)

# Generate boundary evaluation points.
make_eval_grid <- function(neval = 40) {
  half <- ceiling(neval / 2)
  rbind(
    data.frame(x.1 = rep(0, half),
               x.2 = 40 - (seq_len(half) - 1) * 40 / half),
    data.frame(x.1 = (seq_len(neval - half) - 1) * 56 / half,
               x.2 = rep(0, neval - half))
  )
}

# Generate signed distances to each boundary evaluation point.
make_signed_distances <- function(X, eval, assignment) {
  distance <- sapply(seq_len(nrow(eval)), function(j) {
    sqrt((X$x.1 - eval$x.1[j])^2 + (X$x.2 - eval$x.2[j])^2)
  })
  distance * matrix(2 * assignment - 1, nrow = nrow(distance),
                    ncol = ncol(distance))
}

# Load simulated data.
dat <- read.csv("rd2d_data.csv")
X <- dat[, c("x.1", "x.2")]
Y <- dat$Y
A <- dat$assignment
D <- dat$fuzzy

# Set illustration inputs.
neval <- as.integer(Sys.getenv("RD2D_ILLUSTRATION_NEVAL", "40"))
repp <- as.integer(Sys.getenv("RD2D_ILLUSTRATION_REPP", "499"))
eval <- make_eval_grid(neval)
distance <- make_signed_distances(X, eval, A)
wbate_weights <- rep(1, neval)
selected <- c(1, 5, 10, 15, 21, 25, 30, 35, 40)
selected <- selected[selected <= neval]

# Location-based bandwidth selection.
bw_location <- rdbw2d(Y, X, A, eval, masspoints = "off")
summary(bw_location, subset = selected)

# Location-based fuzzy estimation.
fit_location <- rd2d(
  Y, X, A, eval,
  fuzzy = D,
  params.other = "itt.0",
  params.cov = c("main", "itt", "fs", "itt.0"),
  masspoints = "off"
)

# Location-based fuzzy main effect.
summary(
  fit_location, output = "main", subset = selected, cbands = "main",
  WBATE = wbate_weights, LBATE = TRUE, repp = repp
)

# Location-based reduced-form effect.
summary(
  fit_location, output = "itt", subset = selected, cbands = "itt",
  WBATE = wbate_weights, LBATE = TRUE, repp = repp
)

# Location-based first-stage effect.
summary(
  fit_location, output = "fs", subset = selected, cbands = "fs",
  WBATE = wbate_weights, LBATE = TRUE, repp = repp
)

# Location-based control-side reduced-form effect.
summary(
  fit_location, output = "itt.0", subset = selected, cbands = "itt.0",
  WBATE = wbate_weights, LBATE = TRUE, repp = repp
)

# Distance-based bandwidth selection.
bw_distance <- rdbw2d.distance(Y, distance, b = eval, masspoints = "off")
summary(bw_distance, subset = selected)

# Distance-based sharp estimation.
fit_distance <- rd2d.distance(
  Y, distance = distance, b = eval, masspoints = "off", cbands = TRUE
)
summary(
  fit_distance, output = "main", subset = selected, cbands = "main",
  repp = repp
)

# Distance-based fuzzy bandwidth selection.
bw_distance_fuzzy <- rdbw2d.distance(
  Y, distance = distance, b = eval, fuzzy = D,
  bwparam = "itt", masspoints = "off"
)
summary(bw_distance_fuzzy, subset = selected)

# Distance-based fuzzy estimation.
fit_distance_fuzzy <- rd2d.distance(
  Y, distance = distance, b = eval,
  fuzzy = D, bwparam = "itt",
  params.cov = c("main", "itt", "fs"),
  masspoints = "off"
)

# Distance-based fuzzy main effect.
summary(
  fit_distance_fuzzy, output = "main", subset = selected, cbands = "main",
  WBATE = wbate_weights, LBATE = TRUE, repp = repp
)

# Distance-based reduced-form effect.
summary(
  fit_distance_fuzzy, output = "itt", subset = selected, cbands = "itt",
  WBATE = wbate_weights, LBATE = TRUE, repp = repp
)

# Distance-based first-stage effect.
summary(
  fit_distance_fuzzy, output = "fs", subset = selected, cbands = "fs",
  WBATE = wbate_weights, LBATE = TRUE, repp = repp
)
