################################################################################
# RD2D Package
# Data Simulation
################################################################################

rm(list = ls(all = TRUE))

# Build the SPP-calibrated cubic design matrix.
design_matrix <- function(dat) {
  x1 <- dat$x.1
  x2 <- dat$x.2
  cbind(
    1, x1, x2, x1^2, x1 * x2, x2^2,
    x1^3, x1^2 * x2, x1 * x2^2, x2^3
  )
}

# Store the SPP-calibrated cubic DGP.
dgp <- list(
  beta_y_0 = c(
    0.369916579111109, 0.00430768720228995, -0.00245733885625568,
    1.62105590793036e-05, 7.94581926007163e-06, 4.24074450908172e-05,
    2.29593705661776e-08, 1.59000624539961e-07, 3.45504841426239e-07,
    2.81256828567388e-07
  ),
  beta_y_1 = c(
    0.736166509744787, 0.000756347351213138, -0.00154115603887117,
    3.49990029700921e-05, 8.61468650013817e-05, -0.000155449166992341,
    -2.82846014355806e-07, 1.41889707426739e-07, -1.6593205900769e-06,
    3.94207017318509e-06
  ),
  beta_fuzzy_0 = c(
    -26.5660685226902, 3.23372932760228e-14, 1.40086045914661e-13,
    -6.50288763217046e-16, 7.6156958476755e-16, -1.30170689826067e-14,
    -2.72635230089794e-18, -3.95577524930275e-18, -8.48892391384158e-17,
    8.50472091483012e-17
  ),
  beta_fuzzy_1 = c(
    0.328585902510212, 0.00259026946365757, -0.00265595841237584,
    0.000215463378801299, -6.62666277106809e-06, -0.000563004965776261,
    -1.56069812328922e-06, 1.21170156753277e-07, 2.88676468236169e-06,
    1.28517906890237e-05
  ),
  lambda_0 = 0.0442625515378338,
  lambda_1 = 0.581534010672373,
  sigma_y_0 = 0.330524283558143,
  sigma_y_1 = 0.329017540357162
)

# Simulate one SPP-calibrated sample.
simulate_spp_cubic <- function(n, seed) {
  set.seed(seed)
  X <- data.frame(x.1 = 100 * rbeta(n, 3, 4) - 25,
                  x.2 = 100 * rbeta(n, 3, 4) - 25)
  assignment <- as.numeric(X$x.1 >= 0 & X$x.2 >= 0)
  R <- design_matrix(X)
  mu_y_0 <- as.numeric(R %*% dgp$beta_y_0)
  mu_y_1 <- as.numeric(R %*% dgp$beta_y_1)
  mu_fuzzy_0 <- plogis(as.numeric(R %*% dgp$beta_fuzzy_0))
  mu_fuzzy_1 <- plogis(as.numeric(R %*% dgp$beta_fuzzy_1))
  fuzzy0 <- as.numeric(runif(n) <= mu_fuzzy_0)
  fuzzy1 <- as.numeric(runif(n) <= mu_fuzzy_1)
  y0 <- mu_y_0 + dgp$lambda_0 * (fuzzy0 - mu_fuzzy_0) +
    rnorm(n, sd = dgp$sigma_y_0)
  y1 <- mu_y_1 + dgp$lambda_1 * (fuzzy1 - mu_fuzzy_1) +
    rnorm(n, sd = dgp$sigma_y_1)
  data.frame(x.1 = X$x.1, x.2 = X$x.2, assignment = assignment,
             fuzzy = ifelse(assignment == 1, fuzzy1, fuzzy0),
             Y = ifelse(assignment == 1, y1, y0))
}

# Set simulation inputs.
n <- as.integer(Sys.getenv("RD2D_ILLUSTRATION_N", "6000"))
seed <- 20260508

# Save the simulated data.
dat <- simulate_spp_cubic(n = n, seed = seed)
write.csv(dat, "rd2d_data.csv", row.names = FALSE)
