rm(list=ls(all=TRUE))

############################## Model parameters ################################

a0 <- 3.350e-01
b0 <- 2.518e-03
c0 <- -1.272e-03
sd0 <- 0.3321

a1 <- 0.6979117
b1 <- 0.0027384
c1 <- -0.0006046
sd1 <- 0.4352

############################## Generate Data ###################################

set.seed(127)

alpha_1 <- 3
beta_1 <- 4
alpha_2 <- 3
beta_2 <- 4

n <- 20 * 20 * 25 * 2

x.1 <- 100 * rbeta(n, alpha_1, beta_1) - 25
x.2 <- 100 * rbeta(n, alpha_2, beta_2) - 25

X <- cbind(x.1, x.2)
X <- as.data.frame(X)
t <- x.1 >= 0 & x.2 >= 0
y <- (1 - t) * (a0 + b0 * x.1 + c0 * x.2) + t * (a1 + b1 * x.1 + c1 * x.2)
y <- y * 2 + rnorm(n,mean = 0, sd = sd0) * (1 - t) + rnorm(n,mean = 0, sd = sd1) * t

data_rd2d <- data.frame(x.1 = x.1, x.2 = x.2, t = t, y = y)

# Evaluation points
neval <- 40
eval <- matrix(nrow = neval, ncol = 2)
for (i in 1: ceiling(neval * 0.5)){
  eval[i,] <- c(0, 50 - (i-1) * 50 / ceiling(neval * 0.5))
}
for (i in (ceiling(neval * 0.5)+1): neval){
  eval[i,] <- c((i - ceiling(neval * 0.5) - 1) *50 / (ceiling(neval * 0.5)),0)
}
eval <- data.frame(eval)
colnames(eval) <- c("x.1", "x.2")

# Distances to evaluation points
D <- proxy::dist(X, eval, method = "euclidean")  # Use "euclidean" for Euclidean distances
t_expanded <- matrix(rep(2 * t - 1, times = ncol(D)), nrow = nrow(D), ncol = ncol(D))
D <- D * t_expanded

############################## Save Data #######################################

write.csv(data_rd2d, "Data/data_rd2d.csv", row.names = FALSE)
write.csv(eval, "Data/eval.csv", row.names = FALSE)
write.table(D, file = "Data/D.csv", sep = ",", row.names = FALSE, col.names = FALSE)


