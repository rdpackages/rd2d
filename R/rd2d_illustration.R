# rd2d: illustration file
# Authors: M. D. Cattaneo, R. Titiunik, R. R. Yu
# Last update: 2025/5/14

rm(list=ls(all=TRUE))

if (!requireNamespace("rd2d", quietly = TRUE)) install.packages("rd2d")
library(rd2d)

data_dir <- if (dir.exists("Data")) "Data" else file.path("R", "Data")

################################ Load data #####################################

data_rd2d <- read.csv(file.path(data_dir, "data_rd2d.csv"))
eval <- read.csv(file.path(data_dir, "eval.csv"))
D <- as.matrix(read.table(file.path(data_dir, "D.csv"), sep = ",", header = FALSE))

y <- data_rd2d$y
X <- data_rd2d[,c("x.1", "x.2")]
t <- data_rd2d$t

############################### Bivariate method ###############################

result.rd2d <- rd2d(y, X, t, eval, stdvars = FALSE)

print(result.rd2d)
summary(result.rd2d)

# options for summary function
summary(result.rd2d, subset = c(1,5,10,15,21,25,30,35,40)) # show selected evaluation points
summary(result.rd2d, subset = c(1,5,10,15,21,25,30,35,40), CBuniform = TRUE) # show confidence bands
neval <- nrow(eval)
summary(result.rd2d, subset = c(1,5,10,15,21,25,30,35,40), WBATE = rep(1,neval)) # show weighted boundary average treatment effect
summary(result.rd2d, subset = c(1,5,10,15,21,25,30,35,40), output = "bw") # show bandwidth information

# bandwidth selection
bws.rd2d <- rdbw2d(y, X, t, eval)
summary(bws.rd2d,subset = c(1,5,10,15,21,25,30,35,40))

############################## Distance method #################################

result.dist <- rd2d.dist(y,D, b = eval)

print(result.dist)
summary(result.dist, subset = c(1,5,10,15,21,25,30,35,40))
summary(result.dist, subset = c(1,5,10,15,21,25,30,35,40), CBuniform = TRUE)
summary(result.dist, subset = c(1,5,10,15,21,25,30,35,40), WBATE = rep(1,neval))
summary(result.dist, subset = c(1,5,10,15,21,25,30,35,40), output = "bw")

result.dist.kinkon <- rd2d.dist(y,D, b = eval, kink = "on") # kink adjustment
summary(result.dist.kinkon, subset = c(1,5,10,15,21,25,30,35,40))

# bandwidth selection

bws.dist <- rdbw2d.dist(y,D)
summary(bws.dist,subset = c(1,5,10,15,21,25,30,35,40))

bws.dist.kinkon <- rdbw2d.dist(y,D, kink = "on")
summary(bws.dist.kinkon,subset = c(1,5,10,15,21,25,30,35,40))
