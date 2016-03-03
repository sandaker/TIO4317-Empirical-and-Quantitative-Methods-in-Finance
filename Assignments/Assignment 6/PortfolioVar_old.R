library("Rsolnp")
library("lpSolve")
library("fitdistrplus")

# Set correct path in R
setwd("/Users/Kristian/Kristian/Dropbox/Studier/8. semester/TIO4317 Empirisk finans/Assignments/Assignment 6")

# Load the data set of stock and bond indices
df <- read.table("StocksBonds.csv", header=TRUE, sep=",")

# Assets to be compared
assetNames <- c("MSDUUS.Index", "MSEUSAG.Index")
P <- cbind(df[assetNames]) # Absolute price vector
L <- dim(P)[1] # Number of scenarios
r <- log(P[2:L,]) - log(P[1:(L-1),]) # Log returns
r[mapply(is.infinite, r)] <- 0 # Fix problematic values

# Sample 100-day log returns, 100 times, scale up
n_days <- 100
n_times <- 100
R <- matrix(0, n_times, 2)
for (i in 1:n_times) {
  r_1 <- sample(r[,1], n_days, replace = TRUE)
  r_2 <- sample(r[,2], n_days, replace = TRUE)
  R[i,] <- cbind(exp(sum(r_1)), exp(sum(r_2)))
}


### Variables and parameters
n <- dim(df)[2]-1     # Number of unique stock/bond indices
df <- df[1:100,]
L <- dim(df)[1]       # Number of scenarios

# Vector of two index assets
indexNames <- colnames(df)
A <- cbind(df[,1], df[,2], df[,8])
colnames(A) <- c(indexNames[1], indexNames[2], indexNames[8])
P <- A


# Lambda linspace, [0,1]
Lambda <- seq(0,1, length = 100)

# Alpha level VaR
alpha <- 0.95

# Set up the initial decision variable vector
objective.in <- t(t(c(1, rep(0, L))))
binary.vec <- t(seq(2,(L+1)))

# Initialize V_0
V_0 <- P[1,2]*0.5 + P[1,3]*0.5

# Calculating the correct M
M <- V_0 + sum(apply(P, 2, max)[2:3])

# Set up the coefficient matrix
const.mat <- t(c(0, rep(1/L, L)))
for (l in 1:L) {
  constrVec <- t(c(1, rep(0, L)))
  constrVec[l+1] <- -M
  const.mat <- rbind(const.mat, constrVec)
}

# Setting up the equality vector
const.dir <- t(rep(">=", (L+1)))

# Optimization loop for every lambda

Zeta <- t(rep(0, length(Lambda)))

for (i in 1:length(Lambda)) {
  lambda <- Lambda[i]
  # Set up the constant RHS
  const.rhs <- t(alpha)
  for (l in 1:L) {
    const.rhs <- rbind(const.rhs, t(V_0 - M 
        - lambda*P[l,2] - (1-lambda)*P[l,3]))
  }
  zeta <- lp(direction = "min", objective.in = objective.in, const.mat = const.mat,
     const.dir = const.dir, const.rhs = const.rhs, binary.vec = binary.vec)$solution[1]
  Zeta[i] <- zeta
}

plot(c(Zeta) ~ c(Lambda), type = "l")


### Functions 

  # V(x; \tilde(P)) function
  VP <- function(x, P) {
    return (t(P) %*% cbind(x[1:n]))
  }

  # V(x; \bar{P}) function
  # Also constr1 in model
  VPmean <- function(x) {
    return (t(Pmean) %*% t(t(x[1:n])))
  }

  # Objective function
  objective <- function(x) {
    return (t(pl) %*% t(t(x[(n+1):(L+n)])))
  }

### Equality constraints

  eqfun <- function(x) {
    return ((as.numeric(df[1, 2:(n+1)]) %*% cbind(x[1:n])) - V0)
  }

# RHS of constraints
eqB <- 0
ineqLB <- cbind(rep(0,1+2*L))
ineqUB <- rep(10000375, 1+2*L)

LB <- rep(0, n+L)
UB <- rep(100000000, n+L)

mus <- seq(0, 1, length=20)

# Mean absolute deviation optimizations loop
for (i in 1:length(mus)) {
  cat("Optimizing...")
  mu <- mus[i]
  
  # All inequalities in the model
  ineqfun <- function(x) {
    z <- rep(0, 1 + 2*L)
    z[1] <- VPmean(x) - mu*V0
    for (i in 1:L) {
      Pl <- cbind(as.numeric(df[i,2:(n+1)]))
      z[1+i] <- VPmean(x) - VP(x, Pl) + x[n+i]
      z[1+L+i] <- VP(x, Pl) - VPmean(x) + x[n+i]
    }
    return (cbind(z))
  }
  
  cat("Starting.")
  opt <- gosolnp(pars = x0, fun = objective, eqfun = eqfun, eqB = eqB,
               ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB, LB = LB, UB = UB)
  pars <- opt$pars
  if (i == 1) {
    parMatrix <- rbind(pars)
  } else {
    parMatrix <- rbind(parMatrix, pars)
  }
}