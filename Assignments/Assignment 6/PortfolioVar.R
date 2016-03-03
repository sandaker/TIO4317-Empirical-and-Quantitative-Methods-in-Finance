library("Rsolnp")
library("lpSolve")
library("fitdistrplus")

# Set correct path in R
setwd("/Users/Kristian/Kristian/Dropbox/Studier/8. semester/TIO4317 Empirisk finans/Assignments/Assignment 6")

# Load the data set of stock and bond indices
df <- read.table("StocksBonds.csv", header=TRUE, sep=",")

### PARAMETRIC INPUT
assetNames <- c("MSDUUS.Index", "MSEUSAG.Index") # Assets to be compared
alpha <- 0.95 # Alpha level VaR
n_days <- 100 # Number of days to sample
L <- 100 # Number of scenarios
n_optimizations <- 6 # Number of optimizations and plots

### Model

P <- cbind(df[assetNames]) # Absolute price vector
r <- log(P[2:L,]) - log(P[1:(L-1),]) # Log returns
r[mapply(is.infinite, r)] <- 0 # Fix problematic values

# Do several optimizations
pdf(file = paste(getwd(), "/plots/var_plot_", n_optimizations, ".pdf", sep=""))
par(mfcol=c(n_optimizations/2,2)) # Plotting grid
par(mar=c(4,5,4,5))
par(cex=0.5)

for (k in 1:n_optimizations) {
  
  # Sample n-days log returns, L times, scale up
  R <- matrix(0, L, 2)
  for (i in 1:L) {
    r_1 <- sample(r[,1], n_days, replace = TRUE)
    r_2 <- sample(r[,2], n_days, replace = TRUE)
    R[i,] <- cbind(exp(sum(r_1)), exp(sum(r_2)))
  }
  
  
  # Lambda linspace, [0,1]
  Lambda <- seq(0,1, length = 100)
  
  # Set up the initial decision variable vector
  objective.in <- t(t(c(1, rep(0, L))))
  binary.vec <- t(seq(2,(L+1)))
  
  # Initialize V_0
  V_0 <- 1
  
  # Calculating the correct M
  M <- V_0 + V_0*sum(apply(R, 2, min))
  
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
          - lambda*V_0*R[l,1] - (1-lambda)*V_0*R[l,2]))
    }
    zeta <- lp(direction = "min", objective.in = objective.in, const.mat = const.mat,
       const.dir = const.dir, const.rhs = const.rhs, binary.vec = binary.vec)$solution[1]
    Zeta[i] <- zeta
  }
  plot(c(Zeta) ~ c(Lambda), type = "l", 
       main = bquote("VaR at "*.(alpha*100)*"% "*.(assetNames[1])*" vs. "*.(assetNames[2])),
       xlab = expression(lambda),
       ylab = expression(zeta))
  par(new=F)
}

dev.off()








