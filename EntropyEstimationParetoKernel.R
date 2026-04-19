# -------------------------------------------------------
# "Differential Entropy Estimation with a Paretian Kernel:
#  Tail Heaviness and Smoothing"
# 
# R functions for estimating differential entropy
# of heavy-tailed unidimensional data
#
# -------------------------------------------------------
# last update: 24/Apr/2024
# 
# Raul Matsushita
# Sergio Da Silva
# Helena S. Brandao
# Iuri Ribeiro Nobre
# ------------------------------------------------------------

# ========================================================
# Compute a leave one point out distance
# last update: 30/nov/2023
# --------------------------------------------------------

dlopout <- function(index, x){
  min(abs(x[index] - x[-index]))
}
# ========================================================

# ========================================================
# Gets differential entropy estimates based on 
# nearest neighbor distances
# last update: 30/nov/2023
# --------------------------------------------------------

H.NND <- function(x) {
  
  # data validation ---------------
  n.size  <- length(x)
  tie.chk <- ifelse(length(unique(x))<n.size,1,0)
  
  # tie.chk = 0 == no ties
  # tie.chk = 1 == ties
  if (tie.chk == 1) warning("ties are not allowed in this function", call. = TRUE)
  
  # -------------------------------
  index  <- 1:n.size
  rho    <- as.numeric(lapply(index,dlopout, x = x))
  H      <- mean(log(n.size*rho)) + log(2) - digamma(1)
  
  return(H)
}

# ========================================================
# Gets differential entropy estimates based on 
# sample-spacings of order m
# last update: 30/nov/2023
# --------------------------------------------------------
H.SSM <- function(x,m = 1){
  
  # data validation ---------------
  n.size  <- length(x)
  tie.chk <- ifelse(length(unique(x))<n.size,1,0)
  # tie.chk = 0 == no ties
  # tie.chk = 1 == ties
  if (tie.chk == 1) warning("ties are not allowed in this function", call. = TRUE)
  
  # -------------------------------
  x.s    <- sort(x)
  delta  <- x.s[(1+m):n.size] - x.s[1:(n.size-m)]
  f.hist <- m/(delta*n.size)
  H      <- -mean(log(f.hist)) + log(m) - digamma(m)
  
  return(H)
}

# ========================================================

# ========================================================
# Compute a leave-one-point-out Univariate Pareto's kernel 
# density estimate
# last update: 02/Feb/2024
# note: beta = smoothing parameter h
# --------------------------------------------------------

Klopout.P <- function(index, alpha, beta, x.var) {
  
  u          <- x.var[index] - x.var[-index]
  abs.u      <- abs(u)
  K.u        <- (0.5*alpha*(beta^alpha))/(beta + abs.u[abs.u>0])^(alpha+1)
  f.hat      <- mean(K.u,na.rm = TRUE)
  
  return(f.hat)
}
# ========================================================

Klopout.P <- Vectorize(Klopout.P, "index")

# ========================================================
# Gets univariate differential entropy estimates based on 
# the Pareto's kernel function
# last update: 08/Feb/2024
# note: beta = smoothing parameter h
# --------------------------------------------------------

Entropy.P <- function(x.var, par){
  {
    index         <- 1:length(x.var)
    f.hat         <- Klopout.P(index, alpha = par[1], beta = par[2], x.var = x.var)
    log.f.hat     <- log(f.hat) 
    H.xx          <- -mean(log.f.hat[log.f.hat>-Inf],na.rm = TRUE)
  }
  return(H.xx)
}
# ========================================================

# ========================================================
# Compute a leave one point out Univariate Gaussian kernel 
# density estimate
# last update: 08/Feb/2024
# --------------------------------------------------------

Klopout.G <- function(index, h, x.var) {
  
  u          <- (x.var[index] - x.var[-index])/h
  K.u        <- dnorm(u)
  f.hat      <- mean(K.u,na.rm = TRUE)/h
  
  return(f.hat)
}

# ========================================================

Klopout.G <- Vectorize(Klopout.G, "index")

# ========================================================
# Gets univariate differential entropy estimates based on 
# the Gaussian kernel function
# last update: 08/Feb/2024
# --------------------------------------------------------

Entropy.G <- function(x.var, h){
  {
    index         <- 1:length(x.var)
    f.hat         <- Klopout.G(index, h = h, x.var = x.var)
    log.f.hat     <- log(f.hat) 
    H.xx          <- -mean(log.f.hat[log.f.hat>-Inf],na.rm = TRUE)
  }
  
  return(H.xx)
}
# ========================================================

# ========================================================
# Compute a leave one point out Univariate Cauchy's kernel 
# density estimate
# last update: 08/Feb/2024
# --------------------------------------------------------

Klopout.C <- function(index, h, x.var){
  
  u          <- (x.var[index] - x.var[-index])/h
  K.u        <- dcauchy(u)
  f.hat      <- mean(K.u,na.rm = TRUE)/h
  
  return(f.hat)
}
# ========================================================

Klopout.C <- Vectorize(Klopout.C, "index")

# ========================================================
# Gets univariate differential entropy estimates based on 
# the Cauchy's kernel function
# last update: 08/Feb/2024
# --------------------------------------------------------

Entropy.C <- function(x.var, h){
  {
    index         <- 1:length(x.var)
    f.hat         <- Klopout.C(index, h = h, x.var = x.var)
    log.f.hat     <- log(f.hat) 
    H.xx          <- -mean(log.f.hat[log.f.hat>-Inf],na.rm = TRUE)
  }
  
  return(H.xx)
}
# ========================================================

# ========================================================
# Compute a leave-one-point-out Univariate Laplacian kernel 
# density estimate
# last update: 08/Feb/2024
# --------------------------------------------------------

Klopout.L <- function(index, h, x.var){
  
  u          <- (x.var[index] - x.var[-index])/h
  K.u        <- 0.5*exp(-abs(u))
  f.hat      <- mean(K.u,na.rm = TRUE)/h
  
  return(f.hat)
}
# ========================================================

Klopout.L <- Vectorize(Klopout.L, "index")

# ========================================================
# Gets univariate differential entropy estimates based on 
# the Laplacian kernel function
# last update: 08/Feb/2024
# --------------------------------------------------------

Entropy.L <- function(x.var, h){
  {
    index         <- 1:length(x.var)
    f.hat         <- Klopout.L(index, h = h, x.var = x.var)
    log.f.hat     <- log(f.hat) 
    H.xx          <- -mean(log.f.hat[log.f.hat>-Inf],na.rm = TRUE)
  }
  
  return(H.xx)
}
# ========================================================
