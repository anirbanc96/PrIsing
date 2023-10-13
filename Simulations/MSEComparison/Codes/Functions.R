################################################################################
################################ Packages ######################################

library(foreach)
library(doParallel)
library(IsingSampler)
library(igraph)
library(tidyverse)
library(latex2exp)
library(VGAM)

################################################################################

############# Evaluating Pseudo-Likelihood (Usual and Private) #################

L.sigma.beta = function(beta, N, J.N, sigma,
                        priv = F, priv.ty = NULL,
                        Delta = NULL, gamma = NULL, seed = NULL){
  
  ##############################################################################
  # INPUTS: beta    <- value of beta parameter
  #         N       <- size of the underlying graph
  #         J.N.    <- adjacency matrix in the Hamiltonian
  #         sigma   <- size n vector coming from Ising Model
  #         priv    <- indicator for private estimation
  #         priv.ty <- 1 or 0 indicating (epsilon,delta) or (epsilon,0) privacy
  #         Delta   <- Capital Delta used in Algorithm 1
  #         gamma   <- standard deviation of the additive noise
  #         seed    <- seed to generate gaussian noise
  ##############################################################################
  # OUTPUT: L.val      <- value of non-private pseudo-likelihood
  #         L.priv.val <- value of private pseudo-likelihood with addition of 
  #                       Delta and noise
  ##############################################################################
  
  sigma <- t(sigma)
  m.sigma.vec <- Rfast::mat.mult(J.N, sigma)
  L.val <- sum(m.sigma.vec * (sigma - tanh(beta * m.sigma.vec)))/N
  
  if (priv == F){
    
    return (L.val)
    
  }
  else{
    
    set.seed(seed)
    if (priv.ty == 1){
      b <- rnorm(1, mean = 0, sd = gamma)
    }
    else{
      b <- VGAM::rlaplace(1, location = 0, scale = gamma)
    }
    L.priv.val <- -L.val + ((Delta * beta) + b)/N
    return(L.priv.val)
    
  }
}

################################################################################

###################### MPLE (private and non-private) ##########################
beta.hat.N <- function(sigma, N, J.N,
                       priv = F, priv.ty = NULL,
                       Delta = NULL, gamma = NULL, seed = NULL,
                       range.up = 30, range.low = -200){
  
  ##############################################################################
  # INPUTS: sigma     <- the spin vector sampled from Ising Model
  #         N         <- size of the underlying graph
  #         J.N       <- matrix in the Hamiltonian
  #         priv      <- indicator for private estimation
  #         priv.ty   <- 1 or 0 indicating (epsilon,delta) or (epsilon,0)
  #                      privacy
  #         Delta     <- Capital Delta used in Algorithm 1
  #         gamma     <- standard deviation of the additive noise
  #         seed      <- seed to generate gaussian noise
  #         range.up  <- the upper range of uniroot
  #         range.low <- the lower range of uniroot
  ##############################################################################
  # OUTPUT: beta.solve <- the solution of beta found by solving (private) log-
  #                       likelihood equation
  ##############################################################################
  pseudo.likelihood.solve <- uniroot(L.sigma.beta,
                                     interval = c(range.low,range.up),
                                     N = N, J.N = J.N, sigma = sigma, 
                                     priv = priv, priv.ty = priv.ty,
                                     Delta = Delta, gamma = gamma, seed = seed)
  
  beta.solve <- pseudo.likelihood.solve$root
  return (beta.solve)
  
}

################################################################################

######################## Degree vector in Algo 1################################
d.func <- function(N, J.N){
  
  ##############################################################################
  # INPUTS: N   <- size of underlying graph
  #         J.N <- matrix in the Hamiltonian
  ##############################################################################
  # OUTPUT: vector of degrees - (d_{1},..,d_{N}) from Algorithm 1
  ##############################################################################
  
  return (N * rowSums(abs(J.N)))
  
}

################################################################################

######################## Zeta parameter from Algo 1 ############################
zeta.func <- function(N, J.N){
  
  ##############################################################################
  # INPUTS: N   <- size of underlying graph
  #         J.N <- matrix in the Hamiltonian
  ##############################################################################
  # OUTPUT: The parameter $\zeta$ from Algorithm 1
  ##############################################################################
  return ((8/N) * max(d.func(N, J.N)))
  
}

################################################################################

####################### The Delta parameter in Algo 1 ##########################
Delta.func <- function(epsilon, N, J.N){
  
  ##############################################################################
  # INPUTS: epsilon <- epsilon parameter of privacy
  #         N       <- size of underlying graph
  #         J.N     <- matrix in the Hamiltonian
  ##############################################################################
  # OUTPUT: Delta <- The parameter $\Delta$ from Algorithm 1
  ##############################################################################
  
  d.vec <- as.matrix(d.func(N, J.N))
  Delta <- (24/(epsilon * N)) * max(Rfast::mat.mult(abs(J.N), d.vec))
  return (Delta)
  
}

################################################################################

######################## Noise parameter gamma #################################
gamma.func <- function(N, J.N, epsilon, delta){
  
  ##############################################################################
  # INPUTS: N       <- size of underlying graph
  #         J.N     <- matrix in the Hamiltonian
  #         epsilon <- epsilon parameter for privacy
  #         delta   <- delta parameter for privacy
  ##############################################################################
  # OUTPUT: gamma <- noise parameter from Algo 1
  ##############################################################################
  
  if (delta == 0){
    
    zeta <- zeta.func(N, J.N)
    gamma <- 2*(zeta/epsilon)
    
    return(gamma)
  }
  else{
    
    zeta <- zeta.func(N, J.N)
    gamma <- (zeta * sqrt(8 * log(2/delta) + 4 * epsilon))/epsilon
    
    return (gamma)
  }
  
}

################################################################################

################## Estimating beta for multiple experiments#####################

sim.beta.est <- function(n.rep, beta.true, N, J.N,
                         priv = F, epsilon = NULL, delta = NULL){
  
  ##############################################################################
  # INPUTS: n.rep     <- number of expreiments done
  #         beta.true <- the value of true beta
  #         N         <- size of underlying graph
  #         J.N       <- matrix in the Hamiltonian
  #         priv      <- indicator for private estimation
  #         epsilon   <- epsilon parameter for privacy
  #         delta     <- delta parameter for privacy
  ##############################################################################
  # OUPUT: out.beta.hat <- If priv = F -> a n.rep length vector of non-private
  #                                       beta estimates
  #                        If priv = T -> a n.rep length vector of private
  #                                       beta estimates
  ##############################################################################
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  
  registerDoParallel(cl)
  
  writeLines(c("Running Function To Estimate beta"), "log.txt")
  
  if (priv == T){
    
    Delta <- Delta.func(epsilon, N, J.N)
    gamma <- gamma.func(N, J.N, epsilon, delta)
    
  }
  
  out.beta.hat <- foreach(k=1:n.rep, .combine = rbind, .export = ls(envir=globalenv())) %dopar% {
    
    cat(paste("Starting iteration",k,"\n"), file = "log.txt", append = T)
    
    sigma = IsingSampler::IsingSampler(n = 1, J.N, thresholds = rep(0,N),
                                       beta = beta.true, 
                                       responses = c(-1L,1L), method = "MH") 
    
    if (priv == F) {
    
      beta.hat <- beta.hat.N(sigma, N, J.N)
    
    }
    
    else {
      
      seed <- 2*k
      
      priv.ty = (delta>0)
      
      beta.hat.priv <- beta.hat.N(sigma, N, J.N,
                                  priv = T, priv.ty = priv.ty,
                                  Delta = Delta, gamma = gamma,
                                  seed = seed)
      
    }
    
    cat(paste("Completed iteration",k ,"\n"), file = "log.txt", append = T)
    
    if (priv == T){
      
      return (beta.hat.priv)
      
    }
    else{
      
      return (beta.hat)
      
    }
    
  }
  
  stopCluster(cl)
  return (out.beta.hat)
  
}

