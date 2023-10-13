################################################################################
############################# Required Packages ################################

require(tidyverse)

# Load required functions

source("Functions.R")

# Load data and Hamiltonian of Ising Model

source("DataManage.R")

################################################################################
############################## Additional Functions ############################

# Function for providing a vector of private MPLE estimates based on user given
# seeds.

priv.beta <- function(sigma, N, J.N, epsilon, delta, n.rep = 500,
                      seed.vec = 1:n.rep){
  
  ##############################################################################
  # INPUTS: sigma     <- size n vector coming from Ising Model
  #         N         <- size of the underlying graph
  #         J.N.      <- adjacency matrix in the Hamiltonian
  #         epsilon   <- epsilon parameter for privacy
  #         delta     <- delta parameter for privacy
  #         n.rep     <- number of repetitions performed
  #         seed.vec  <- n.rep sized vector of seeds needed to noise generation 
  ##############################################################################
  # OUTPUT: beta.est.priv <- n.rep sized vector of private MPLE estimates
  ##############################################################################
  
  # n.rep sized vector of seeds needed to noise generation 
  seed.vec <- (1:n.rep)
  
  # Value of Delta from Algorithm
  Delta <- Delta.func(epsilon, N, J.N)
  
  # Value of gamma (noise variance) in Algorithm
  gamma <- gamma.func(N, J.N, epsilon, delta)
  
  # Type of privacy -> (epsilon, delta) or (epsilon,0)
  priv.ty <- (delta>0)
  
  # parallelising for faster code
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  
  registerDoParallel(cl)
  
  writeLines(c("Running Function To Estimate private beta"), "log.txt")
  
  # running loop for n.rep many private MPLE estimate of beta
  
  beta.est.priv <- foreach(k=1:n.rep, .combine = rbind, .export = ls(envir=globalenv())) %dopar% {
    
    cat(paste("Starting iteration",k,"\n"), file = "log.txt", append = T)
    
    beta.est.priv.iter <- beta.hat.N(sigma, N, J.N,
                                     priv = T, priv.ty = priv.ty,
                                     Delta = Delta, gamma = gamma,
                                     seed = seed.vec[k])
    
    cat(paste("Completed iteration",k ,"\n"), file = "log.txt", append = T)
    
    beta.est.priv.iter
    
  }
  stopCluster(cl)
  
  return (as.vector(beta.est.priv))
}

################################################################################

# Function for producing n.errorbar many MSE values

MSE.func <- function(n.iter = 500, vec, mean, n.errorbar = 100){
  
  ##############################################################################
  # INPUTS: n.iter     <- number of values taken to get a single value of MSE
  #         vec        <- the vector of values (n.iter x n.errorbar length)
  #         mean       <- the value around which we calculate MSE
  #         n.errorbar <- number of MSE values needed to estimate sd of MSE's
  ##############################################################################
  # OUTPUT: MSE.vec <- n.errorbar sized vector of values of MSE around mean
  ##############################################################################
  
  MSE.vec <- rep(0, n.errorbar)
  for (i in 1:n.errorbar){
    
    vec.iter <- vec[(n.iter*(i-1) + 1):(n.iter * i)]
    MSE.vec[i] <- sum((vec.iter - mean)^2)/n.iter
    
  }
  
  return (MSE.vec)
  
}

################################################################################
# MPLE estimate of beta

beta.hat.nonpriv <- beta.hat.N((2*HIV.status.df$HIV-1), N, J.N)

# number of values of MSE needed to build errorbar
n.errorbar <- 100

# number of estimates needed to estimate MSE
n.iter <- 500

# total number of estimates needed to estimate and build errorbar of MSE
n.rep <- n.iter * n.errorbar


# Privacy Parameters
epsilon.vec <- c(5, 7.5, 10, 20, 50, 75, 100)
epsilon.len <- length(epsilon.vec)
delta <- 1/N

# Matrix to store values of MSE for each epsilon
MSE.mat <- matrix(0, epsilon.len, n.errorbar)

for (i in 1:epsilon.len){
  
  beta.est.priv <- priv.beta((2*HIV.status.df$HIV-1), N, J.N,
                                     epsilon.vec[i], delta, n.rep = n.rep)
  
  MSE.mat[i,] <- MSE.func(n.iter, beta.est.priv, beta.hat.nonpriv,
                          n.errorbar)
  
  cat(c("Finished computing cost for epsilon = ", epsilon.vec[i], "\n"))
  
}

# saving values of MSE for each epsilon
MSE.dat <- as_tibble(cbind(epsilon.vec, MSE.mat))

write.csv(MSE.dat, file = "MSEvsEpsilonBetahatPNP.csv")

################################################################################

# Run the following lines for re-plotting from old saved data.

# MSE.dat <- read.csv("MSEvsEpsilonBetahatPNP.csv")[,-1]
# epsilon.vec <- MSE.dat$epsilon.vec
# epsilon.len <- length(epsilon.vec)

################################################################################

# Creating matrix with values of epsilon, mean and mean +- sd of cost

MSE.plot.mat <- as_tibble(matrix(0, nrow = epsilon.len, ncol = 4))
colnames(MSE.plot.mat) <- c("epsilon", "mean", "up", "down")

MSE.plot.mat <- MSE.plot.mat %>%
  mutate(epsilon = epsilon.vec) %>%
  mutate(mean = apply(MSE.dat[,-1], 1, mean)) %>%
  mutate(up = apply(MSE.dat[,-1], 1, sd)) %>%
  mutate(down = apply(MSE.dat[,-1], 1, sd))

# Making epsilon values as factor for making plots

MSE.plot.mat[,1] <- as.factor(epsilon.vec)

# Plotting cost of privacy over epsilon

MSE.plot <- MSE.plot.mat %>%
  ggplot(aes(x = epsilon, y = mean)) +
  geom_ribbon(aes(ymin = mean - down, ymax = mean + up, group = epsilon[1]),
              fill = "thistle2") +
  geom_line(aes(x = epsilon, y = mean,
                group = epsilon[1]), col = "mediumpurple3", size = 1) +
  geom_point(aes(x = epsilon, y = mean), col = "purple4", size = 2) + 
  labs(title = TeX(r'(Cost of Privacy)'),
       x = TeX(r'($\epsilon$)'),
       y = TeX(r'(Cost)'), color = "Legend Title\n") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.25,0.7),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

# saving the plot

ggsave(plot = MSE.plot, filename = "MSEvsEpsilonBetahatPNP.pdf", device = "pdf",
       width = 4, height = 3)