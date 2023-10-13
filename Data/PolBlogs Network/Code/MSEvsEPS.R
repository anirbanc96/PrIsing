################################################################################
############################# Required Packages ################################

require(latex2exp)
require(tidyverse)

# Load required functions
source("Functions.R")

# Load data and Hamiltonian of Ising Model
source("DataManage.R")
################################################################################

# MPLE estimate of beta
beta.hat.nonpriv <- beta.hat.N((2*opinion-1), N, J.N)

# Selecting range of beta around beta.hat.nonpriv
beta.range <- seq(beta.hat.nonpriv-0.5, beta.hat.nonpriv+0.5, length.out = 5)

# value of privacy parameters
epsilon.vec <- c(5, 7.5, 10, 20, 50, 100)

# number of experiments done
n.rep <- 500

# matrix of MSE values
MSE.mat <- matrix(0, length(beta.range), (length(epsilon.vec) + 1))

# running a loop to calculate MSE values over beta and epsilon

# MSE values for private estimator
for (i in 1:length(epsilon.vec)){
  
  
    epsilon <- epsilon.vec[i]
    delta <- 1/N
  
    for (j in 1:length(beta.range)){
    
      beta.est.iter <- sim.beta.est(n.rep, beta.range[j], N, J.N,
                                    priv = T, epsilon, delta)[,2]
    
      MSE.mat[j,i] <- sum((beta.est.iter - beta.range[j])^2)/n.rep
    
      print (c(j,i))
    
    }
}

# MSE values for non-private estimator
for (j in 1:length(beta.range)){
  
  beta.est.iter <- sim.beta.est(n.rep, beta.range[j], N, J.N,
                                priv = F)
  
  MSE.mat[j,(length(epsilon.vec) + 1)] <- sum((beta.est.iter - beta.range[j])^2)/n.rep
  
  print (c(j,(length(epsilon.vec) + 1)))
  
}

# Storing MSE values in a dataframe
MSE.dat <- cbind(rep(0, (length(epsilon.vec) + 1)),
                 t(MSE.mat))



colnames(MSE.dat) <- c("epsilon", "beta1", "beta2", "beta3", "beta4", "beta5")
MSE.dat <- as_tibble(MSE.dat)

MSE.dat[,1] <- as.factor(c(epsilon.vec,Inf))

write.csv(MSE.dat, file = "MSEvsEpsilonValuesPol.csv")

################################################################################

# Run the following lines for re-plotting from old saved data.

# MSE.dat <- read.csv("MSEvsEpsilonValuesPol.csv")[,-1]
# MSE.dat[,1] <- as.factor(MSE.dat[,1])

################################################################################

# Plotting MSE vs epsilon for multiple beta values
MSE.plot <- MSE.dat %>%
  ggplot(aes(x = epsilon, y = beta1)) +
  geom_point(aes(x = epsilon, y = beta1), col = "navy") + 
  geom_line(aes(x = epsilon, y = beta1, col = "Beta 1",
                group = epsilon[1])) +
  geom_point(aes(x = epsilon, y = beta2), col = "darkorange") + 
  geom_line(aes(x = epsilon, y = beta2, col = "Beta 2",
                group = epsilon[2])) +
  geom_point(aes(x = epsilon, y = beta3), col = "darkred") + 
  geom_line(aes(x = epsilon, y = beta3, col = "Beta 3",
                group = epsilon[3])) +
  geom_point(aes(x = epsilon, y = beta4), col = "darkmagenta") + 
  geom_line(aes(x = epsilon, y = beta4, col = "Beta 4",
                group = epsilon[4])) +
  geom_point(aes(x = epsilon, y = beta5), col = "forestgreen") + 
  geom_line(aes(x = epsilon, y = beta5, col = "Beta 5",
                group = epsilon[5])) +
  scale_colour_manual("Beta", breaks = c("Beta 1",
                                            "Beta 2", "Beta 3",
                                            "Beta 4", "Beta 5"),
                      values = c("cornflowerblue", "orange", "red", "plum",
                                 "darkseagreen"),
                      labels = c(expression(paste(beta, " = ", hat(beta), " - ", 0.5)),
                                 expression(paste(beta, " = ", hat(beta), " - ", 0.25)),
                                 expression(paste(beta, " = ", hat(beta))),
                                 expression(paste(beta, " = ", hat(beta), " + ", 0.25)),
                                 expression(paste(beta, " = ", hat(beta), " + ", 0.5)))) +
  labs(title = TeX(r'(MSE vs $\epsilon$)'), x = TeX(r'($\epsilon$)'),
       y = TeX(r'(MSE)'), color = "Legend Title\n") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.6),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

# saving the plot
ggsave(plot = MSE.plot, filename = "MSEvsEpsilonPol.pdf", device = "pdf",
       width = 4, height = 3)

