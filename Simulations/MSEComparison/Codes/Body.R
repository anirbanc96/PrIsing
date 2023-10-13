################################################################################

# Loading required functions and packages
source("Functions.R")

################################################################################

# vector of sample sizes
N.vec <- c(100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)
N.len <- length(N.vec)

# high and low temperature regimes
beta.vec <- c(0.5,1.5)

# values of epsilon parameter
epsilon <- c(5, 7.5, 10)

# number of experiments done
n.rep <- 500

################################################################################

# HIGH TEMPERATURE REGIME

# matrix of MSE values for private MPLE
MSE.p.le1 <- matrix(0, length(epsilon), N.len)

# vector of MSE values for non-private MSE
MSE.np.le1 <- rep(0, N.len)

for (i in 1:N.len){
  
  N = N.vec[i]
  
  # delta parameter required for privacy
  delta <- 1/N
  
  # probability of erdos renyi graph
  p = N^{-1/3}
  
  # Matrix in the Hamiltonian
  J.N <- as_adjacency_matrix(erdos.renyi.game(N, p, type="gnp"), sparse = FALSE)
  
  # Appropriate normalisation
  J.N = J.N/(N*p)
  
  beta.est.iter <- sim.beta.est(n.rep, beta.vec[1], N, J.N,
                                priv = F)
  
  # non-private MSE for the current iteration
  MSE.np.le1[i] <- sum((beta.est.iter - beta.vec[1])^2)/n.rep
  
  # private MSE for each value of epsilon in the current iteration
  for (j in 1:length(epsilon)){
  
    beta.est.iter.priv <- sim.beta.est(n.rep, beta.vec[1], N, J.N,
                                  priv = T, epsilon[j], delta)
    MSE.p.le1[j,i] <- sum((beta.est.iter.priv - beta.vec[1])^2)/n.rep
  
  }
  
  print (i)
}

################################################################################

# LOW TEMPERATURE REGIME

# matrix of MSE values for private MPLE
MSE.p.ge1 <- matrix(0, length(epsilon), N.len)

# vector of MSE values for non-private MSE
MSE.np.ge1 <- rep(0, N.len)

for (i in 1:N.len){
  
  N = N.vec[i]
  
  # delta parameter required for privacy
  delta <- 1/N
  
  # probability of erdos renyi graph
  p = N^{-1/3}
  
  # Matrix in the Hamiltonian
  J.N <- as_adjacency_matrix(erdos.renyi.game(N, p, type="gnp"), sparse = FALSE)
  
  # Appropriate normalisation
  J.N = J.N/(N*p)
  
  beta.est.iter <- sim.beta.est(n.rep, beta.vec[2], N, J.N,
                                priv = F)
  
  # non-private MSE for the current iteration
  MSE.np.ge1[i] <- sum((beta.est.iter - beta.vec[2])^2)/n.rep
  
  # private MSE for each value of epsilon in the current iteration
  for (j in 1:length(epsilon)){
    
    beta.est.iter.priv <- sim.beta.est(n.rep, beta.vec[2], N, J.N,
                                       priv = T, epsilon[j], delta)
    MSE.p.ge1[j,i] <- sum((beta.est.iter.priv - beta.vec[2])^2)/n.rep
    
  }
  
  print (i)
}

################################################################################

# Concatinate all MSE values forming a dataframe of MSE for different values
# of N

MSE.dat <- cbind(N.vec, MSE.np.le1, t(MSE.p.le1), MSE.np.ge1, t(MSE.p.ge1))
colnames(MSE.dat) <- c("N", "np.le1", "p.le1.e1", "p.le1.e2", "p.le1.e3",
                       "np.ge1", "p.ge1.e1", "p.ge1.e2", "p.ge1.e3")
MSE.dat <- as_tibble(MSE.dat)

# Saving MSE dataframe
write.csv(MSE.dat, file = "MSEValues.csv")

################################################################################

# Run the following lines for re-plotting from old saved data.

# MSE.dat <- read.csv("MSEvalues.csv")
# N.vec <- MSE.dat[,1]

################################################################################

# Plotting MSE vs N over different epsilon values in high temperature 
# regime

MSE.le1.plot <- MSE.dat %>%
  ggplot(aes(x = N.vec, y = np.le1)) +
  geom_point(aes(x = N.vec, y = np.le1), col = "navy") + 
  geom_line(aes(x = N.vec, y = np.le1, col = "Non-Private")) +
  geom_point(aes(x = N.vec, y = p.le1.e1), col = "darkorange") + 
  geom_line(aes(x = N.vec, y = p.le1.e1, col = "Private 1")) +
  geom_point(aes(x = N.vec, y = p.le1.e2), col = "darkred") + 
  geom_line(aes(x = N.vec, y = p.le1.e2, col = "Private 2")) +
  geom_point(aes(x = N.vec, y = p.le1.e3), col = "darkmagenta") + 
  geom_line(aes(x = N.vec, y = p.le1.e3, col = "Private 3")) +
  scale_colour_manual("Privacy", breaks = c("Non-Private",
                                            "Private 1", "Private 2",
                                            "Private 3"),
                      values = c("cornflowerblue", "orange", "red", "plum"),
                      labels = c(expression(paste(epsilon, " = ", infinity)),
                                 expression(paste(epsilon, " = ", 5)),
                                 expression(paste(epsilon, " = ", 7.5)),
                                 expression(paste(epsilon, " = ", 10)))) +
  labs(title = TeX(r'(MSE of MPLE at $\beta=0.5$)'), x = TeX(r'($n$)'),
       y = TeX(r'(MSE)'), color = "Legend Title\n") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.65),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

# saving the plot
ggsave(plot = MSE.le1.plot, filename = "MSE05.pdf", device = "pdf",
       width = 4, height = 3)

################################################################################

# Plotting MSE vs N over different epsilon values in high temperature 
# regime

MSE.ge1.plot <- MSE.dat %>%
  ggplot(aes(x = N.vec, y = np.ge1)) +
  geom_point(aes(x = N.vec, y = np.ge1), col = "navy") + 
  geom_line(aes(x = N.vec, y = np.ge1, col = "Non-Private")) +
  geom_point(aes(x = N.vec, y = p.ge1.e1), col = "darkorange") + 
  geom_line(aes(x = N.vec, y = p.ge1.e1, col = "Private 1")) +
  geom_point(aes(x = N.vec, y = p.ge1.e2), col = "darkred") + 
  geom_line(aes(x = N.vec, y = p.ge1.e2, col = "Private 2")) +
  geom_point(aes(x = N.vec, y = p.ge1.e3), col = "darkmagenta") + 
  geom_line(aes(x = N.vec, y = p.ge1.e3, col = "Private 3")) +
  scale_colour_manual("Privacy", breaks = c("Non-Private",
                                            "Private 1", "Private 2",
                                            "Private 3"),
                      values = c("cornflowerblue", "orange", "red", "plum"),
                      labels = c(expression(paste(epsilon, " = ", infinity)),
                                 expression(paste(epsilon, " = ", 5)),
                                 expression(paste(epsilon, " = ", 7.5)),
                                 expression(paste(epsilon, " = ", 10)))) +
  labs(title = TeX(r'(MSE of MPLE at $\beta=1.5$)'), x = TeX(r'($n$)'),
       y = TeX(r'(MSE)'), color = "Legend Title\n") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,0.65),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

# saving the plot
ggsave(plot = MSE.ge1.plot, filename = "MSE15.pdf", device = "pdf",
       width = 4, height = 3)
