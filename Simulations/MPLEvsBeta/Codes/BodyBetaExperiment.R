################################################################################

# loading required functions and packages
source("Functions.R")

################################################################################

# Input Parameters

# Size of underlying graph
N = 2000

# probability of erdos renyi graph
p = N^{-1/3}

# value of privacy parameters
epsilon <- 5
delta <- 1/N

# values of true beta
beta.range <- (1:20)/10

# number of experiments done
n.rep <- 500

# Matrix in the Hamiltonian
J.N <- as_adjacency_matrix(erdos.renyi.game(N, p, type="gnp"), sparse = FALSE)

# Appropriate normalisation
J.N = J.N/(N*p)


beta.len <- length(beta.range)

# matrix of non-private beta estimates for different true beta values
beta.est.nonpriv <- matrix(0, nrow = beta.len, ncol = n.rep)

# matrix of non-private beta estimates for different true beta values
beta.est.priv <- matrix(0, nrow = beta.len, ncol = n.rep)

# Running experiments for MPLE of beta
start <- Sys.time()
for (i in 1:beta.len){
  
  beta.est.iter <- sim.beta.est(n.rep, beta.range[i], N, J.N,
                                priv = T, epsilon, delta)
  beta.est.nonpriv[i, ] <- beta.est.iter[, 1]
  beta.est.priv[i, ] <- beta.est.iter[, 2]
  print (i)
}
Sys.time() - start

# non-private and private estimates of beta
# each row the dataframe contains n.rep many estimates of beta given the
# underlying graph J_N.

beta.est.nonpriv <- as.data.frame(cbind(beta.range, beta.est.nonpriv))
beta.est.priv <- as.data.frame(cbind(beta.range, beta.est.priv))

# saving estimated private and non-private beta values

write.csv(beta.est.nonpriv, file = "nonprivate.csv")
write.csv(beta.est.priv, file = "private.csv")

################################################################################

# Run the following lines for re-plotting from old saved data.

# beta.est.nonpriv <- read.csv("nonprivate.csv")
# beta.est.priv <- read.csv("private.csv")
# beta.range <- beta.est.nonpriv[,2]
# beta.len <- length(beta.range)

################################################################################

# dataframe of non-private MPLE with errobar bounds
betaest.np <- matrix(0, nrow = beta.len, ncol = 4)
colnames(betaest.np) <- c("beta", "mean", "up", "down")
betaest.np <- as.data.frame(betaest.np)
betaest.np <- betaest.np %>%
  mutate(beta = beta.range) %>%
  mutate(mean = apply(beta.est.nonpriv[,-c(1:2)], 1, mean)) %>%
  mutate(up = mean + apply(beta.est.nonpriv[,-c(1:2)], 1, sd)) %>%
  mutate(down = mean - apply(beta.est.nonpriv[,-c(1:2)], 1, sd))

# dataframe of private MPLE with errobar bounds
betaest.p <- matrix(0, nrow = beta.len, ncol = 4)
colnames(betaest.p) <- c("beta", "mean", "up", "down")
betaest.p <- as.data.frame(betaest.p)
betaest.p <- betaest.p %>%
  mutate(beta = beta.range) %>%
  mutate(mean = apply(beta.est.priv[,-c(1:2)], 1, mean)) %>%
  mutate(up = mean + apply(beta.est.priv[,-c(1:2)], 1, sd)) %>%
  mutate(down = mean - apply(beta.est.priv[,-c(1:2)], 1, sd))

# plotting MPLE with errorbar for private and non-private estimator

est.plot <- betaest.np %>%
  ggplot(aes(x = beta, y = beta)) +
  geom_line(aes(x = beta, y = beta), lty = "dashed") +
  geom_point(aes(x = beta, y = betaest.p$mean), color = "darkorange2",
             size = 2) +
  geom_errorbar(data = betaest.p, aes(x = beta, ymin = down, ymax = up,
                                      color = "Private"), size = 1) +
  geom_point(aes(x = beta, y = mean), color = "navy", size = 2) +
  geom_errorbar(aes(x = beta, ymin = down, ymax = up, color = "Non-Private"),
                size = 1) +
  labs(title = TeX(r'(MPLE and Error Bars)'), x = TeX(r'($\beta$)'),
       y = TeX(r'(MPLE)'), color = "Legend Title\n") +
  theme_bw() +
  scale_colour_manual("Privacy", breaks = c("Non-Private", "Private"),
                      values = c("cornflowerblue", "orange"),
                      labels = c(expression(paste(epsilon, " = ", infinity)),
                                 expression(paste(epsilon, " = ", 5)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.77,0.25),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

# saving the plot
ggsave(plot = est.plot, filename = "MPLEvsBeta.pdf", device = "pdf",
       width = 4, height = 3)

################################################################################