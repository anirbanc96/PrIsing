################################################################################
############################# Required Packages ################################

require(tidyverse)

# Load required functions

source("Functions.R")

# Load data and Hamiltonian of Ising Model

source("DataManage.R")

################################################################################
############################# Setting Parameters ###############################
# value of privacy parameters

epsilon <- 5
delta <- 1/(N)

# range of beta

beta.range <- seq(0,3.2, by = 0.2)

# number of experiments

n.rep <- 500

################################################################################
########################### Running Experiments ################################

beta.len <- length(beta.range)
beta.est.nonpriv <- matrix(0, nrow = beta.len, ncol = n.rep)
beta.est.priv <- matrix(0, nrow = beta.len, ncol = n.rep)

# Running Experiments

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

################################################################################

# Run the following lines for re-plotting from old saved data.

# beta.est.nonpriv <- read.csv("Non-Private.csv")[,-1]
# beta.est.priv <- read.csv("Private.csv")[,-1]
# 
# beta.range <- as.vector(beta.est.priv$beta.range)
# beta.len <- length(beta.range)

################################################################################

# Plotting true beta vs MPLE estimate for private and non-private estimates

betaest.np <- matrix(0, nrow = beta.len, ncol = 4)
colnames(betaest.np) <- c("beta", "mean", "up", "down")
betaest.np <- as.data.frame(betaest.np)

# Making data frame of non-private estimates

betaest.np <- betaest.np %>%
  mutate(beta = beta.range) %>%
  mutate(mean = apply(beta.est.nonpriv[,-c(1:2)], 1, mean)) %>%
  mutate(up = mean + apply(beta.est.nonpriv[,-c(1:2)], 1, sd)) %>%
  mutate(down = mean - apply(beta.est.nonpriv[,-c(1:2)], 1, sd))

# Making data frame of private estimates

betaest.p <- matrix(0, nrow = beta.len, ncol = 4)
colnames(betaest.p) <- c("beta", "mean", "up", "down")
betaest.p <- as.data.frame(betaest.p)
betaest.p <- betaest.p %>%
  mutate(beta = beta.range) %>%
  mutate(mean = apply(beta.est.priv[,-c(1:2)], 1, mean)) %>%
  mutate(up = mean + apply(beta.est.priv[,-c(1:2)], 1, sd)) %>%
  mutate(down = mean - apply(beta.est.priv[,-c(1:2)], 1, sd))

# Plotting comparison of non-private and private MPLE

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
  labs(title = TeX(r'(Private and Non-Private MPLE)'), x = TeX(r'($\beta$)'),
       y = TeX(r'(MPLE)'), color = "Legend Title\n") +
  theme_bw() +
  scale_colour_manual("Privacy", breaks = c("Non-Private", "Private"),
                      values = c("cornflowerblue", "orange"),
                      labels = c(expression(paste(epsilon, " = ", infinity)),
                                 expression(paste(epsilon, " = ", 5)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.25,0.7),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

write.csv(beta.est.nonpriv, "Non-Private.csv")
write.csv(beta.est.priv, "Private.csv")

# saving the plot

ggsave(plot = est.plot, filename = "PolBlogsSocial.pdf", device = "pdf",
        width = 4, height = 3)

################################################################################