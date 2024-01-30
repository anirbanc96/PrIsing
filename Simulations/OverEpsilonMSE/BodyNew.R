source("Functions.R")
require(ggpubr)

MSE.dat.f <- function(N, J.N, beta1, beta2, epsilon.vec, delta,
                       priv = T, n.rep){
  
  delta <- 1/N
  
  MSE.p.le1 <- rep(0, length(epsilon.vec))
  MSE.flip.le1 <- rep(0, length(epsilon.vec))
  MSE.p.ge1 <- rep(0, length(epsilon.vec))
  MSE.flip.ge1 <- rep(0, length(epsilon.vec))
  
  for (i in 1:length(epsilon.vec)){
    
    beta.est.iter.priv <- sim.beta.est(n.rep, beta1, N, J.N,
                                       priv = T, epsilon.vec[i], delta)
    
    MSE.p.le1[i] <- sum((beta.est.iter.priv - beta1)^2)/n.rep
    
    beta.est.iter.priv <- sim.beta.est(n.rep, beta2, N, J.N,
                                       priv = T, epsilon.vec[i], delta)
    
    MSE.p.ge1[i] <- sum((beta.est.iter.priv - beta2)^2)/n.rep
    
    print (i)
  }
  
  MSE.dat <- cbind(epsilon.vec, MSE.p.le1, MSE.p.ge1)
  colnames(MSE.dat) <- c("Epsilon", "p.le1", "p.ge1")
  
  MSE.dat <- as_tibble(MSE.dat)
  MSE.dat[,1] <- as.factor(epsilon.vec)
  
  return (MSE.dat)
  
}

MSE.plot.f <- function(MSE.dat, name){
  
  MSE.plot <- MSE.dat %>%
    ggplot(aes(x = Epsilon, y = p.le1)) +
    geom_point(aes(x = Epsilon, y = p.le1), col = "navy") +
    geom_line(aes(x = Epsilon, y = p.le1, col = "PBL1", group = Epsilon[1])) +
    geom_point(aes(x = Epsilon, y = p.ge1), col = "darkorange") +
    geom_line(aes(x = Epsilon, y = p.ge1, col = "PBG1", group = Epsilon[2])) +
    scale_colour_manual(breaks = c("PBL1", "PBG1"),
                        values = c("cornflowerblue", "orange"),
                        labels = c(expression(paste(beta, " = ", 0.5)),
                                   expression(paste(beta, " = ", 1.5)))) +
    labs(title = name, x = TeX(r'($\epsilon$)'),
         y = TeX(r'(MSE)'), color = "Legend Title\n") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title=element_blank(),
          legend.position = c(0.75,0.65),
          legend.text.align = 0,
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"))
  
  MSE.plot
  
  return (MSE.plot)
  
}

epsilon.vec <- c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 25, 100)

beta1 <- 0.5; beta2 <- 1.5

n.rep <- 500

# Erdos-Renyi MSE

N <- 2000
p <- N^{-1/3}

# Matrix in the Hamiltonian
J.N <- as_adjacency_matrix(erdos.renyi.game(N, p, type="gnp"), sparse = FALSE)

# Appropriate normalisation
J.N <- J.N/(N*p)

MSE.Erdos.Renyi <- MSE.dat.f(N, J.N, beta1, beta2, epsilon.vec, delta,
                              priv = T, n.rep)
MSE.Erdos.Renyi.plot <- MSE.plot.f(MSE.Erdos.Renyi,
                                   name = "Erdos-Renyi Graph")

# HIV MSE

source("DataHIV.R")

MSE.HIV <- MSE.dat.f(N, J.N, beta1, beta2, epsilon.vec, delta,
                      priv = T, n.rep)
MSE.HIV.plot <- MSE.plot.f(MSE.HIV,
                           name = "HIV Network")

# PolBlogs MSE

source("DataPolBlogs.R")

MSE.PolBlogs <- MSE.dat.f(N, J.N, beta1, beta2, epsilon.vec, delta,
                      priv = T, n.rep)
MSE.PolBlogs.plot <- MSE.plot.f(MSE.PolBlogs,
                                name = "PolBlog Network")

MSE.Erdos.Renyi.plot
MSE.HIV.plot
MSE.PolBlogs.plot

MSE.plot <- ggarrange(MSE.Erdos.Renyi.plot, MSE.HIV.plot, MSE.PolBlogs.plot, ncol = 3,
          common.legend = T, legend = "right")

write.csv(MSE.Erdos.Renyi, "MSEErdosRenyi.csv")
write.csv(MSE.HIV, "MSEHIV.csv")
write.csv(MSE.PolBlogs, "MSEPolBlogs.csv")

ggsave(plot = MSE.plot, filename = "MPLEvsEpsilon.pdf", device = "pdf",
       width = 7, height = 4)

source("AllinOne.R")