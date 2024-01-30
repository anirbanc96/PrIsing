require(tidyverse)
require(latex2exp)

MSE.Erdos.Renyi <- read.csv("MSEErdosRenyi.csv")[,-1]
MSE.HIV <- read.csv("MSEHIV.csv")[,-1]
MSE.PolBlogs <- read.csv("MSEPolBlogs.csv")[,-1]

MSE <- cbind(MSE.Erdos.Renyi, MSE.HIV[,-1], MSE.PolBlogs[,-1])
colnames(MSE) <- c("Epsilon", "ERL1", "ERG1", "HIVL1", "HIVG1", "POLL1", "POLG1")
MSE.plot <- MSE %>%
  ggplot(aes(x = as.factor(Epsilon), y = ERL1)) + 
  geom_line(aes(x = as.factor(Epsilon), y = ERL1, col = "ERL1", group = 1)) + 
  geom_point(aes(x = as.factor(Epsilon), y = ERL1), col = "navy") +
  geom_line(aes(x = as.factor(Epsilon), y = ERG1, col = "ERG1", group = 2)) +
  geom_point(aes(x = as.factor(Epsilon), y = ERG1), col = "deepskyblue") +
  geom_line(aes(x = as.factor(Epsilon), y = HIVL1, col = "HIVL1", group = 3)) +
  geom_point(aes(x = as.factor(Epsilon), y = HIVL1), col = "darkred") +
  geom_line(aes(x = as.factor(Epsilon), y = HIVG1, col = "HIVG1", group = 4)) +
  geom_point(aes(x = as.factor(Epsilon), y = HIVG1), col = "deeppink3") +
  geom_line(aes(x = as.factor(Epsilon), y = POLL1, col = "POLL1", group = 5)) +
  geom_point(aes(x = as.factor(Epsilon), y = POLL1), col = "darkolivegreen") +
  geom_line(aes(x = as.factor(Epsilon), y = POLG1, col = "POLG1", group = 6)) +
  geom_point(aes(x = as.factor(Epsilon), y = POLG1), col = "forestgreen") +
  scale_colour_manual(breaks = c("ERL1", "ERG1", "HIVL1", "HIVG1", "POLL1", "POLG1"),
                      values = c("cornflowerblue", "skyblue1", "indianred", "plum1",
                                 "palegreen4", "seagreen2"),
                      labels = c(expression(paste("ER (",beta, " = 0.5)")),
                                 expression(paste("ER (",beta, " = 1.5)")),
                                 expression(paste("HIV (",beta, " = 0.5)")),
                                 expression(paste("HIV (",beta, " = 1.5)")),
                                 expression(paste("POL (",beta, " = 0.5)")),
                                 expression(paste("POL (",beta, " = 1.5)")))) +
  labs(title = "MSE of MPLE", x = TeX(r'($\epsilon$)'),
       y = TeX(r'(MSE)'), color = "Legend Title\n") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        legend.position = c(0.82,0.62),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))


ggsave(plot = MSE.plot, filename = "AllinOne.pdf", device = "pdf",
       width = 4, height = 3)

