################################################################################
############################# Required Packages ################################

require(tidyverse)
require(ggraph)

source("Functions.R")

################################################################################
############################## Loading Data ####################################

source("DataManage.R")

################################################################################
############################# Data Wrangling ###################################

# build the graph object

del <- which(rowSums(A.data) == 0)
A.data <- A.data[-del,-del]

network <- graph_from_adjacency_matrix(A.data, mode = "undirected")
network <- simplify(network, remove.multiple = T, remove.loops = T) 

# opinion status of individuals
# from documentation opinion == 1 if conservative and 0 is liberal

opinion <- opinion[-del]

opinion.status <- rep(0,length(opinion))

for (i in 1:length(opinion)){
  
  if (opinion[i] == 1){
    opinion.status[i] <- "Conservative"
  }
  else{
    opinion.status[i] <- "Liberal"
  }
  
}

# Drawing the political blogs networks

pol.network <- ggraph(network, layout = "graphopt") +
  geom_edge_link(color = "lightgray", width = 0.1) +
  geom_node_point(aes(color = opinion.status), size = 0.7) + 
  scale_colour_manual("Political Leaning", breaks = c("Conservative", "Liberal"),
                      values = c("coral2", "cornflowerblue")) +
  theme_bw(base_size = 18) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.position = "bottom",
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

ggsave(plot = pol.network, filename = "PolBlogsNetwork.pdf", device = "pdf",
       width = 6, height = 5)

################################################################################