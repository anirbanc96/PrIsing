################################################################################
############################# Required Packages ################################

require(tidyverse)
require(ggraph)

source("Functions.R")

################################################################################
############################## Loading Data ####################################

# Loading Data from study 22140-0001

data1 <- read_tsv("22140-0001-Data.tsv")

# Loading Data from study 22140-0002

data2 <- read_tsv("22140-0002-Data.tsv")

# Loading Data from study 22140-0003

data3 <- read_tsv("22140-0003-Data.tsv")
################################################################################
############################# Data Wrangling ###################################

############################# Social Network ###################################

data <- data2 %>%
  filter(HIV1 %in% c(0,1), HIV2 %in% c(0,1)) %>%
  
  # TIETYPE = 1 <- Social Network
  # TIETYPE = 3 <- Sexual Network
  
  # STUDYNUM = 1 <- Colorado Springs Study
  # DYADTYPE = 0 <- Ego-Alter Dyads
  
  filter(STUDYNUM == 1, TIETYPE == 1, DYADTYPE == 0) %>%
  select(ID1, ID2, HIV1, HIV2)

# Building dataframe with respondent id and HIV status

HIV.status.df <- data.frame(ID = c(data$ID1, data$ID2),
                            HIV = c(data$HIV1, data$HIV2))
HIV.status.df <- HIV.status.df[!duplicated(HIV.status.df),]

# Constructing dataframe having edges as vertex pairs

data.edges <- graph.data.frame(d = data[, 1:2], directed = F)

# Constructing adjacency matrix from edges

A.data <- as_adjacency_matrix(data.edges) %>%
  as.matrix() 

# Removing duplicate edges

A.data[A.data>1] = 1

# Number of vertices

N <- length(HIV.status.df$HIV)

# build the graph object

network <- graph_from_adjacency_matrix(A.data, mode = "undirected")
network <- simplify(network, remove.multiple = T, remove.loops = T) 

# HIV status of individuals

HIV.status <- rep(0,N)

for (i in 1:N){
  
  if (HIV.status.df[i,2] == 1){
    HIV.status[i] <- "Positive"
  }
  else{
    HIV.status[i] <- "Negative"
  }
  
}

# Drawing the social network

social.network <- ggraph(network, layout = "graphopt") +
  geom_edge_link(color = "lightgray", width = 0.3) +
  geom_node_point(aes(color = HIV.status), size = 1.2) + 
  scale_colour_manual("HIV Status", breaks = c("Positive", "Negative"),
                      values = c("firebrick2", "limegreen")) +
  theme_bw(base_size = 18) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.key = element_rect(fill = "white", colour = "black"),
        # legend.text = element_text(size = 14),
        legend.position = "bottom",
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

ggsave(plot = social.network, filename = "HIVSocialNetwork.pdf", device = "pdf",
       width = 6, height = 5)

################################################################################