################################################################################
############################# Required Packages ################################

require(tidyverse)

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

# Sub-setting connections and corresponding HIV status

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

# Normalised Adjacency matrix

row.ij <- rowSums(A.data)%o%rowSums(A.data)

# Hamiltonian for the Ising Model

J.N <- A.data/sqrt(row.ij)
################################################################################