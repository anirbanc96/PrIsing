################################################################################
############################# Required Packages ################################

library(Matrix)

library(igraph)

library(gsbm)

require(tidyverse)

source("Functions.R")

################################################################################
############################## Loading Data ####################################

#Loading the data "blogosphere"

data(blogosphere)

# Getting the Adjacency Matrix

A.data <- blogosphere$A

# Political Leaning of the Blogs

opinion <- blogosphere$opinion

################################################################################
############################# Data Wrangling ###################################

# Number of vertices

bigdeg = which(rowSums(A.data) > 50)

A.data = A.data[-bigdeg,-bigdeg]
opinion = opinion[-bigdeg]

degmax1 = sapply(rowSums(A.data), function(x) max(x,1))
row.ij <- degmax1%o%degmax1

N <- length(opinion)

# Hamiltonian for the Ising Model

J.N <- A.data/sqrt(row.ij)
################################################################################

liberal.large <- which(opinion[bigdeg] == 0)
conservative.large <- which(opinion[bigdeg] == 1)

write.csv(blogosphere$names[liberal.large], "Large Liberals.csv")
write.csv(blogosphere$names[conservative.large], "Large Conservatives.csv")

################################################################################