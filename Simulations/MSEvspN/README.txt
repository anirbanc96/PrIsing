Experiment: MSE vs alpha of MPLE estimate under Private and Non-Private Setup

Parameters: True Beta = c(0.5, 1.5) #Reflecting the phase transition at 1
	    N <- 2000
	    A.N <- Erdos Renyi graph on N vertices generated with edge probability p_N = N^{-alpha} with alpha = seq(0.1, 0.9, length.out = 9)
	    J.N <- A.N/Np_N
	    Privacy Parameters: Epsilon = c(5,7.5, 10), Delta = 1/N

Simulation Details: Number of repititions done: 500