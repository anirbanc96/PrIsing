Experiment: MSE vs N of MPLE estimate under Private and Non-Private Setup

Parameters: True Beta = c(0.5, 1.5) #Reflecting the phase transition at 1
	    N <- c(100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)
	    A.N <- Erdos Renyi graph on N vertices generated with edge probability p_N = N^{-1/3}
	    J.N <- A.N/Np_N
	    Privacy Parameters: Epsilon = c(5,7.5, 10), Delta = 1/N

Simulation Details: Number of repititions done: 500