Experiment: MSE vs Epsilon Plot for beta around non-private beta hat

Parameters: Betahat <- Estimated from HIV status in the data
	    A.N <- Adjacency Matrix from Data Network
	    J.N <- A.N/(rowSums(A.N) %o% rowSums(A.N)
	    Privacy Parameters: Epsilon = c(5, 7.5, 10, 20, 50, 100, InF), Delta = 1/N
	    BetaRange <- Betahat +- 0.5

Simulation Details: Number of repititions done: 500