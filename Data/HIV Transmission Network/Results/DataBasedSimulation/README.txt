Experiment: Data Based simulation using HIV Transmission network.
	    Errorbar plots of Private and Non-Private MPLE over a range of beta.

Parameters: Beta = seq(0, 2, by = 0.1)
	    A.N <- Adjacency Matrix from Data Network
	    J.N <- A.N/(rowSums(A.N) %o% rowSums(A.N)
	    Privacy Parameters: Epsilon = 5, Delta = 1/N

Simulation Details: Number of repititions done: 500