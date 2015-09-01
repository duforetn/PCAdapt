Factors_list <- NULL
nrep = 5
nK = 3

for (rep in 1:nrep){
	Factors_list <- append(Factors_list, list(t(as.matrix(read.table(paste("results4pops_K", nK, "_", rep, '.scores', sep = ''))))))
}

cor_runs <- matrix(0, (nrep/2)*(nrep - 1), 3)
names_cor <- rep(0, (nrep/2)*(nrep - 1))
line = 0

for (i in 1:(nrep - 1)){

  	for (j in (i + 1):nrep){
		line = line + 1
	    	cor_ij <- cor(Factors_list[[i]], Factors_list[[j]])
		cor_runs[line, ] = apply(cor_ij**2, 1, max)
		names_cor[line] = paste(i, j)
  	}
}
rownames(cor_runs) <- names_cor

cor_runs
