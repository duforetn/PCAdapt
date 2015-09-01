file = "results4pops_K"

nK = 6
system("rm tmp")
for (k in 1:nK){
	system(paste("head -n 1 ", file, k, ".stats >> tmp", sep = ""))
}

err <- read.table("tmp")
plot(err[, 2], pch = 18, xlab = 'K', ylab = 'errors')
