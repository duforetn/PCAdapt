file = "results4pops_K3_fast"

Factors <- as.matrix(read.table(paste(file, '.scores', sep = '')))
K = nrow(Factors)
nFiles = 2
nsnp = 10000
par(mfrow = c(K, 1))
for (k in 1:K) plot(Factors[k, ], xlab = 'individuals', ylab = paste("PC", k))

Scores <- NULL
for (f in 1:nFiles){

Scores <- rbind(Scores, as.matrix(read.table(paste(file, '_file', f, sep = ""), header = TRUE)))

}

assign <- Scores[, 2]
top5p <- quantile(Scores[, 1], probs = .95)
x11(); plot(Scores[, 1], col = assign, ylab = "h' statistic", pch = 19, cex = .4)
abline(top5p, 0, lty = 2)
legend('topr', pch = c(19, 19, 19, 45), col = c(1:3, 1), legend = c('PC 1', 'PC 2', 'PC 3', 'top 5%'))

