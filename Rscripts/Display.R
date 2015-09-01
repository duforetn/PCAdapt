file = "results4pops_K3"

Factors <- as.matrix(read.table(paste(file, '.scores', sep = '')))
K = nrow(Factors)
nsnp = 6000
par(mfrow = c(K, 1))
for (k in 1:K) plot(Factors[k, ], xlab = 'individuals')

system(paste("tail -n ", nsnp, file, "> BF"))
BF <- as.matrix(read.table('BF'))

assign <- apply(BF[, 3:(2 + K)], 1, which.max)
top5p <- quantile(BF[, 1], probs = .95)
x11(); plot(BF[, 1], col = assign, ylab = "log10(BF)", pch = 19, cex = .4)
abline(top5p, 0, lty = 2)
legend('topr', pch = c(19, 19, 19, 45), col = c(1:3, 1), legend = c('Factor 1', 'Factor 2', 'Factor 3', 'top 5%'))

