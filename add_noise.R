library(data.table)
library(RMThreshold)

dat.expression <- fread("fpkm_uq_purityCorr_cnCorr.txt", stringsAsFactors = FALSE)
#dat.expression <- na.omit(dat.expression)
dat.expression.noisy <- add.Gaussian.noise(as.matrix(dat.expression[,3:ncol(dat.expression)]),
                                           mean = 0.000000001,
                                           stddev = 0.000000001,
                                           symm = FALSE)
dat.expression.noisy <- as.data.frame(dat.expression.noisy)
dat.expression.noisy <- cbind(dat.expression[,1:2], dat.expression.noisy)

write.table(dat.expression.noisy, "fpkm_uq_purityCorr_cnCorr_noisy.txt", sep = "\t", quote = FALSE, row.names = FALSE)
