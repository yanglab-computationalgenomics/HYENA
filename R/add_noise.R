#!/bin/env

oldw <- getOption("warn")
options(warn = -1)


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MatrixGenerics))
suppressPackageStartupMessages(library(RMThreshold))

option_list <- list( 
  make_option(c("-e", "--exp"), action="store",
              help="Gene expression input (output file from quantNorm.R)"),
  make_option(c("-w", "--write"), action="store", default="./intermediate/",
              help="Folder for intermediate files [default %default]"),
  make_option(c("-x", "--prefix"), action="store", default="DAT",
              help="Prefix for output files [default %default]"),
  make_option(c("-v","--verbose"), action="store_true", default=FALSE,
              help="Verbose mode for troubleshooting [default %default]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="")
args <- parse_args(parser, positional_arguments=0)
opt <- args$options

# load expression data
dat.expression <- fread(opt$exp, stringsAsFactors = FALSE)
dat.expression <- as.matrix(dat.expression)
print("Expression data loaded:")
if(opt$verbose) {
  print(dat.expression[1:5,1:5])
  }

# annotate rows
rownames(dat.expression) <- dat.expression[,1]
dat.expression <- dat.expression[,-1]
class(dat.expression) <- "numeric"

# remove genes that are 0 for all patients
dat.expression <- dat.expression[rowSums(dat.expression[]) > 0 ,]

if(opt$verbose) {
  print(dat.expression[1:5,1:5])
  }


# add noise to the entire expression matrix
print("Adding Gaussian noise...")
dat.expression.noisy <- add.Gaussian.noise(as.matrix(dat.expression),
                                           mean = 0.000000001,
                                           stddev = 0.000000001,
                                           symm = FALSE)
dat.expression.noisy <- as.data.frame(dat.expression.noisy)

if(opt$verbose) {
  print(dat.expression.noisy[1:5,1:5])
  }


if (file.exists(opt$write)) {
 write.table(dat.expression.noisy, file = paste0(opt$write, opt$prefix, ".exp_quant_noisy.txt"), quote = FALSE, sep = "\t")
} else {
 dir.create(opt$write)
 write.table(dat.expression.noisy, file = paste0(opt$write, opt$prefix, ".exp_quant_noisy.txt"), quote = FALSE, sep = "\t")
}


fpkm.qn <- dat.expression.noisy

# transpose the expression matrix so that the variables (genes) are in columns
fpkm.qn.t <- t(fpkm.qn)
if(opt$verbose){
  print("Expression matrix transposed.")
  print(fpkm.qn.t[1:5,1:5])
}

# prepare RLE values
fpkm.qn.log <- log(fpkm.qn + 1)
fpkm.qn.log <- as.matrix(fpkm.qn.log)
features_meds <- rowMedians(fpkm.qn.log)
med_devs <- fpkm.qn.log - features_meds
med_devs <- t(med_devs) # necessary for downstream analysis (plotting)

# Transform expression data to normal scores
print("Transforming expression data to normal scores...")

ns_transform <- function (y) {
  y <- qqnorm(y, plot.it = FALSE)$x
}

temp <- rownames(fpkm.qn.t)
fpkm.qn.t.ns <- apply(fpkm.qn.t, 2, ns_transform)
rownames(fpkm.qn.t.ns) <- temp
#if(opt$verbose){print(head(fpkm.qn.t.ns))}

# Perform  PCA using prcomp() ####
pca.prcomp <- prcomp(fpkm.qn.t.ns, center = FALSE, scale. = FALSE)
var_explained <- (pca.prcomp$sdev^2/sum(pca.prcomp$sdev^2))*100
print("PCA analysis completed.")

if (!(file.exists(paste0(opt$write, "plots")))) {dir.create(paste0(opt$write, "plots"), recursive = TRUE)}
pdf(paste0(opt$write, "plots/", opt$prefix, ".var_explained.pdf"))
barplot(var_explained[1:10], ylab = "% var explained", names = colnames(pca.prcomp$x)[1:10])
invisible(dev.off())

pdf(paste0(opt$write, "plots/", opt$prefix, ".pca.pdf"))
plot(pca.prcomp$x[,1:2])
invisible(dev.off())

save(fpkm.qn.t,fpkm.qn.t.ns,pca.prcomp,med_devs, file = paste0(opt$write, opt$prefix, ".exp.Rdata"))

#write.table(as.data.frame(pca.prcomp$x), file = paste0(opt$write, opt$prefix, ".prcomp_centered_scaled_x.txt"), quote = FALSE, sep = "\t")
write.table(pca.prcomp$sdev, file = paste0(opt$write, opt$prefix, ".prcomp_centered_scaled_sdev.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#write.table(pca.prcomp$center, file = paste0(opt$write, opt$prefix, ".prcomp_centered_scaled_center.txt"), quote = FALSE, sep = "\t", col.names = FALSE)
#write.table(pca.prcomp$scale, file = paste0(opt$write, opt$prefix, ".prcomp_centered_scaled_scale.txt"), quote = FALSE, sep = "\t", col.names = FALSE)
#write.table(as.data.frame(pca.prcomp$rotation), file = paste0(opt$write, opt$prefix, ".prcomp_centered_scaled_rotation.txt"), quote = FALSE, sep = "\t")
write.table(summary(pca.prcomp)$importance, file = paste0(opt$write, opt$prefix, ".prcomp_centered_scaled_importance.txt"), quote = FALSE, sep = "\t")



options(warn = oldw)