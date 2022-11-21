oldw <- getOption("warn")
options(warn = -1)


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(EDASeq))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))

option_list <- list( 
  make_option(c("-e", "--exp"), action="store",
              help="Gene expression matrix to be used as the input"),
  make_option(c("-w", "--write"), action="store", default="./intermediate/",
              help="Folder for intermediate files [default %default]"),
  make_option(c("-x", "--prefix"), action="store", default="DAT",
              help="Prefix for output files [default %default]"),
  make_option(c("-i", "--id"), action="store",
              help="Sample id file")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="")
args <- parse_args(parser, positional_arguments=0)
opt <- args$options

#####################################################################
#####################################################################

samples <- fread(opt$id, stringsAsFactors = FALSE)
aliqid_rnaseq <- samples$aliquot_id_rnaseq

exp <- fread(opt$exp, stringsAsFactors = FALSE)
class(exp) <- "data.frame"
header <- c(colnames(exp)[1],aliqid_rnaseq)
exp <- exp[,header]

sample_num <- ncol(exp)-1

# prep the matrix
exp <- as.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,-1]
class(exp) <- "numeric"


# plot RLE BEFORE quantile normalization

if (!(file.exists(opt$write))) {
  dir.create(opt$write, recursive = TRUE)
}


pdf(paste0(opt$write, opt$prefix, ".RLE_before.pdf"))
plotRLE(as.matrix(exp), outline=FALSE, xaxt = "n", xlab = paste0("samples(",sample_num,")"), ylab = "relative log expression", main = "RLE before quantile normalization")
invisible(dev.off())


# perform quantile normalization
print("Quantile normalizing expression data...")
exp.quant <- normalize.quantiles(exp)
rownames(exp.quant) <- rownames(exp)
colnames(exp.quant) <- colnames(exp)


# plot RLE AFTER quantile normalization
pdf(paste0(opt$write, opt$prefix, ".RLE_after.pdf"))
plotRLE(as.matrix(exp.quant), outline=FALSE, xaxt = "n", xlab = paste0("samples(",sample_num,")"), ylab = "relative log expression", main = "RLE after quantile normalization")
invisible(dev.off())


write.table(exp.quant, file = paste0(opt$write, opt$prefix, ".exp_quant.txt"), quote = FALSE, sep = "\t")




options(warn = oldw)